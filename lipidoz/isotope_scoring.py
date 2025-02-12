"""
lipidoz/isotope_scoring.py
Dylan Ross (dylan.ross@pnnl.gov)

    module for performing isotope scoring of aldehyde/criegee ions

    general method:
        1. use target RT to fit LC peak for precursor -> observed RT, peak width, intensity
        2. determine if the peak is saturated (either by some threshold or another method based on peak shape), if so:
            2.1 determine scan indices for leading edge of LC peak
            2.2 use new scan indices to extract mass spectrum (hopefully not saturated)
        3. perform peak finding on mass spectrum and match peaks to predicted isotopes
        4. score based on agreement between mass and abundance
        5. repeat steps 1-4 for all possible aldehyde and criegee fragments corresponding to all possible 
            double bond positions
"""


import io
from typing import Optional, Any, Dict, Tuple, List
from collections.abc import Generator
import threading

from matplotlib import pyplot as plt, rcParams, use as mpuse
import numpy as np
from scipy import spatial
from mzapy.isotopes import predict_m_m1_m2
from mzapy.peaks import lerp_1d, find_peaks_1d_localmax, find_peaks_1d_gauss, _gauss
from mzapy._util import _ppm_error

from lipidoz._util import (
    OzData,
    Formula,
    _polyunsat_ald_crg_formula, 
    _calc_dbp_bounds, 
    _debug_handler
)


_XIC_MIN_ABS_INTENSITY = 100
_MS1_MIN_ABS_INTENSITY = 100


# TODO: Type annotation for results, or better yet convert the results to a dataclass instead of nested dicts


def _plot_xic_with_peak(xic_rt, xic_int, prt, pht, pwt,
                        restrict_xlim=False, add_cmp_range=False):
    """
    plots XIC with fitted peak(s), returns png image as bytes

    Parameters
    ----------
    xic_rt : ``numpy.ndarray(float)``
        retention time component of XIC
    xic_int : ``numpy.ndarray(int)``
        intensity component of XIC
    target_rt : ``float``
        target retention time
    peak1_rt : ``float``
        closest fitted peak RT
    rt_tol : ``float``
        retention time tolerance
    peak_rts : ``list(float)``
        retention time of annotated peak(s)
    peak_heights : ``list(float)``
        height of annotated peak(s)
    peak_fwhms : ``list(float)``
        FWHM of annotated peak(s)

    Returns
    -------
    img : ``byes``
        png image of plot as bytes
    """
    rcParams['font.size'] = 6
    fig = plt.figure(figsize=(2.25, 1.5))
    ax = fig.add_subplot()
    ax.plot(xic_rt, xic_int, ls='-', c='#777777', zorder=-2, lw=1.5)
    ax.bar(prt, max(xic_int), pwt, color='b', alpha=0.1, zorder=-2)
    if pht is not None:
        # peak height can be None, in which case we just plot the trace + extraction window but NOT the peak fit
        x = xic_rt[(xic_rt >= prt - 2 * pwt) & (xic_rt <= prt + 2 * pwt)]
        ax.plot(x, _gauss(x, prt, pht, pwt), "b-", lw=1)
    # optionally add some light lines showing XIC comparison range (+/- FWHM x3)
    if add_cmp_range:
        ax.axvline(prt - 3 * pwt, c="#666666", lw=0.5, ls="--", zorder=-3)
        ax.axvline(prt + 3 * pwt, c="#666666", lw=0.5, ls="--", zorder=-3)
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    if restrict_xlim:
        ax.set_xlim([max(prt - 5 * pwt, 0), min(prt + 5 * pwt, max(xic_rt))])
    ax.set_xlabel('RT (min)')
    ax.set_ylabel('intensity')
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=250, bbox_inches='tight')
    plt.close()
    buf.seek(0)
    img = buf.read()
    buf.close()
    return img


def _find_closest_peak(target_x, peak_xs, peak_hts, peak_wts, x_threshold):
    """
    finds the closes peak to a target peak

    Parameters
    ----------
    target_x : float
        target x value
    peak_xs : list(float)
        peak x values
    peak_hts : list(float)
        peak heights
    peak_wts : list(float)
        peak FWHMs
    x_threshold : float
        maximal acceptable difference between target x value and a matched peak

    Returns
    -------
    x : float
        closest peak x value
    ht : float
        closest peak height
    wt : float
        closest peak FWHM
    """
    if len(peak_xs) == 0:
        # peaks list is empty
        return None, None, None
    i_best = 0
    dx = abs(peak_xs[0] - target_x)
    for i in range(1, len(peak_xs)):
        if abs(peak_xs[i] - target_x) < dx:
            i_best = i
            dx = abs(peak_xs[i] - target_x)
    if dx > x_threshold:
        # no peak was close enough to the target value
        return None, None, None
    return peak_xs[i_best], peak_hts[i_best], peak_wts[i_best]


def _fit_xic_rt(xic_rt, xic_int, rt_fit_method, debug_flag, debug_cb):
    """
    performs peak fitting on an XIC, returning the fitted peak parameters. Returns (None, None, None), None if anything 
    bad happens

    Parameters
    ----------
    xic_rt : ``np.ndarray(float)``
        retention time component of XIC
    xic_int : ``np.ndarray(float)``
        intensity component of XIC
    rt_fit_method : ``str``
        method to use for fitting XIC peaks, "gauss" or "localmax"
    debug_flag : ``str``
        specifies how to dispatch the message and/or plot, None to do nothing
    debug_cb : ``func``
        callback function that takes the debugging message as an argument, can be None if
        debug_flag is not set to 'textcb'

    Returns
    -------
    peak_info : ``tuple(float, float, float)``
        peak mean RT, height, and FWHM
    img : ``bytes``
        png image of plot as bytes
    """
    if len(xic_rt) < 2:  # two points needed for interpolation (right? or only 1?)
        # not enough data, cannot fit RT peak
        return (None, None, None), None
    # perform linear interpolation on the data (200 points/min)
    xic_i_rt, xic_i_int = lerp_1d(xic_rt, xic_int, min(xic_rt), max(xic_rt), 200)
    # fit xic peak
    if rt_fit_method == 'gauss':
        peaks = find_peaks_1d_gauss(xic_i_rt, xic_i_int, 0.1, _XIC_MIN_ABS_INTENSITY, 0.1, 1.0, 2, True)
    elif rt_fit_method == 'localmax':
        peaks = find_peaks_1d_localmax(xic_i_rt, xic_i_int, 0.1, _XIC_MIN_ABS_INTENSITY, 0.1, 1.0, 0.1)
    else:
        msg = '_fit_xic_rt: rt_fit_method must be "localmax" or "gauss" (was: "{}")'
        raise ValueError(msg.format(rt_fit_method))
    for prt, pht, pwt in zip(*peaks):
        img = _plot_xic_with_peak(xic_i_rt, xic_i_int, prt, pht, pwt, restrict_xlim=True)
        # (debug) report the difference between fitted and target RT, plot the peak fit
        msg = f"peak RT: {prt:.2f}"
        _debug_handler(debug_flag, debug_cb, msg=msg, img=img)
        yield (prt, pht, pwt), img


def _find_leading_edge_rt(xic_rt, xic_int, peak_rt, peak_ht, le_range, debug_flag, debug_cb):
    """
    ?

    Parameters
    ----------
    
    le_range : tuple(float)

    debug_flag : ``str``
        specifies how to dispatch the message and/or plot, None to do nothing
    debug_cb : ``func``
        callback function that takes the debugging message as an argument, can be None if
        debug_flag is not set to 'textcb'

    Returns
    -------
    le_rt_min : 
        minimum RT for leading edge of peak
    le_rt_max
        maximum RT for leading edge of peak
    """
    
    if peak_rt is None:
        # peak fitting failed
        return None, None
    # determine the RT where signal is <= (le_threshold[1] * peak height) on the leading edge of the peak
    le_t_min, le_t_max = le_range
    le_t_min *= peak_ht
    le_t_max *= peak_ht
    idx = np.abs(xic_rt - peak_rt).argmin()
    while xic_int[idx] > le_t_max:
        idx -= 1
        if idx <= 0:
            # leading edge is off the xic
            return None, None
    # when the loop is complete, set max RT 
    le_rt_max = xic_rt[idx]
    while xic_int[idx] > le_t_min and idx > 0:
        idx -= 1
    # when the loop is complete, set min RT 
    le_rt_min = xic_rt[idx]
    # (debug) print the RT bounds of the leading edge
    msg = 'RT peak leading edge bounds: {:.2f}, {:.2f}'.format(le_rt_min, le_rt_max)
    _debug_handler(debug_flag, debug_cb, msg=msg)
    return le_rt_min, le_rt_max


def _plot_ms1_isotope_dist(ms1_i_mz, ms1_i_int, mz_targets, abun_targets, peak_mzs, peak_heights, peak_fwhms):
    """
    plots the MS1 isotope distribution with fitted peaks and comparison against theoretical isotope distribution, 
    returns a png image of plot as bytes

    Parameters
    ----------
    ms1_i_mz : ``numpy.ndarray(float)``
        m/z component of MS1 spectrum
    ms1_i_int : ``numpy.ndarray(int)``
        intensity component of MS1 spectrum
    mz_targets : ``tuple(float)``
        m/z values expected from theoretical isotope distribution
    abun_targets : ``tuple(float)``
        abundances expected from theoretical isotope distribution
    peak_mzs : ``list(float)``
        m/z annotated peak(s)
    peak_heights : ``list(float)``
        height of annotated peak(s)
    peak_fwhms : ``list(float)``
        FWHM of annotated peak(s)

    Returns
    -------
    img : ``bytes``
        png image of plot as bytes
    """
    rcParams['font.size'] = 6
    fig = plt.figure(figsize=(6.25, 1.5))
    ax = fig.add_subplot()
    ax.plot(ms1_i_mz, ms1_i_int, ls='-', c='#777777', zorder=-2, lw=1.5)
    for m, i in zip(mz_targets, abun_targets,):
        ax.axvline(m, ls='--', c='r', lw=1)
        ax.plot([m - 0.15, m + 0.15], [i, i], 'r--', lw=1)
    for pmz, pht, pwt in zip(peak_mzs, peak_heights, peak_fwhms):
        if pmz is not None:
            x = ms1_i_mz[(ms1_i_mz >= pmz - 2 * pwt) & (ms1_i_mz <= pmz + 2 * pwt)]
            ax.plot(x, _gauss(x, pmz, pht, pwt), "b-", lw=1)
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.set_xlabel('m/z')
    ax.set_ylabel('intensity')
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=250, bbox_inches='tight')
    buf.seek(0)
    plt.close()
    img = buf.read()
    buf.close()
    return img


def _mz_int_isotope_score(target_mz, target_abun_scaled, mpeak_mzs, mpeak_hts):
    """
    """
    ppme = []
    for tmz, pmz in zip(target_mz, mpeak_mzs):
        if pmz is not None:
            ppme.append(_ppm_error(tmz, pmz))
    ire = []
    for tht, pht in zip(target_abun_scaled, mpeak_hts):
        if pht is not None:
            # peak error relative to the intensity of the base peak
            ire.append(100. * (pht - tht) / target_abun_scaled[0])
    return np.mean(np.abs(ppme)), np.mean(np.abs(ire))


def _calc_cosine_from_recon_spectrum(ms1_i_mz, ms1_i_int, target_mz, target_abun_scaled):
    """ 
    computes cosine similarity score between observed ms1 spectrum and a spectrum reconstituted from the theoretical
    isotopic distribution
    """
    # theoretical FWHM is a function of the system's resolving power and m/z (fwhm = mz / R)
    # typical resolving power for the systems we have been using is ~15-20k -> use 17.5k
    theo_fwhm = target_mz[0] / 17500.
    # add gaussians for each theoretical peak
    theo_int = _gauss(ms1_i_mz, target_mz[0], target_abun_scaled[0], theo_fwhm)
    theo_int += _gauss(ms1_i_mz, target_mz[1], target_abun_scaled[1], theo_fwhm)
    theo_int += _gauss(ms1_i_mz, target_mz[2], target_abun_scaled[2], theo_fwhm)
    # return the cosine distance
    return spatial.distance.cosine(ms1_i_int, theo_int)


def _calc_cosine_from_xic(pre_xic_rt, pre_xic_int, frg_xic_rt, frg_xic_int):
    """
    """
    # make sure that both XICs have more than 2 points and some of them are > 0
    # if not, then just return cosine distance of 1
    if (
        len(pre_xic_int) < 2
        or len(frg_xic_int) < 2
        or sum(pre_xic_int) == 0 
        or sum(frg_xic_int) == 0
    ):
        return 1. 
    # lerp both signals to ensure the same sampling (200 points/min)
    _, pre_xic_i_int = lerp_1d(pre_xic_rt, pre_xic_int, min(pre_xic_rt), max(pre_xic_rt), 200)
    _, frg_xic_i_int = lerp_1d(frg_xic_rt, frg_xic_int, min(pre_xic_rt), max(pre_xic_rt), 200)
    # normalize both signals
    pre_xic_i_int /= max(pre_xic_i_int)
    frg_xic_i_int /= max(frg_xic_i_int)
    return spatial.distance.cosine(pre_xic_i_int, frg_xic_i_int)


def _extract_ms1_and_calc_isotope_score(oz_data, target_mz, target_abun, rt_min, rt_max, ms1_fit_method, 
                                        debug_flag, debug_cb):
    """
    gets the mass spectrum from a specified RT range and compares theoretical and observed isotopic distribution

    Parameters
    ----------

    debug_flag : ``str``
        specifies how to dispatch the message and/or plot, None to do nothing
    debug_cb : ``func``
        callback function that takes the debugging message as an argument, can be None if
        debug_flag is not set to 'textcb'

    Returns
    -------
    scores : ``tuple(float, float, float)``
        m/z, intensity, and MS1 cosine distance scores
    iso_dist_img : ``bytes``
        png image of plot as bytes
    """
    # m/z bounds: M isotope - 1.5, M+2 isotope + 0.5
    mz_min, mz_max = target_mz[0] - 1.5, target_mz[2] + 0.5
    # load MS1 spectrum
    ms1_mz, ms1_int = oz_data.collect_ms1_arrays_by_rt(rt_min, rt_max, mz_bounds=(mz_min, mz_max))
    # make sure the MS1 spectrum actually has data in it (technically need at least 2 points for interpolation)
    if len(ms1_mz) < 2:
        # not enough data
        return (None, None, None), None
    # interpolate (500 points/Da)
    ms1_i_mz, ms1_i_int = lerp_1d(ms1_mz, ms1_int, mz_min, mz_max, 500)
    if ms1_fit_method == 'gauss':
        peak_mzs, peak_hts, peak_wts = find_peaks_1d_gauss(ms1_i_mz, ms1_i_int, 0.01, 
                                                           _MS1_MIN_ABS_INTENSITY, 0.01, 0.25, 20, False)
    elif ms1_fit_method == 'localmax':
        peak_mzs, peak_hts, peak_wts = find_peaks_1d_localmax(ms1_i_mz, ms1_i_int, 0.01, 
                                                              _MS1_MIN_ABS_INTENSITY, 0.01, 0.25, 0.05)
    else:
        msg = (
            "_extract_ms1_and_calc_isotope_score: ms1_fit_method must be "
            f"'localmax' or 'gauss' (was: '{ms1_fit_method}')"
        )
        raise ValueError(msg)
    # select the closest matching peaks for each isotope
    mpeak_mzs, mpeak_hts, mpeak_wts = [None, None, None], [None, None, None], [None, None, None]
    for i in [0, 1, 2]:
        mpeak_mzs[i], mpeak_hts[i], mpeak_wts[i] = _find_closest_peak(target_mz[i], peak_mzs, peak_hts, peak_wts, 0.1)
        # set intensity to 0 for any peaks that were not found (helps with scoring)
        mpeak_hts[i] = 0 if mpeak_hts[i] is None else mpeak_hts[i]
    # scale target abundance by M, M+1, or M+2 observed abundance, whichever is present first
    if mpeak_mzs[0] is not None:
        target_abun_scaled = [_ * mpeak_hts[0] for _ in target_abun]
    elif mpeak_mzs[1] is not None:
        target_abun_scaled = [
            mpeak_hts[1] / target_abun[1], 
            mpeak_hts[1], 
            mpeak_hts[1] * (target_abun[2] / target_abun[1])
        ]
    else:
        # not worth considering if M or M+1 were not found
        # no matching peaks found
        return (None, None, None), None
    # (debug) print out peak params, plot fit compared to isotopes
    for i in [0, 1, 2]:
        m_name = f"M{["  ", f"+{i}", f"+{i}"][i]}"
        if mpeak_mzs[i] is not None:
            mz = mpeak_mzs[i]
            dmz = mz - target_mz[i]
            ht = mpeak_hts[i]
            dht = ht - target_abun_scaled[i]
            msg = f"{m_name} m/z: {mz:9.4f} ({dmz:+6.4f}) intensity: {ht:.2e} ({dht:+.2e})"
            _debug_handler(debug_flag, debug_cb, msg=msg)
        else:
            _debug_handler(debug_flag, debug_cb, msg=m_name + ' peak not found')
    isotope_dist_img = _plot_ms1_isotope_dist(ms1_i_mz, ms1_i_int, target_mz, target_abun_scaled, mpeak_mzs, mpeak_hts, mpeak_wts)
    mz_score, int_score = _mz_int_isotope_score(target_mz, target_abun_scaled, mpeak_mzs, mpeak_hts)
    cos_score = _calc_cosine_from_recon_spectrum(ms1_i_mz, ms1_i_int, target_mz, target_abun_scaled)
    # (debug)
    msg = f"mz ppm: {mz_score:.3f} abundance percent: {int_score:.3f} cos distance: {cos_score:.3f}"
    _debug_handler(debug_flag, debug_cb, msg=msg, img=isotope_dist_img)
    return (mz_score, int_score, cos_score), isotope_dist_img


def _do_fragment_scoring(oz_data: OzData, 
                         fragment: str, 
                         formula: Formula, 
                         mz_tol: float, 
                         ms1_fit_method: str, 
                         pre_rt: float,
                         pre_wt: float,
                         pre_xic_rt: Any, 
                         pre_xic_int: Any, 
                         debug_flag: Optional[str], 
                         debug_cb: Optional[Any]
                         ) -> Optional[Dict[str, Any]] :
    """
    Perform scoring analysis for a single fragment

    Parameters
    ----------
    oz_data
        interface for OzID data in either MZA or UIMF
    fragment
        "aldehyde" or "criegee"
    formula
        molecular formula of the fragment
    mz_tol
        m/z tolerance for data extraction
    ms1_fit_method
        "gauss" or "localmax"
    pre_rt, pre_wt
        precursor retention time and peak width
    pre_xic_rt, pre_xic_int
        precursor XIC arrays
    debug_flag 
        specifies how to dispatch the message and/or plot, None to do nothing
    debug_cb 
        callback function that takes the debugging message as an argument, can be None if
        debug_flag is not set to 'textcb'

    Returns
    -------
    results
        individual fragment scoring results, or None 
    """
    fragment_results = {}
    # scoring for aldehyde
    _debug_handler(debug_flag, debug_cb, msg=fragment)
    target_mz, target_abun = predict_m_m1_m2(formula)
    # get XIC using fragment m/z +/- mz_tol and RT bounds from the precursor peak
    mz_min, mz_max = target_mz[0] - mz_tol, target_mz[0] + mz_tol
    # RT range for extracting XIC (precursor RT +/- FWHM x4)
    xic_rtb = (pre_rt - pre_wt * 4, pre_rt + pre_wt * 4) 
    # RT range for extracting MS1 (precursor RT +/- FWHM /2)
    ms1_rtb = (pre_rt - pre_wt / 2, pre_rt + pre_wt / 2) 
    xic = oz_data.collect_xic_arrays_by_mz(mz_min, mz_max, rt_bounds=xic_rtb)
    if sum(xic[1]) > 0 and (pre_scoring := _extract_ms1_and_calc_isotope_score(oz_data, target_mz, target_abun, 
                                                                               *ms1_rtb, ms1_fit_method, 
                                                                               debug_flag, debug_cb)) is not None:
        # unpack and add into results
        (mz_ppm, abun_pct, mz_cos_dist), iso_dist_img = pre_scoring
        # determine the XIC comparison range based on peak width (+/- FWHM x3)
        pre_xic_cmp_idx = (pre_xic_rt >= pre_rt - pre_wt * 3) & (pre_xic_rt <= pre_rt + pre_wt * 3)
        xic_cmp_idx = (xic[0] >= pre_rt - pre_wt * 3) & (xic[0] <= pre_rt + pre_wt * 3)
        xic_peak_ht = (
            max(a) 
            if len(a := xic[1][(xic[0] >= ms1_rtb[0]) & (xic[0] <= ms1_rtb[1])]) > 0
            else 0.
        )
        fragment_results = {
            'target_mz': target_mz[0],
            'target_rt': None,              # * no longer doing XIC peak fitting for fragments
            'xic_peak_rt': None,            # * 
            'xic_peak_ht': xic_peak_ht,
            'xic_peak_fwhm': None,          # *
            'mz_ppm': mz_ppm,               
            'abun_percent': abun_pct,
            'mz_cos_dist': mz_cos_dist if mz_cos_dist is not None else 1.,
            'rt_cos_dist': _calc_cosine_from_xic(
                pre_xic_rt[pre_xic_cmp_idx], 
                pre_xic_int[pre_xic_cmp_idx],
                xic[0][xic_cmp_idx], 
                xic[1][xic_cmp_idx]
            ),
            'isotope_dist_img': iso_dist_img,
            'xic_fit_img': _plot_xic_with_peak(*xic, pre_rt, None, pre_wt, add_cmp_range=True),
            'saturation_corrected': False,  # no more saturation correction
        }
        return fragment_results
    else:  # if unable to do isotope scoring, return None
        return None


def _do_fragment_scoring_infusion(oz_data: OzData, 
                                  fragment: str, 
                                  formula: Formula, 
                                  ms1_fit_method: str, 
                                  debug_flag: Optional[str], 
                                  debug_cb: Optional[Any]
                                  ) -> Dict[str, Any] :
    """
    """
    fragment_results = {}
    # scoring for aldehyde
    _debug_handler(debug_flag, debug_cb, fragment)
    target_mz, target_abun = predict_m_m1_m2(formula)
    scores, iso_dist_img = _extract_ms1_and_calc_isotope_score(oz_data, 
                                                               target_mz, 
                                                               target_abun, 
                                                               oz_data.min_rt, 
                                                               oz_data.max_rt, 
                                                               ms1_fit_method, 
                                                               debug_flag, 
                                                               debug_cb)
    if scores is not None:
        fragment_results = {
            'target_mz': target_mz[0],
            'mz_ppm': scores[0],
            'abun_percent': scores[1],
            'mz_cos_dist': scores[2],
            'isotope_dist_img': iso_dist_img,
        }
        return fragment_results


# TODO (Dylan Ross): Look for second-order OzID fragments if at least two chains were specified and at least 
#                    two contain unsaturations. Like with first-order OzID fragments, compute the unique set 
#                    of all second order fragments, each defined by a pair of db_idx values and corresponding 
#                    pair of db_pos values. Probably good to make a helper function that performs this calculation 
#                    that takes the lists of fa_carbons and fa_unsats as input and outputs a set of these unique
#                    pairs. Also need to implement the actual fragment formula calculations, but I am pretty sure 
#                    it should be just as simple as sequentially calling the predict OzID fragment formulas, using
#                    either the aldehyde or criegee formula as input on the second round. But there will have to
#                    be some pruning of these formulas to avoid double counting of the AC/CA species.


def score_db_pos_isotope_dist_polyunsat(oz_data: OzData, 
                                        precursor_formula: Formula, 
                                        fa_nc: Tuple[int], 
                                        fa_nu: Tuple[int], 
                                        mz_tol: float, 
                                        rt_fit_method: str = 'gauss', 
                                        ms1_fit_method: str = 'localmax', 
                                        check_saturation: bool = False, 
                                        saturation_threshold: float = 2e6, 
                                        remove_d: Optional[int] = None, 
                                        debug_flag: Optional[str] = None,
                                        debug_cb: Optional[Any] = None, 
                                        info_cb: Optional[Any] = None, 
                                        early_stop_event: Optional[threading.Event] = None
                                        ) -> Generator[Dict[str, Any], Any, Any] :
    """
    performs isotope distribution scoring for a range of potential double-bond positions for polyunsaturated lipids, 
    also works for monounsaturated lipids, and essentially does nothing for completely saturated lipids

    .. warning:: 
        Checking for and correcting signal saturation is not implemented since the 
        updates to XIC processing for precursor/fragments. Enabling the `check_saturation` 
        option will cause an exception (`AssertionError`).

    Parameters
    ----------
    oz_data 
        UIMF data interface instance for OzID data
    precursor_formula
        chemical formula of the precursor ion
    fa_nc 
        number of carbons in precursor fatty acids
    fa_nu 
        number of DB in each precursor fatty acid, in same order as precursor_nc
    mz_tol 
        m/z tolerance for extracting XICs
    rt_fit_method 
        specify method to use for fitting the RT peak ('gauss' works best in testing)
    ms1_fit_method 
        specify method to use for fitting MS1 spectrum ('localmax' works best in testing)
    check_saturation 
        whether to check for signal saturation and use leading edge strategy if necessary
        NOTE: saturation correction is not currently implemented since last round of updates
    saturation_threshold 
        specify a threshold intensity for determining peak saturation
    remove_d 
        adjust molecular formulas to get rid of D labels on fatty acid tail that are part of the neutral loss
        (specific to SPLASH lipids)
    debug_flag 
        specifies how to dispatch the message and/or plot, None to do nothing
    debug_cb 
        callback function that takes the debugging message as an argument, can be None if
        debug_flag is not set to 'textcb'
    info_cb 
        optional callback function that gets called at several intermediate steps and gives information about data
        processing details. Callback function takes a single argument which is a ``str`` info message
    early_stop_event 
        When the workflow is running in its own thread and this event gets set, processing is stopped gracefully

    Yields
    ------
    result
        dictionary containing analysis results
    """
    # use non GUI backend if debug is set to False (required for lipidoz_gui)
    if debug_flag != 'full':
        mpuse('Agg') 
    _debug_handler(debug_flag, debug_cb, msg='----------precursor----------')
    # predict precursor isotope distribution
    pre_target_mz, pre_target_abun = predict_m_m1_m2(precursor_formula)
    _debug_handler(debug_flag, debug_cb, msg=f"mz={pre_target_mz[0]:.4f}")
    # compute bounds of db indices and positions, keep only unique pairs
    # this does not need to be repeated for every XIC peak, so just do it once ahead of time
    dbidx_dbpos = set()
    for fac, fau in zip(fa_nc, fa_nu):
        if fac > 0:
            for db_idx in range(1, fau + 1):
                dbp_min, dbp_max = _calc_dbp_bounds(fac, fau, db_idx)
                for db_pos in range(dbp_min, dbp_max + 1):
                    dbidx_dbpos.add((db_idx, db_pos))
    n_combos = len(dbidx_dbpos)
    if info_cb is not None:
        msg = 'INFO: there are {} combinations of (db_idx, db_pos) to examine'
        info_cb(msg.format(n_combos))
    # get XIC for  target_rt +/- 0.5 using precursor m/z +/- mz_tol
    mz_min, mz_max = pre_target_mz[0] - mz_tol, pre_target_mz[0] + mz_tol
    pre_xic = oz_data.collect_xic_arrays_by_mz(mz_min, mz_max)
    # flag: was a peak found for the precursor?
    found_peak = False
    # get the fitted peak info from the xic
    for (peak_rt, peak_ht, peak_wt), pre_xic_fit_img in _fit_xic_rt(*pre_xic, rt_fit_method, 
                                                                    debug_flag, debug_cb):
        if info_cb is not None:
            msg = f"INFO: precursor peak RT = {peak_rt:.2f} min"
            info_cb(msg)
        # check for a stop event and break the loop if the early stop event gets set
        if early_stop_event is not None and early_stop_event.is_set():
            break
        # set flag: peak was found
        found_peak = True
        # NOTE: Results starts empty, then different pieces of info get added to it and if a sufficient
        #       portion of the total processing is successful it gets yielded. Otherwise nothing is yielded. 
        #       separate results get yielded for each XIC peak that was found.
        results = {}
        # ----------------------------------------------------
        # TODO (Dylan Ross): For now, disable checking/correcting signal saturation. This process is made more 
        #                    complicated by the fact that now fragments will be extracted from static RT ranges
        #                    instead of from separate fits on their XICs. Because of this, the logic for checking
        #                    and correcting signal saturation in the precursor would need to be modified.
        if check_saturation:
            msg = (
                "score_db_pos_isotope_dist_polyunsat: signal saturation checking/correction has not been "
                "reimplemented since changing how XICs are processed for precursor/fragments"
            )
            assert False, msg
        # check for and potentially correct signal saturation
        # TODO (Dylan Ross): make this check more sophisticated, detect peak saturation from peak shape or something
        # is_saturated = peak_ht > saturation_threshold
        # if check_saturation and is_saturated:
        #     msg = 'peak saturation detected, only using scans from leading edge of RT peak for isotope scoring'
        #     _debug_handler(debug_flag, debug_cb, msg=msg)
        #     # find the RT at the leading edge of the peak (RT scans from between 5% and 25% peak height)
        #     rt_min, rt_max = _find_leading_edge_rt(pre_xic_rt, pre_xic_int, peak_rt, peak_ht, peak_wt, (0.05, 0.25), 
        #                                         debug_flag, debug_cb)
        # else:
        #     rt_min, rt_max = peak_rt - rt_tol, peak_rt + rt_tol
        # if rt_min is None or rt_max is None:
        #     # finding the leading edge failed, do not apply correction
        #     msg = 'failed to select scans from leading edge of RT peak, reverting to original RT bounds'
        #     _debug_handler(debug_flag, debug_cb, msg=msg)
        #     rt_min, rt_max = peak_rt - rt_tol, peak_rt + rt_tol
        #     is_saturated = False
        # ----------------------------------------------------
        # extract MS1 and score (using RT bounds)
        if (pre_scoring := _extract_ms1_and_calc_isotope_score(oz_data, pre_target_mz, pre_target_abun, 
                                                               peak_rt - peak_wt / 2, peak_rt + peak_wt / 2, 
                                                               ms1_fit_method, 
                                                               debug_flag, debug_cb)) is not None:
            # unpack and add into results
            (pre_mz_ppm, pre_abun_pct, pre_mz_cos_dist), pre_iso_dist_img = pre_scoring
            results['precursor'] = {
                'target_mz': pre_target_mz[0],
                'target_rt': None,                  # not using target RT anymore
                'xic_peak_rt': peak_rt,
                'xic_peak_ht': peak_ht,
                'xic_peak_fwhm': peak_wt,
                'mz_ppm': pre_mz_ppm, 
                'abun_percent': pre_abun_pct, 
                'mz_cos_dist': pre_mz_cos_dist,
                'isotope_dist_img': pre_iso_dist_img,
                'xic_fit_img': pre_xic_fit_img,
                'saturation_corrected': False  # as noted above, this is currently disabled
            }
            # iterate over (db_idx, db_pos) pairs and perform analysis on the corresponding fragments
            results['fragments'] = {}
            i = 1
            for i, (db_idx, db_pos) in enumerate(dbidx_dbpos):
                # check for a stop event and break the loop if the early stop event gets set
                # this will eventually fall through and not yield anything so that should be
                # sufficient, right? Yes.
                if early_stop_event is not None and early_stop_event.is_set():
                    break
                # ITERATE THROUGH DB POSITIONS AND LOOK FOR FRAGMENTS
                #----------------------------------------------------
                if db_idx not in results['fragments']:
                    results['fragments'][db_idx] = {}
                results['fragments'][db_idx][db_pos] = {}
                msg = f"----------{db_idx=},{db_pos=}----------"
                _debug_handler(debug_flag, debug_cb, msg=msg)
                # get aldehyde and criegee formulas
                ald_formula, crg_formula = _polyunsat_ald_crg_formula(precursor_formula, db_pos, db_idx)
                if remove_d is not None:
                    # adjust the molecular formulas to get rid of D7 labels in neutral loss fragments
                    # applicable to SPLASH lipid standards only
                    ald_formula.pop('D')
                    ald_formula['H'] += remove_d
                    crg_formula.pop('D')
                    crg_formula['H'] += remove_d
                # scoring for aldehyde
                ald_results = _do_fragment_scoring(oz_data, 'aldehyde', ald_formula, 
                                                   mz_tol, ms1_fit_method,
                                                   peak_rt, peak_wt, *pre_xic, 
                                                   debug_flag, debug_cb)
                # scoring for criegee
                crg_results = _do_fragment_scoring(oz_data, 'criegee', crg_formula, 
                                                   mz_tol, ms1_fit_method,
                                                   peak_rt, peak_wt, *pre_xic, 
                                                   debug_flag, debug_cb)
                # NOTE: There is no use keeping fragment results that are incomplete (i.e., no
                #       XIC fit image or no isotope distribution image)
                ald_results = ald_results if (ald_results is not None
                                              and ald_results["isotope_dist_img"] is not None 
                                              and ald_results["xic_fit_img"] is not None) else None
                crg_results = crg_results if (crg_results is not None
                                              and crg_results["isotope_dist_img"] is not None 
                                              and crg_results["xic_fit_img"] is not None) else None
                results['fragments'][db_idx][db_pos]['aldehyde'] = ald_results 
                results['fragments'][db_idx][db_pos]['criegee'] = crg_results
                if info_cb is not None:
                    msg = f"INFO: {i + 1} of {n_combos} done"
                    info_cb(msg)
                #----------------------------------------------------
            # finally, yield the results for this target + XIC peak
            yield results
        else:  # if unable to do precursor scoring, report that 
            _debug_handler(debug_flag, debug_cb, msg='unable to determine precursor scores')
            if info_cb is not None and debug_flag is None:  # avoid duplication of messages
                info_cb('INFO: unable to determine precursor scores')
            #return None - now we yield results, so yield nothing instead
    # report that no peak was found for the precursor if the flag was not set
    if not found_peak:
        _debug_handler(debug_flag, debug_cb, msg='unable to find RT peak for precursor')
        if info_cb is not None and debug_flag is None:  # avoid duplication of messages
            info_cb('INFO: unable to find RT peak for precursor')
        #return None - now we yield results, so yield nothing instead
    _debug_handler(debug_flag, debug_cb, msg='------------------------------')


def score_db_pos_isotope_dist_polyunsat_infusion(oz_data, precursor_formula, fa_nc, fa_nu, mz_tol, 
                                                 ms1_fit_method='localmax', remove_d=None, 
                                                 debug_flag=None, debug_cb=None):
    """
    works the same as score_db_pos_isotope_dist_polyunsat but for direct infusion data (*i.e.* without retention 
    time information). All components of the analysis related to retention time are omitted.

    Parameters
    ----------
    oz_data : ``mzapy.MZA``
        mza data interface instance for OzID data
    precursor_formula : ``dict(str:int)``
        chemical formula of the precursor ion
    fa_nc : ``tuple(int)``
        number of carbons in precursor fatty acids
    fa_nu : ``tuple(int)``
        number of DB in each precursor fatty acid, in same order as precursor_nc
    mz_tol : ``float``
        m/z tolerance for extracting XICs
    ms1_fit_method : ``str``, default='localmax'
        specify method to use for fitting MS1 spectrum ('localmax' works best in testing)
    remove_d : ``int``, optional
        adjust molecular formulas to get rid of D labels on fatty acid tail that are part of the neutral loss
        (specific to SPLASH lipids)
    debug_flag : ``str``, optional
        specifies how to dispatch the message and/or plot, None to do nothing
    debug_cb : ``func``, optional
        callback function that takes the debugging message as an argument, can be None if
        debug_flag is not set to 'textcb'

    Returns
    -------
    result : dict(...)
        dictionary containing analysis results
    """
    # has not been updated to account for changes from working with UIMF data
    assert False, "not implemented"
    # use non GUI backend if debug is set to False (required for lipidoz_gui)
    if debug_flag != 'full':
        mpuse('Agg') 
    results = {}
    _debug_handler(debug_flag, debug_cb, msg='----------precursor----------')
    # predict precursor isotope distribution
    pre_target_mz, pre_target_abun = predict_m_m1_m2(precursor_formula)
    # extract MS1 and score 
    pre_scores, pre_iso_dist_img = _extract_ms1_and_calc_isotope_score(oz_data, pre_target_mz, pre_target_abun, 
                                                                       oz_data.min_rt, oz_data.max_rt, 
                                                                       ms1_fit_method, 
                                                                       debug_flag, debug_cb)
    if pre_scores is None:
        _debug_handler(debug_flag, debug_cb, msg='unable to determine precursor scores')
        return None
    results['precursor'] = {
        'target_mz': pre_target_mz[0],
        'mz_ppm': pre_scores[0], 
        'abun_percent': pre_scores[1], 
        'mz_cos_dist': pre_scores[2],
        'isotope_dist_img': pre_iso_dist_img,
    }
    # compute bounds of db indices and positions, keep only unique pairs
    dbidx_dbpos = set()
    for fac, fau in zip(fa_nc, fa_nu):
        if fac > 0:
            for db_idx in range(1, fau + 1):
                dbp_min, dbp_max = _calc_dbp_bounds(fac, fau, db_idx)
                for db_pos in range(dbp_min, dbp_max + 1):
                    dbidx_dbpos.add((db_idx, db_pos))
    # iterate over (db_idx, db_pos) pairs and perform analysis
    results['fragments'] = {}
    for db_idx, db_pos in dbidx_dbpos:
        # ITERATE THROUGH DB POSITIONS AND LOOK FOR FRAGMENTS
        #---------------------------------------------------------------------------------------------------
        if db_idx not in results['fragments']:
            results['fragments'][db_idx] = {}
        results['fragments'][db_idx][db_pos] = {}
        _debug_handler(debug_flag, debug_cb, msg='\n----------db_idx={},db_pos={}----------'.format(db_idx, db_pos))
        # get aldehyde and criegee formulas
        ald_formula, crg_formula = _polyunsat_ald_crg_formula(precursor_formula, db_pos, db_idx)
        if remove_d is not None:
            # adjust the molecular formulas to get rid of D7 labels in neutral loss fragments
            # applicable to SPLASH lipid standards only
            ald_formula.pop('D')
            ald_formula['H'] += remove_d
            crg_formula.pop('D')
            crg_formula['H'] += remove_d
        # scoring for aldehyde
        ald_results = _do_fragment_scoring_infusion(oz_data, 'aldehyde', ald_formula, mz_tol, ms1_fit_method, 
                                                    debug_flag, debug_cb)
        # scoring for criegee
        crg_results = _do_fragment_scoring_infusion(oz_data, 'criegee', crg_formula, mz_tol, ms1_fit_method, 
                                                    debug_flag, debug_cb)
        results['fragments'][db_idx][db_pos]['aldehyde'] = ald_results
        results['fragments'][db_idx][db_pos]['criegee'] = crg_results
        #---------------------------------------------------------------------------------------------------
    _debug_handler(debug_flag, debug_cb, msg='------------------------------\n\n')
    # TODO (Dylan Ross): Look for second-order OzID fragments if at least two chains were specified and at least 
    #                    two contain unsaturations. Like with first-order OzID fragments, compute the unique set 
    #                    of all second order fragments, each defined by a pair of db_idx values and corresponding 
    #                    pair of db_pos values. Probably good to make a helper function that performs this calculation 
    #                    that takes the lists of fa_carbons and fa_unsats as input and outputs a set of these unique
    #                    pairs. Also need to implement the actual fragment formula calculations, but I am pretty sure 
    #                    it should be just as simple as sequentially calling the predict OzID fragment formulas, using
    #                    either the aldehyde or criegee formula as input on the second round. But there will have to
    #                    be some pruning of these formulas to avoid double counting of the AC/CA species.
    return results


# TODO: This function seems like it might not have the updates to XIC extraction that the untargeted variant has?
#       Maybe the logic after determining double bond positions from both functions should be factored into a
#       separate internal function that both can call down into?

def score_db_pos_isotope_dist_targeted(oz_data: OzData, 
                                       precursor_formula: Formula, 
                                       db_idxs: List[int], 
                                       db_posns: List[int], 
                                       target_rt: float, 
                                       rt_tol: float, 
                                       rt_peak_win: float, 
                                       mz_tol: float, 
                                       rt_fit_method: str = 'gauss', 
                                       ms1_fit_method: str = 'localmax', 
                                       check_saturation: bool = True, 
                                       saturation_threshold: float = 1e5, 
                                       remove_d: Optional[int] = None, 
                                       debug_flag: Optional[str] = None, 
                                       debug_cb: Optional[Any] = None, 
                                       info_cb: Optional[Any] = None
                                       ) -> Dict[str, Any] :
    """
    performs isotope distribution scoring for targeted double bond positions

    Parameters
    ----------
    oz_data 
        UIMF data interface instance for OzID data
    precursor_formula
        chemical formula of the precursor ion
    db_idxs 
        list of targeted double bond indices
    db_posns 
        list of targeted double bond positions
    target_rt 
        precursor target retention time
    rt_tol 
        retention time tolerance
    rt_peak_win 
        size of RT window to extract for peak fitting
    mz_tol 
        m/z tolerance for extracting XICs
    rt_fit_method
        specify method to use for fitting the RT peak ('gauss' works best in testing)
    ms1_fit_method 
        specify method to use for fitting MS1 spectrum ('localmax' works best in testing)
    check_saturation 
        whether to check for signal saturation and use leading edge strategy if necessary
    saturation_threshold 
        specify a threshold intensity for determining peak saturation
    remove_d 
        adjust molecular formulas to get rid of D labels on fatty acid tail that are part of the neutral loss
        (specific to SPLASH lipids)
    debug_flag 
        specifies how to dispatch the message and/or plot, None to do nothing
    debug_cb 
        callback function that takes the debugging message as an argument, can be None if
        debug_flag is not set to 'textcb'
    info_cb 
        optional callback function that gets called at several intermediate steps and gives information about data
        processing details. Callback function takes a single argument which is a ``str`` info message

    Returns
    -------
    result 
        dictionary containing analysis results
    """
    # use non GUI backend if debug is set to False (required for lipidoz_gui)
    if debug_flag != 'full':
        mpuse('Agg') 
    results = {}
    _debug_handler(debug_flag, debug_cb, msg='----------precursor----------')
    # predict precursor isotope distribution
    pre_target_mz, pre_target_abun = predict_m_m1_m2(precursor_formula)
    # get XIC for  target_rt +/- 0.5 using precursor m/z +/- mz_tol
    mz_min, mz_max = pre_target_mz[0] - mz_tol, pre_target_mz[0] + mz_tol
    rtb = (target_rt - rt_peak_win, target_rt + rt_peak_win)
    pre_xic_rt, pre_xic_int = oz_data.collect_xic_arrays_by_mz(mz_min, mz_max, rt_bounds=rtb)
    # get the fitted peak info from the xic
    peak_info, pre_xic_fit_img = _fit_xic_rt(pre_xic_rt, pre_xic_int, pre_target_mz[0], mz_tol, 
                                             target_rt, rt_tol, rt_fit_method, 1., \
                                             debug_flag, debug_cb)
    peak_rt, peak_ht, peak_wt = peak_info
    if peak_rt is None:
        _debug_handler(debug_flag, debug_cb, msg='unable to find RT peak for precursor')
        if info_cb is not None and debug_flag is None:  # avoid duplication of messages
            info_cb('INFO: unable to find RT peak for precursor')
        return None
    # determine whether peak is saturated, change RT bounds used in scoring if so
    # TODO (Dylan Ross) make this check more sophisticated, detect peak saturation from peak shape or something
    is_saturated = peak_ht > saturation_threshold
    if check_saturation and is_saturated:
        msg = 'peak saturation detected, only using scans from leading edge of RT peak for isotope scoring'
        _debug_handler(debug_flag, debug_cb, msg=msg)
        # find the RT at the leading edge of the peak (RT scans from between 5% and 25% peak height)
        rt_min, rt_max = _find_leading_edge_rt(pre_xic_rt, pre_xic_int, peak_rt, peak_ht, peak_wt, (0.05, 0.25), 
                                               debug_flag, debug_cb)
    else:
        rt_min, rt_max = peak_rt - rt_tol, peak_rt + rt_tol
    if rt_min is None or rt_max is None:
        # finding the leading edge failed, do not apply correction
        msg = 'failed to select scans from leading edge of RT peak, reverting to original RT bounds'
        _debug_handler(debug_flag, debug_cb, msg=msg)
        rt_min, rt_max = peak_rt - rt_tol, peak_rt + rt_tol
        is_saturated = False
    # extract MS1 and score (using RT bounds)
    pre_scores, pre_iso_dist_img = _extract_ms1_and_calc_isotope_score(oz_data, pre_target_mz, pre_target_abun, 
                                                                       rt_min, rt_max, ms1_fit_method, 
                                                                       debug_flag, debug_cb)
    if pre_scores is None:
        _debug_handler(debug_flag, debug_cb, msg='unable to determine precursor scores')
        if info_cb is not None and debug_flag is None:  # avoid duplication of messages
            info_cb('INFO: unable to determine precursor scores')
        return None
    results['precursor'] = {
        'target_mz': pre_target_mz[0],
        'target_rt': target_rt,
        'xic_peak_rt': peak_rt,
        'xic_peak_ht': peak_ht,
        'xic_peak_fwhm': peak_wt,
        'mz_ppm': pre_scores[0], 
        'abun_percent': pre_scores[1], 
        'mz_cos_dist': pre_scores[2],
        'isotope_dist_img': pre_iso_dist_img,
        'xic_fit_img': pre_xic_fit_img,
        'saturation_corrected': check_saturation and is_saturated
    }
    results['fragments'] = {}
    n_combos = len(db_idxs)
    for i, (db_idx, db_pos) in enumerate(zip(db_idxs, db_posns)):
        # ITERATE THROUGH DB POSITIONS AND LOOK FOR FRAGMENTS
        #---------------------------------------------------------------------------------------------------
        if db_idx not in results['fragments']:
            results['fragments'][db_idx] = {}
        results['fragments'][db_idx][db_pos] = {}
        msg = '----------db_idx={},db_pos={}----------'.format(db_idx, db_pos)
        _debug_handler(debug_flag, debug_cb, msg=msg)
        # get aldehyde and criegee formulas
        ald_formula, crg_formula = _polyunsat_ald_crg_formula(precursor_formula, db_pos, db_idx)
        if remove_d is not None:
            # adjust the molecular formulas to get rid of D7 labels in neutral loss fragments
            # applicable to SPLASH lipid standards only
            ald_formula.pop('D')
            ald_formula['H'] += remove_d
            crg_formula.pop('D')
            crg_formula['H'] += remove_d
        # scoring for aldehyde
        ald_results = _do_fragment_scoring(oz_data, 'aldehyde', ald_formula, mz_tol, ms1_fit_method, rtb,
                                           rt_fit_method, peak_rt, rt_tol, pre_xic_rt, pre_xic_int, 
                                           check_saturation, saturation_threshold, 
                                           debug_flag, debug_cb)
        # scoring for criegee
        crg_results = _do_fragment_scoring(oz_data, 'criegee', crg_formula, mz_tol, ms1_fit_method, rtb, 
                                           rt_fit_method, peak_rt, rt_tol, pre_xic_rt, pre_xic_int, 
                                           check_saturation, saturation_threshold, 
                                           debug_flag, debug_cb)
        results['fragments'][db_idx][db_pos]['aldehyde'] = ald_results
        results['fragments'][db_idx][db_pos]['criegee'] = crg_results
        if info_cb is not None:
            msg = 'INFO: {} of {} done'.format(i + 1, n_combos)
            info_cb(msg)
            i += 1
        #---------------------------------------------------------------------------------------------------
    _debug_handler(debug_flag, debug_cb, msg='------------------------------')
    return results

