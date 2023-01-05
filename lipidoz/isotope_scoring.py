"""
lipidoz/isotope_scoring.py

Dylan Ross (dylan.ross@pnnl.gov)

    module for performing isotope scoring of aldehyde/criegee ions

    method:
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

from matplotlib import pyplot as plt, rcParams, image as mpimage, use as mpuse
import numpy as np
from scipy import spatial

from mzapy.isotopes import predict_m_m1_m2
from mzapy.peaks import lerp_1d, find_peaks_1d_localmax, find_peaks_1d_gauss, _gauss
from mzapy._util import _ppm_error

from lipidoz._util import _polyunsat_ald_crg_formula, _calc_dbp_bounds, _debug_handler


def _plot_xic_peaks(xic_rt, xic_int, target_rt, peak1_rt, rt_tol, peak_rts, peak_heights, peak_fwhms):
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
    rcParams['font.size'] = 10
    fig = plt.figure(figsize=(2.25, 1.5))
    ax = fig.add_subplot()
    ax.plot(xic_rt, xic_int, ls='-', c='#777777', zorder=-2, lw=1.5)
    ax.axvline(target_rt, ls='--', c='r', lw=1, zorder=-1)
    ax.bar(target_rt, max(xic_int), rt_tol * 2, color='r', alpha=0.1, zorder=-2)
    if peak1_rt is not None:
        ax.bar(peak1_rt, max(xic_int), rt_tol * 2, color='b', alpha=0.1, zorder=-2)
    for prt, pht, pwt in zip(peak_rts, peak_heights, peak_fwhms):
        ax.plot([prt, prt], [0, pht], 'b-', lw=1)
        ax.plot([prt - (pwt / 2), prt + (pwt / 2)], [pht / 2, pht / 2], 'b-', lw=1)
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.set_xlabel('RT (min)')
    ax.set_ylabel('intensity')
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=300, bbox_inches='tight')
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


def _fit_xic_rt(xic_rt, xic_int, target_mz, mz_tol, target_rt, rt_tol, rt_fit_method, rt_win, debug_flag, debug_cb):
    """
    performs peak fitting on an XIC, returning the fitted peak parameters. Returns (None, None, None), None if anything 
    bad happens

    Paramters
    ---------
    xic_rt : ``np.ndarray(float)``
        retention time component of XIC
    xic_int : ``np.ndarray(float)``
        intensity component of XIC
    target_mz : ``float``
        target m/z
    mz_tol : ``float``
        m/z tolerance
    target_rt : ``float``
        target retention time
    rt_tol : ``float``
        retention time tolerance
    rt_fit_method : ``str``
        specify method to use for fitting the RT peak
    rt_win : ``float``
        size of window around target RT to extract for peak fitting
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
    # perform linear interpolation on the data (400 points/min)
    xic_i_rt, xic_i_int = lerp_1d(xic_rt, xic_int, target_rt - rt_win / 2., target_rt + rt_win / 2., 400)
    # fit xic peak
    if rt_fit_method == 'gauss':
        peak_rts, peak_hts, peak_wts = find_peaks_1d_gauss(xic_i_rt, xic_i_int, 0.25, 1000, 0.1, 1.0, 5, False)
    elif rt_fit_method == 'localmax':
        peak_rts, peak_hts, peak_wts = find_peaks_1d_localmax(xic_i_rt, xic_i_int, 0.25, 1000, 0.1, 1.0, 0.1)
    else:
        msg = '_fit_xic_rt: rt_fit_method must be "localmax" or "gauss" (was: "{}")'
        raise ValueError(msg.format(rt_fit_method))
    # select the peak with RT that is closest to target_rt
    peak_rt, peak_ht, peak_wt = _find_closest_peak(target_rt, peak_rts, peak_hts, peak_wts, 0.5)
    img = _plot_xic_peaks(xic_i_rt, xic_i_int, target_rt, peak_rt, rt_tol, peak_rts, peak_hts, peak_wts)
    if peak_rt is None:
        # unable to fit RT peak
        return (None, None, None), None
    # (debug) report the difference between fitted and target RT, plot the peak fit
    msg = 'target RT: {:.2f} peak RT: {:.2f} difference: {:+.2f}'.format(target_rt, peak_rt, peak_rt - target_rt)
    _debug_handler(debug_flag, debug_cb, msg=msg, img=img)
    return (peak_rt, peak_ht, peak_wt), img


def _find_leading_edge_rt(xic_rt, xic_int, peak_rt, peak_ht, peak_wt, le_range, debug_flag, debug_cb):
    """
    ?

    Paramters
    ---------
    
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
    rcParams['font.size'] = 10
    fig = plt.figure(figsize=(6.25, 1.5))
    ax = fig.add_subplot()
    ax.plot(ms1_i_mz, ms1_i_int, ls='-', c='#777777', zorder=-2, lw=1.5)
    for m, i in zip(mz_targets, abun_targets,):
        ax.axvline(m, ls='--', c='r', lw=1)
        ax.plot([m - 0.15, m + 0.15], [i, i], 'r--', lw=1)
    for pmz, pht, pwt in zip(peak_mzs, peak_heights, peak_fwhms):
        if pmz is not None:
            ax.plot([pmz, pmz], [0, pht], 'b-', lw=1)
            ax.plot([pmz - (pwt / 2), pmz + (pwt / 2)], [pht / 2, pht / 2], 'b-', lw=1)
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.set_xlabel('m/z')
    ax.set_ylabel('intensity')
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=300, bbox_inches='tight')
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
    # lerp both signals to ensure the same sampling (400 points/min)
    pre_xic_i_rt, pre_xic_i_int = lerp_1d(pre_xic_rt, pre_xic_int, min(pre_xic_rt), max(pre_xic_rt), 400)
    frg_xic_i_rt, frg_xic_i_int = lerp_1d(frg_xic_rt, frg_xic_int, min(pre_xic_rt), max(pre_xic_rt), 400)
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
        return None, None
    # interpolate (1000 points/Da)
    ms1_i_mz, ms1_i_int = lerp_1d(ms1_mz, ms1_int, mz_min, mz_max, 1000)
    if ms1_fit_method == 'gauss':
        peak_mzs, peak_hts, peak_wts = find_peaks_1d_gauss(ms1_i_mz, ms1_i_int, 0.01, 1000, 0.01, 0.25, 20, False)
    elif ms1_fit_method == 'localmax':
        peak_mzs, peak_hts, peak_wts = find_peaks_1d_localmax(ms1_i_mz, ms1_i_int, 0.01, 1000, 0.01, 0.25, 0.05)
    else:
        msg = '_extract_ms1_and_calc_isotope_score: ms1_fit_method must be "localmax" or "gauss" (was: "{}")'
        raise ValueError(msg.format(ms1_fit_method))
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
        target_abun_scaled = [mpeak_hts[1] / target_abun[1], mpeak_hts[1], mpeak_hts[1] * (target_abun[2] / target_abun[1])]
    elif mpeak_mzs[2] is not None:
        target_abun_scaled = [mpeak_hts[2] / target_abun[2], mpeak_hts[2] * (target_abun[1] / target_abun[2]), mpeak_hts[2]]
    else:
        # no matching peaks found
        return None, None
    # (debug) print out peak params, plot fit compared to isotopes
    for i in [0, 1, 2]:
        m_name = 'M  ' if i == 0 else 'M+{}'.format(i)
        if mpeak_mzs[i] is not None:
            mz = mpeak_mzs[i]
            dmz = mz - target_mz[i]
            ht = mpeak_hts[i]
            dht = ht - target_abun_scaled[i]
            msg = '{} m/z: {:9.4f} ({:+6.4f}) intensity: {:.2e} ({:+.2e})'.format(m_name, mz, dmz, ht, dht)
            _debug_handler(debug_flag, debug_cb, msg=msg)
        else:
            _debug_handler(debug_flag, debug_cb, msg=m_name + ' peak not found')
    isotope_dist_img = _plot_ms1_isotope_dist(ms1_i_mz, ms1_i_int, target_mz, target_abun_scaled, mpeak_mzs, mpeak_hts, mpeak_wts)
    mz_score, int_score = _mz_int_isotope_score(target_mz, target_abun_scaled, mpeak_mzs, mpeak_hts)
    cos_score = _calc_cosine_from_recon_spectrum(ms1_i_mz, ms1_i_int, target_mz, target_abun_scaled)
    # (debug)
    msg = 'mz ppm: {:.3f} abundance percent: {:.3f} cos distance: {:.3f}'.format(mz_score, int_score, cos_score)
    _debug_handler(debug_flag, debug_cb, msg=msg, img=isotope_dist_img)
    return (mz_score, int_score, cos_score), isotope_dist_img


def _do_fragment_scoring(oz_data, fragment, formula, mz_tol, ms1_fit_method, rt_bounds, rt_fit_method, precursor_rt, 
                         rt_tol, pre_xic_rt, pre_xic_int, check_saturation, saturation_threshold, 
                         debug_flag, debug_cb):
    """
    ?

    Parameters
    ----------

    debug_flag : ``str``
        specifies how to dispatch the message and/or plot, None to do nothing
    debug_cb : ``func``
        callback function that takes the debugging message as an argument, can be None if
        debug_flag is not set to 'textcb'

    Returns
    -------

    """
    fragment_results = {}
    # scoring for aldehyde
    _debug_handler(debug_flag, debug_cb, msg=fragment)
    target_mz, target_abun = predict_m_m1_m2(formula)
    # get XIC for  target_rt +/- 0.5 using aldehyde m/z +/- mz_tol
    mz_min, mz_max = target_mz[0] - mz_tol, target_mz[0] + mz_tol
    xic_rt, xic_int = oz_data.collect_xic_arrays_by_mz(mz_min, mz_max, rt_bounds=rt_bounds)
    # get the fitted peak info from the xic
    peak_info, xic_fit_img = _fit_xic_rt(xic_rt, xic_int, target_mz[0], mz_tol, 
                                         precursor_rt, rt_tol, rt_fit_method, 1., 
                                         debug_flag, debug_cb)
    peak_rt, peak_ht, peak_wt = peak_info
    if peak_rt is not None:
        # TODO (Dylan Ross) make this check more sophisticated, detect peak saturation from peak shape or something
        is_saturated = peak_ht > saturation_threshold
        if check_saturation and is_saturated:
            msg = 'peak saturation detected, only using scans from leading edge of RT peak for isotope scoring'
            _debug_handler(debug_flag, debug_cb, msg=msg)
            # find the RT at the leading edge of the peak (RT scans from between 5% and 25% peak height)
            rt_min, rt_max = _find_leading_edge_rt(xic_rt, xic_int, peak_rt, peak_ht, peak_wt, (0.05, 0.25), 
                                                   debug_flag, debug_cb)
        else:
            rt_min, rt_max = peak_rt - rt_tol, peak_rt + rt_tol
        # check for case where saturation correction fails, use original RT bounds
        if rt_min is None or rt_max is None:
            msg = 'failed to select scans from leading edge of RT peak, reverting to original RT bounds'
            _debug_handler(debug_flag, debug_cb, msg=msg)
            # unset is_saturated
            is_saturated = False
            rt_min, rt_max = peak_rt - rt_tol, peak_rt + rt_tol
        scores, iso_dist_img = _extract_ms1_and_calc_isotope_score(oz_data, target_mz, target_abun, 
                                                                   rt_min, rt_max, ms1_fit_method, 
                                                                   debug_flag, debug_cb)
        if scores is not None:
            fragment_results = {
                'target_mz': target_mz[0],
                'target_rt': precursor_rt,
                'xic_peak_rt': peak_rt,
                'xic_peak_ht': peak_ht,
                'xic_peak_fwhm': peak_wt,
                'mz_ppm': scores[0],
                'abun_percent': scores[1],
                'mz_cos_dist': scores[2],
                'rt_cos_dist': _calc_cosine_from_xic(pre_xic_rt, pre_xic_int, xic_rt, xic_int),
                'isotope_dist_img': iso_dist_img,
                'xic_fit_img': xic_fit_img,
                'saturation_corrected': check_saturation and is_saturated,
            }
            return fragment_results
    else:
        _debug_handler(debug_flag, debug_cb, msg='unable to find RT peak for fragment')
        return None


def _do_fragment_scoring_infusion(oz_data, fragment, formula, mz_tol, ms1_fit_method, debug):
    fragment_results = {}
    # scoring for aldehyde
    if debug:
        print(fragment)
    target_mz, target_abun = predict_m_m1_m2(formula)
    scores, iso_dist_img = _extract_ms1_and_calc_isotope_score(oz_data, target_mz, target_abun, 
                                                               oz_data.min_rt, oz_data.max_rt, ms1_fit_method, debug)
    if scores is not None:
        fragment_results = {
            'target_mz': target_mz[0],
            'mz_ppm': scores[0],
            'abun_percent': scores[1],
            'mz_cos_dist': scores[2],
            'isotope_dist_img': iso_dist_img,
        }
        return fragment_results


def score_db_pos_isotope_dist_polyunsat(oz_data, precursor_formula, fa_nc, fa_nu, precursor_rt, rt_tol, 
                                        rt_peak_win, mz_tol, rt_fit_method='gauss', ms1_fit_method='localmax', 
                                        check_saturation=True, saturation_threshold=1e5, remove_d=None, 
                                        debug_flag=None, debug_cb=None, info_cb=None):
    """
    performs isotope distribution scoring for a range of potential double-bond positions for polyunsaturated lipids, 
    also works for monounsaturated lipids, and essentially does nothing for completely saturated lipids

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
    precursor_rt : ``float``
        precursor retention time
    rt_tol : ``float``
        retention time tolerance
    rt_peak_win : ``float``
        size of RT window to extract for peak fitting
    mz_tol : ``float``
        m/z tolerance for extracting XICs
    rt_fit_method : ``str``, default='gauss'
        specify method to use for fitting the RT peak ('gauss' works best in testing)
    ms1_fit_method : ``str``, default='localmax'
        specify method to use for fitting MS1 spectrum ('localmax' works best in testing)
    check_saturation : ``bool``, default=True
        whether to check for signal saturation and use leading edge strategy if necessary
    saturation_threshold : ``float``, default=1e5
        specify a threshold intensity for determining peak saturation
    remove_d : ``int``, optional
        adjust molecular formulas to get rid of D labels on fatty acid tail that are part of the neutral loss
        (specific to SPLASH lipids)
    debug_flag : ``str``, optional
        specifies how to dispatch the message and/or plot, None to do nothing
    debug_cb : ``func``, optional
        callback function that takes the debugging message as an argument, can be None if
        debug_flag is not set to 'textcb'
    info_cb : ``function``, optional
        optional callback function that gets called at several intermediate steps and gives information about data
        processing details. Callback function takes a single argument which is a ``str`` info message

    Returns
    -------
    result : dict(...)
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
    rtb = (precursor_rt - rt_peak_win, precursor_rt + rt_peak_win)
    pre_xic_rt, pre_xic_int = oz_data.collect_xic_arrays_by_mz(mz_min, mz_max, rt_bounds=rtb)
    # get the fitted peak info from the xic
    peak_info, pre_xic_fit_img = _fit_xic_rt(pre_xic_rt, pre_xic_int, pre_target_mz[0], mz_tol, 
                                             precursor_rt, rt_tol, rt_fit_method, 1., \
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
        'target_rt': precursor_rt,
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
    n_combos = len(dbidx_dbpos)
    if info_cb is not None:
        msg = 'INFO: there are {} combinations of (db_idx, db_pos) to examine'
        info_cb(msg.format(n_combos))
    i = 1
    for db_idx, db_pos in dbidx_dbpos:
        # ITERATE THROUGH DB POSITIONS AND LOOK FOR FRAGMENTS
        #---------------------------------------------------------------------------------------------------
        if db_idx not in results['fragments']:
            results['fragments'][db_idx] = {}
        results['fragments'][db_idx][db_pos] = {}
        msg = '\n----------db_idx={},db_pos={}----------'.format(db_idx, db_pos)
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
                                           rt_fit_method, precursor_rt, rt_tol, pre_xic_rt, pre_xic_int, 
                                           check_saturation, saturation_threshold, 
                                           debug_flag, debug_cb)
        # scoring for criegee
        crg_results = _do_fragment_scoring(oz_data, 'criegee', crg_formula, mz_tol, ms1_fit_method, rtb, 
                                           rt_fit_method, precursor_rt, rt_tol, pre_xic_rt, pre_xic_int, 
                                           check_saturation, saturation_threshold, 
                                           debug_flag, debug_cb)
        results['fragments'][db_idx][db_pos]['aldehyde'] = ald_results
        results['fragments'][db_idx][db_pos]['criegee'] = crg_results
        if info_cb is not None:
            msg = 'INFO: {} of {} done'.format(i, n_combos)
            info_cb(msg)
            i += 1
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

