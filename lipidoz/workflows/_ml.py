"""
lipidoz/workflows/_ml.py
Dylan Ross (dylan.ross@pnnl.gov)

    define components of standard high-level workflows for OzID
"""


import os
from tempfile import TemporaryDirectory

import numpy as np
from mzapy import MZA
from mzapy.isotopes import valid_ms_adduct, ms_adduct_formula, monoiso_mass

from lipidoz import __version__ as VER
from lipidoz._util import (
    CustomUReader,
    _debug_handler, 
    _polyunsat_ald_crg_formula, 
    _calc_dbp_bounds,
    new_lipidoz_results
)
from lipidoz.ml.data import (
    load_preml_data, 
    load_ml_targets, 
    split_true_and_false_preml_data, 
    preml_to_ml_data
)
from lipidoz.workflows._isotope_scoring import run_isotope_scoring_workflow_targeted
from lipidoz.ml.models.resnet18 import ResNet18
from lipidoz._pyliquid import parse_lipid_name


def _load_target_list(target_list_file, ignore_preferred_ionization, rt_correction_func):
    """
    loads target list from .csv file, performs some basic validation
    Exceptions raised in this function happen before MZA is initialized which is better than having
    them raised while it is initialized (the IO threads can complicate things), so as much parameter
    validation should be done here as possible

    Parameters
    ----------
    target_list_file : ``str``
        filename and path for target list (.csv format)
    ignore_preferred_ionization : ``bool``
        whether to ignore cases where a lipid/adduct combination violates the lipid class' preferred ionization state
    rt_correction_func : ``function``
        if not None, specifies a function that takes uncorrected retention time as an argument and returns the 
        corrected retention time

    Returns
    -------
    target_lipids : ``list(pyliquid.lipids.LipidWithChains)``
        target lipids as ``pyliquid.lipids.LipidWithChains``
    target_adducts : ``list(str)``
        MS adducts of target lipids
    target_rts : ``list(float)``
        retention times of target lipids (corrected if rt_correction_func is not None)
    """
    names, adducts, rts = np.loadtxt(target_list_file, dtype=str, delimiter=',', skiprows=1, 
                                     unpack=True, comments='#', ndmin=2)
    target_lipids = [parse_lipid_name(name) for name in names]
    adducts = [adduct for adduct in adducts]
    rts = [float(rt) for rt in rts]
    # perform as much input validation as possible ahead of time
    for i, (name, lipid, adduct) in enumerate(zip(names, target_lipids, adducts)):
        line = i + 2
        if lipid is None:
            msg = '_load_target_list: name "{}" was not able to be parsed as lipid (line: {})'
            raise ValueError(msg.format(name, line))
        if lipid.__class__.__name__ != 'LipidWithChains':
            msg = ('_load_target_list: lipid {} has more than 1 acyl chain ({}) but individual chain compositions'
                   ' were not provided (line: {})')
            raise ValueError(msg.format(name, lipid.n_chains, line))
        if not valid_ms_adduct(adduct):
            msg = '_load_target_list: {} is not a recognized MS adduct (line: {})'
            raise ValueError(msg.format(adduct, line))
        if ((adduct[-1] == '+' and lipid.ionization not in ['both', 'pos']) 
            or (adduct[-1] == '-' and lipid.ionization not in ['both', 'neg'])):
            if not ignore_preferred_ionization:
                msg = '_load_target_list: lipid {} ionization is {} but adduct was {} (line: {})'
                raise ValueError(msg.format(name, lipid.ionization, adduct, line))
    # correct RT values if a correction function was provided
    if rt_correction_func is not None:
        rts = [rt_correction_func(_) for _ in rts]
    return target_lipids, adducts, rts


def collect_preml_dataset(oz_data_file, target_list_file, rt_tol, d_label=None, d_label_in_nl=None, 
                          debug_flag=None, debug_cb=None, ignore_preferred_ionization=False, 
                          rt_correction_func=None, mza_version='new'):
    """
    collects a dataset which can be used in training ML models. The dataset is a dictionary with metadata
    and minimally processed RTMZ data. The RTMZ data is extracted in a window with the following bounds:

    * target RT +/- rt_tol -- this should be set wide enough to accomodate the chromatographic peak
    * target m/z (M isotope) - 0.5, target m/z (M isotope) + 2.5 -- this covers the M, M+1, M+2 isotopes 

    Parameters
    ----------
    oz_data_file : ``str``
        filename and path for OzID data (.mza format)
    target_list_file : ``str``
        filename and path for target list (.csv format)
    rt_tol : ``float``
        retention time tolerance, defines data extraction window
    d_label : ``int``, optional
        number of deuteriums in deuterium-labeled standards (*i.e.* SPLASH and Ultimate SPLASH mixes)
    d_label_in_nl : ``bool``, optional
        if deuterium labels are present, indicates whether they are included in the neutral loss during OzID 
        fragmentation, this is False if the deuteriums are on the lipid head group and True if they are at the end of 
        a FA tail (meaning that the aldehyde and criegee fragment formulas must be adjusted to account for loss of the
        label during fragmentataion)
    debug_flag : ``str``, optional
        specifies how to dispatch the message and/or plot, None to do nothing
    debug_cb : ``func``, optional
        callback function that takes the debugging message as an argument, can be None if
        debug_flag is not set to 'textcb'
    ignore_preferred_ionization : ``bool``, default=False
        whether to ignore cases where a lipid/adduct combination violates the lipid class' preferred ionization state
    rt_correction_func : ``function``, optional
        provide a function that takes an uncorrected retention time as an argument
        then returns the corrected retention time
    mza_version : ``str``, default='new'
            temporary measure for indicating whether the the scan indexing needs to account for partitioned
            scan data ('new') or not ('old'). Again, this is only temporary as at some point the mza version
            will be encoded as metadata into the file and this accommodation can be made automatically.

    Returns
    -------
    pre_ml_dataset : ``dict(...)``
        dataset used for assembling ML training data
    """
    # load the data 
    match (ext := os.path.splitext(oz_data_file)[-1]):
        case ".uimf":
            oz_data = CustomUReader(oz_data_file)
            # additional initialization, skip Frame 1 
            oz_data.accum_frame_spectra_allscans(skip_frame_1=True)
        case ".mza":
            oz_data = MZA(oz_data_file, cache_scan_data=True, mza_version="new")
        case _:
            raise ValueError(f"unrecognized raw data file extension: {ext}")
    # store metadata
    pre_ml_dset = {
        'metadata': {
            'workflow': 'pre_ml', 
            'lipidoz_version': VER, 
            'oz_data_file': oz_data_file,
            'target_list_file': target_list_file, 
            'rt_tol': rt_tol, 
            'd_label': d_label, 
            'd_label_in_nl': d_label_in_nl,
        },
        'targets': {},
    }
    # load the target list
    target_lipids, target_adducts, target_rts = _load_target_list(target_list_file, 
                                                                  ignore_preferred_ionization, 
                                                                  rt_correction_func)
    # main data processing
    for tlipid, tadduct, trt in zip(target_lipids, target_adducts, target_rts):
        if d_label is not None:
            # adjust formula by replacing protons with deuterium label
            tlipid.formula['H'] -= d_label
            tlipid.formula['D'] = d_label
        remove_d = d_label if ((d_label is not None) and d_label_in_nl) else None
        adduct_formula = ms_adduct_formula(tlipid.formula, tadduct)
        # run the analysis for individual lipid species
        msg = '\n' + str(tlipid) + ' ' + tadduct + '\n=========================='
        _debug_handler(debug_flag, debug_cb, msg=msg)
        # pre
        rt_min, rt_max = trt - rt_tol, trt + rt_tol
        mz = monoiso_mass(adduct_formula)
        pre_data = oz_data.collect_rtmz_arrays(rt_bounds=(rt_min, rt_max), mz_bounds=(mz - 1.5, mz + 2.5))
        # compute bounds of db indices and positions, keep only unique pairs
        dbidx_dbpos = set()
        for fac, fau in zip(tlipid.fa_carbon_chains, tlipid.fa_unsat_chains):
            if fac > 0:
                for db_idx in range(1, fau + 1):
                    dbp_min, dbp_max = _calc_dbp_bounds(fac, fau, db_idx)
                    for db_pos in range(dbp_min, dbp_max + 1):
                        dbidx_dbpos.add((db_idx, db_pos))
        # iterate over the unique dbidx/dbpos combinations
        for db_idx, db_pos in dbidx_dbpos:
            k = '{}|{}|{:.2f}min|{}|{}'.format(str(tlipid), tadduct, trt, db_idx, db_pos)
            # get aldehyde and criegee formulas
            ald_formula, crg_formula = _polyunsat_ald_crg_formula(adduct_formula, db_pos, db_idx)
            if remove_d is not None:
                # adjust the molecular formulas to get rid of D7 labels in neutral loss fragments
                # applicable to SPLASH lipid standards only
                ald_formula.pop('D')
                ald_formula['H'] += remove_d
                crg_formula.pop('D')
                crg_formula['H'] += remove_d
            # ald
            ald_mz = monoiso_mass(ald_formula)
            ald_data = oz_data.collect_rtmz_arrays(rt_bounds=(rt_min, rt_max), 
                                                    mz_bounds=(ald_mz - 1.5, ald_mz + 2.5))
            # crg
            crg_mz = monoiso_mass(crg_formula)
            crg_data = oz_data.collect_rtmz_arrays(rt_bounds=(rt_min, rt_max), 
                                                    mz_bounds=(crg_mz - 1.5, crg_mz + 2.5))
            # add data to dataset
            pre_ml_dset['targets'][k] = {
                'pre_data': pre_data, 
                'ald_data': ald_data, 
                'crg_data': crg_data,
                'pre_mz': mz,
                'ald_mz': ald_mz,
                'crg_mz': crg_mz,
                'rt': trt,
            }
            _debug_handler(debug_flag, debug_cb, msg=k)
    # clean up
    oz_data.close()
    return pre_ml_dset


def convert_multi_preml_datasets_labeled(preml_files, ml_target_files, 
                                         rt_sampling_augment=True, normalize_intensity=True, rt_corr_funcs=None,
                                         debug_flag=None, debug_cb=None):
    """
    iterates through pairs of pre-ml dataset files and ml target lists, splits the pre-ml datasets into True/False
    examples, then converts to binned ml datasets. True/False examples from all pairs are combined and returned as 
    two arrays of ML data

    Parameters
    ----------
    preml_files : ``list(str)``
        paths to pre-ml dataset files
    ml_target_files : ``list(str)``
        paths to target list .csv files
    rt_sampling_augment : ``bool``, default=True
        re-sample RT dimension from RTMZ data multiple times in order to augment training examples (~10x)
    normalize_intensity : ``bool``, default=True
        normalize the intensities in each 2D RTMZ array so that they are in the range 0->1
    rt_corr_funcs : ``list(function)``, optional
        list of RT correction functions, one per ml_target_file, applies retention time corrections 
    debug_flag : ``str``, optional
        specifies how to dispatch the message and/or plot, None to do nothing
    debug_cb : ``func``, optional
        callback function that takes the debugging message as an argument, can be None if
        debug_flag is not set to 'textcb'

    Returns
    -------
    true_ml_data : ``numpy.ndarray``
    false_ml_data : ``numpy.ndarray``
        arrays of binned data for ML split by annotation from all pairs with shapes: (N, 3, 24, 400), where N is the 
        number of True or False examples in the array
    """
    i = 0
    t_ml_data, f_ml_data = [], []
    rt_corr_funcs = [None for _ in preml_files] if rt_corr_funcs is None else rt_corr_funcs
    for preml_file, ml_target_file, rt_corr_func in zip(preml_files, ml_target_files, rt_corr_funcs):
        msg = 'converting pre-ml data from file: {} and target list: {} ...'.format(preml_file, ml_target_file)
        _debug_handler(debug_flag, debug_cb, msg=msg)
        t_preml_data, f_preml_data = split_true_and_false_preml_data(load_preml_data(preml_file), 
                                                                     load_ml_targets(ml_target_file, 
                                                                                     rt_corr_func=rt_corr_func))
        msg = '  True targets: {} False targets: {}'.format(len(t_preml_data['targets']), 
                                                            len(f_preml_data['targets']))
        _debug_handler(debug_flag, debug_cb, msg=msg)
        t_ml_data.append(preml_to_ml_data(t_preml_data, 
                                          rt_sampling_augment=rt_sampling_augment, 
                                          normalize_intensity=normalize_intensity,
                                          debug_flag=debug_flag, debug_cb=debug_cb))
        f_ml_data.append(preml_to_ml_data(f_preml_data, 
                                          rt_sampling_augment=rt_sampling_augment, 
                                          normalize_intensity=normalize_intensity,
                                          debug_flag=debug_flag, debug_cb=debug_cb))
        
        msg = ('... finished converting pre-ml data from file: {} and target list: {} for {} True and {} False'
                ' training examples')
        n_T, n_F = len(t_ml_data[i]), len(f_ml_data[i])
        i += 1
        msg = msg.format(preml_file, ml_target_file, n_T, n_F)
        _debug_handler(debug_flag, debug_cb, msg=msg)
    return np.concatenate(t_ml_data), np.concatenate(f_ml_data)


def convert_multi_preml_datasets_unlabeled(preml_files, normalize_intensity=True,
                                           debug_flag=None, debug_cb=None):
    """
    iterates through pre-ml dataset files converts to binned ml datasets (for unlabeled data)

    Parameters
    ----------
    preml_files : ``list(str)``
        paths to pre-ml dataset files
    normalize_intensity : ``bool``, default=True
        normalize the intensities in each 2D RTMZ array so that they are in the range 0->1
    debug_flag : ``str``, optional
        specifies how to dispatch the message and/or plot, None to do nothing
    debug_cb : ``func``, optional
        callback function that takes the debugging message as an argument, can be None if
        debug_flag is not set to 'textcb'

    Returns
    -------
    ml_data : ``numpy.ndarray``
        array of binned data for ML with shape: (N, 3, 24, 400), where N is the number of training examples
    """
    i = 0
    ml_data = []
    for preml_file in preml_files:
        msg = 'converting pre-ml data from file: {} ...'.format(preml_file)
        _debug_handler(debug_flag, debug_cb, msg=msg)
        preml_data = load_preml_data(preml_file)
        ml_data.append(preml_to_ml_data(preml_data,
                                        normalize_intensity=normalize_intensity,
                                        debug_flag=debug_flag, debug_cb=debug_cb))
        msg = '... finished converting pre-ml data from file: {} for {} training examples'
        n = len(ml_data[i])
        _debug_handler(debug_flag, debug_cb, msg=msg.format(preml_file, n))
        i += 1
    return np.concatenate(ml_data)


def _convert_preml_data_and_dl_labels_into_targeted_target_list(lipidoz_results):
    """
    use information from the preml_data and DL predicted labels in a lipidoz_results 
    dict to construct a target list that will work with the targeted variant of
    the isotope scoring workflow

    Parameters
    ----------
    lipidoz_results : ``dict(...)``

    Returns
    -------
    target_list_content : ``str``
        contents of the generated target list (.csv) as a string
    """
    target_list_content = "lipid,adduct,retention_time,db_idx,db_pos\n"
    for k, l in zip(lipidoz_results["preml_data"]["targets"].keys(), lipidoz_results["ml_pred_lbls"]):
        if int(l) == 1:
            target_list_content += k.replace("|", ",").replace("min", "") + "\n"
    return target_list_content


def hybrid_deep_learning_and_isotope_scoring(oz_data_file, target_list_file, rt_tol, rt_peak_win, mz_tol,
                                             dl_params_file,
                                             d_label=None, d_label_in_nl=None, 
                                             debug_flag=None, debug_cb=None):
    """
    A hybrid workflow that incorporates deep learning inference as a prefilter
    then performs targeted isotope scoring workflow on predicted True double bond positions

    Parameters
    ----------
    oz_data_file : ``str``
        filename and path for OzID data (.mza format)
    target_list_file : ``str``
        filename and path for target list (.csv format)
    rt_tol : ``float``
        retention time tolerance, defines data extraction window
    rt_peak_win : ``float``
        size of retention time window to extract for fitting retention time peak
    mz_tol : ``float``
        m/z tolerance for extracting XICs
    dl_params_file : ``str``
        pre-trained DL model parameters file
    d_label : ``int``, optional
        number of deuteriums in deuterium-labeled standards (*i.e.* SPLASH and Ultimate SPLASH mixes)
    d_label_in_nl : ``bool``, optional
        if deuterium labels are present, indicates whether they are included in the neutral loss during OzID 
        fragmentation, this is False if the deuteriums are on the lipid head group and True if they are at the end of 
        a FA tail (meaning that the aldehyde and criegee fragment formulas must be adjusted to account for loss of the
        label during fragmentataion)
    debug_flag : ``str``, optional
        specifies how to dispatch the message and/or plot, None to do nothing
    debug_cb : ``func``, optional
        callback function that takes the debugging message as an argument, can be None if
        debug_flag is not set to 'textcb'
    
    Returns
    -------
    lipidoz_results : ``dict(...)``
        results from DL prefiltering and targeted isotope scoring analysis
    """
    # make the lipidoz results dict to store everything in
    lipidoz_results = new_lipidoz_results()
    # extract pre-ml dataset and store in lipidoz_results
    lipidoz_results["preml_data"] = collect_preml_dataset(oz_data_file, target_list_file, rt_tol, 
                                                          d_label=d_label, d_label_in_nl=d_label_in_nl, 
                                                          debug_flag=debug_flag, debug_cb=debug_cb)
    # convert pre-ml data to ML data and store in lipidoz_results
    lipidoz_results["ml_data"] = preml_to_ml_data(lipidoz_results["preml_data"], 
                                                  debug_flag=debug_flag, 
                                                  debug_cb=debug_cb)
    # load the model and its pre-trained parameters
    rn18 = ResNet18()
    rn18.load(dl_params_file)
    # run inference and store predictions
    lipidoz_results["ml_pred_lbls"] = rn18.predict(lipidoz_results["ml_data"])
    lipidoz_results["ml_pred_probs"] = rn18.predict_proba(lipidoz_results["ml_data"])
    # TODO (Dylan Ross): Use the predicted True double bond positions as targets for targeted
    #                    variant of the isotope distribution analysis. For the initial implementation
    #                    just go ahead and create a tempfile in csv format to serve as the input 
    #                    target list for the isotope distribution analysis and run that as normal. The
    #                    problem with this approach though is that raw data extraction will need to be
    #                    repeated for the targets in the isotope distribution analysis while all info
    #                    that is needed should already technically be contained in the pre-ml data. 
    #                    I doubt this will be too big of a deal since the targeted variant of the
    #                    isotope distribution analysis should already cut down on processing time due
    #                    to the greatly reduced amount of data that needs to be considered in detail,
    #                    but in the future it would probably make sense to go ahead and come up with
    #                    another set of isotope distribution function variants that work directly from
    #                    already extracted pre-ml data rather than doing the data extraction from the 
    #                    MZA files themselves as they do now.
    target_list_content = _convert_preml_data_and_dl_labels_into_targeted_target_list(lipidoz_results)
    # run targeted isotope distribution analysis 
    # (use temporary directory with generated target list csv)
    with TemporaryDirectory() as tmp_dir:
        target_list_file = os.path.join(tmp_dir, "filtered_targets.csv")
        with open(target_list_file, "w") as outf:
            outf.write(target_list_content)
        lipidoz_results["isotope_scoring_results"] = run_isotope_scoring_workflow_targeted(oz_data_file, 
                                                                                           target_list_file, 
                                                                                           rt_tol, 
                                                                                           rt_peak_win, 
                                                                                           mz_tol, 
                                                                                           d_label=d_label, 
                                                                                           d_label_in_nl=d_label_in_nl, 
                                                                                           debug_flag=debug_flag, 
                                                                                           debug_cb=debug_cb)
    return lipidoz_results
    
