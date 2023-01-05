"""
lipidoz/workflows.py

Dylan Ross (dylan.ross@pnnl.gov)

    define components of standard high-level workflows for OzID

    TODO (Dylan Ross) add in support for isotope analysis on non-chromatographic (i.e. direct infusion) data
"""


import os
import pickle

import numpy as np
import xlsxwriter

from mzapy import MZA
from mzapy.isotopes import valid_ms_adduct, monoiso_mass, ms_adduct_formula

from lipidoz import __version__ as VER
from lipidoz.isotope_scoring import score_db_pos_isotope_dist_polyunsat, score_db_pos_isotope_dist_polyunsat_infusion
from lipidoz._util import _polyunsat_ald_crg_formula, _calc_dbp_bounds, _debug_handler
from lipidoz.ml.data import load_preml_data, load_ml_targets, split_true_and_false_preml_data, preml_to_ml_data
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


def _load_target_list_infusion(target_list_file, ignore_preferred_ionization):
    """
    similar to _load_target_list function but for infusion data (no retention time)

    Parameters
    ----------
    target_list_file : ``str``
        filename and path for target list (.csv format)
    ignore_preferred_ionization : ``bool``
        whether to ignore cases where a lipid/adduct combination violates the lipid class' preferred ionization state

    Returns
    -------
    target_lipids : ``list(pyliquid.lipids.LipidWithChains)``
        target lipids as ``pyliquid.lipids.LipidWithChains``
    target_adducts : ``list(str)``
        MS adducts of target lipids
    """
    names, adducts = np.loadtxt(target_list_file, dtype=str, delimiter=',', skiprows=1, 
                                unpack=True, comments='#', ndmin=2)
    target_lipids = [parse_lipid_name(name) for name in names]
    adducts = [adduct for adduct in adducts]
    # perform as much input validation as possible ahead of time
    for i, (name, lipid, adduct) in enumerate(zip(names, target_lipids, adducts)):
        line = i + 2
        if lipid is None:
            msg = '_load_target_list_infusion: name "{}" was not able to be parsed as lipid (line: {})'
            raise ValueError(msg.format(name, line))
        if lipid.__class__.__name__ != 'LipidWithChains':
            msg = ('_load_target_list_infusion: lipid {} has more than 1 acyl chain ({}) but individual chain compositions'
                   ' were not provided (line: {})')
            raise ValueError(msg.format(name, lipid.n_chains, line))
        if not valid_ms_adduct(adduct):
            msg = '_load_target_list_infusion: {} is not a recognized MS adduct (line: {})'
            raise ValueError(msg.format(adduct, line))
        if ((adduct[-1] == '+' and lipid.ionization not in ['both', 'pos']) 
            or (adduct[-1] == '-' and lipid.ionization not in ['both', 'neg'])):
            if not ignore_preferred_ionization:
                msg = '_load_target_list_infusion: lipid {} ionization is {} but adduct was {} (line: {})'
                raise ValueError(msg.format(name, lipid.ionization, adduct, line))
    return target_lipids, adducts


def run_isotope_scoring_workflow(oz_data_file, target_list_file, rt_tol, rt_peak_win, mz_tol, 
                                 d_label=None, d_label_in_nl=None, progress_cb=None, info_cb=None, 
                                 early_stop_event=None, debug_flag=None, debug_cb=None, rt_correction_func=None, 
                                 ignore_preferred_ionization=False, mza_version='new'):
    """
    workflow for performing isotope scoring for the determination of db positions. 
    inputs are the data file and target list file, output is a dictionary containing metadata about the analysis and
    the analysis results for all of the lipids in the target list. 
    The target list should have columns containing the following information (in order, 1 header row is skipped):
    
    * lipid name in standard abbreviated format, with FA composition fully specified,
      *e.g.*, PC(18:1_16:0) or TG(16:0/18:1/20:2)
    * MS adduct, *e.g.*, [M+H]+ or [M-2H]2-
    * target retention time

    Parameters
    ----------
    oz_data_file : ``str``
        filename and path for OzID data (.mza format)
    target_list_file : ``str``
        filename and path for target list (.csv format)
    rt_tol : ``float``
        retention time tolerance (for MS1 data extraction)
    rt_peak_win : ``float``
        size of retention time window to extract for fitting retention time peak
    mz_tol : ``float``
        m/z tolerance for extracting XICs
    d_label : ``int``, optional
        number of deuteriums in deuterium-labeled standards (*i.e.* SPLASH and Ultimate SPLASH mixes)
    d_label_in_nl : ``bool``, optional
        if deuterium labels are present, indicates whether they are included in the neutral loss during OzID 
        fragmentation, this is False if the deuteriums are on the lipid head group and True if they are at the end of 
        a FA tail (meaning that the aldehyde and criegee fragment formulas must be adjusted to account for loss of the
        label during fragmentataion)
    progress_cb : ``function``, optional
        option for a callback function that gets called every time an individual lipid species has been processed, this
        callback function should take as arguments (in order):

        * lipid name (``str``)
        * adduct (``str``)
        * current position in target list (``int``)
        * total lipids in target list(``int``)

    info_cb : ``function``, optional
        optional callback function that gets called at several intermediate steps and gives information about data
        processing details. Callback function takes a single argument which is a ``str`` info message
    early_stop_event : ``threading.Event``, optional
        When the workflow is running in its own thread and this event gets set, processing is stopped gracefully
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
    isotope_scoring_results : ``dict(...)``
        results dictionary with metadata and scoring information
    """
    # store metadata
    results = {
        'metadata': {
            'workflow': 'isotope_scoring',
            'lipidoz_version': VER,
            'oz_data_file': oz_data_file,
            'target_list_file': target_list_file, 
            'rt_tol': rt_tol, 
            'rt_peak_win': rt_peak_win, 
            'mz_tol': mz_tol,
            'd_label': d_label, 
            'd_label_in_nl': d_label_in_nl,
        },
        'targets': {},
    }
    # load the target list
    target_lipids, target_adducts, target_rts = _load_target_list(target_list_file, 
                                                                  ignore_preferred_ionization, 
                                                                  rt_correction_func)
    n = len(target_lipids)
    if info_cb is not None:
        msg = 'INFO: loaded target list: {} ({} targets)'.format(target_list_file, n)
    # load the data 
    oz_data = MZA(oz_data_file, cache_scan_data=True, mza_version=mza_version)
    if info_cb is not None:
        msg = 'INFO: loaded OzID data file: {}'.format(oz_data_file)
        info_cb(msg)
    # main data processing
    i = 1
    for tlipid, tadduct, trt in zip(target_lipids, target_adducts, target_rts):
        # check for a stop event
        if early_stop_event is not None and early_stop_event.is_set():
            if info_cb is not None and debug_flag is None:  # avoid duplication of messages
                info_cb('INFO: EARLY STOP EVENT HAS BEEN SET!')
            _debug_handler(debug_flag, debug_cb, msg='early stop event has been set')
            oz_data.close()
            return None
        if info_cb is not None:
            msg = 'INFO: analyzing target lipid {} {}'.format(tlipid, tadduct)
            info_cb(msg)
        if d_label is not None:
            # adjust formula by replacing protons with deuterium label
            tlipid.formula['H'] -= d_label
            tlipid.formula['D'] = d_label
        # check for SPLASH lipids
        remove_d = d_label if ((d_label is not None) and d_label_in_nl) else None
        adduct_formula = ms_adduct_formula(tlipid.formula, tadduct)
        # run the analysis for individual lipid species
        msg = '\n' + str(tlipid) + ' ' + tadduct + '\n=========================='
        _debug_handler(debug_flag, debug_cb, msg=msg)
        lipid_result = score_db_pos_isotope_dist_polyunsat(oz_data, adduct_formula, tlipid.fa_carbon_chains, 
                                                           tlipid.fa_unsat_chains, trt, rt_tol, 
                                                           rt_peak_win, mz_tol, remove_d=remove_d, 
                                                           debug_flag=debug_flag, debug_cb=debug_cb, info_cb=info_cb)
        # add individual result to full results
        rt_str = '{:.2f}min'.format(trt)
        if str(tlipid) in results['targets']:
            if tadduct in results['targets'][str(tlipid)]:
                results['targets'][str(tlipid)][tadduct][rt_str] = lipid_result
            else:
                results['targets'][str(tlipid)][tadduct] = {rt_str: lipid_result}
        else:
            results['targets'][str(tlipid)] = {tadduct: {rt_str: lipid_result}}
        # call progress callback function if provided
        if progress_cb is not None:
            progress_cb(str(tlipid), tadduct, i, n)
            i += 1
    # clean up
    oz_data.close()
    return results


def run_isotope_scoring_workflow_infusion(oz_data_file, target_list_file, mz_tol, 
                                          d_label=None, d_label_in_nl=None, progress_cb=None, early_stop_event=None, 
                                          debug_flag=None, debug_cb=None, 
                                          ignore_preferred_ionization=False, mza_version='new'):
    """
    workflow for performing isotope scoring for the determination of db positions from infusion data
    inputs are the data file and target list file, output is a dictionary containing metadata about the analysis and
    the analysis results for all of the lipids in the target list. 
    The target list should have columns containing the following information (in order, 1 header row is skipped):
    
    * lipid name in standard abbreviated format, with FA composition fully specified,
      *e.g.*, PC(18:1_16:0) or TG(16:0/18:1/20:2)
    * MS adduct, *e.g.*, [M+H]+ or [M-2H]2-

    Parameters
    ----------
    oz_data_file : ``str``
        filename and path for OzID data (.mza format)
    target_list_file : ``str``
        filename and path for target list (.csv format)
    mz_tol : ``float``
        m/z tolerance for extracting XICs
    d_label : ``int``, optional
        number of deuteriums in deuterium-labeled standards (*i.e.* SPLASH and Ultimate SPLASH mixes)
    d_label_in_nl : ``bool``, optional
        if deuterium labels are present, indicates whether they are included in the neutral loss during OzID 
        fragmentation, this is False if the deuteriums are on the lipid head group and True if they are at the end of 
        a FA tail (meaning that the aldehyde and criegee fragment formulas must be adjusted to account for loss of the
        label during fragmentataion)
    progress_cb : ``function``, optional
        option for a callback function that gets called every time an individual lipid species has been processed, this
        callback function should take as arguments (in order):

        * lipid name (``str``)
        * adduct (``str``)
        * current position in target list (``int``)
        * total lipids in target list(``int``)

    early_stop_event : ``threading.Event``, optional
        When the workflow is running in its own thread and this event gets set, processing is stopped gracefully
    debug_flag : ``str``, optional
        specifies how to dispatch the message and/or plot, None to do nothing
    debug_cb : ``func``, optional
        callback function that takes the debugging message as an argument, can be None if
        debug_flag is not set to 'textcb'
    ignore_preferred_ionization : ``bool``, default=False
        whether to ignore cases where a lipid/adduct combination violates the lipid class' preferred ionization state
    mza_version : ``str``, default='new'
            temporary measure for indicating whether the the scan indexing needs to account for partitioned
            scan data ('new') or not ('old'). Again, this is only temporary as at some point the mza version
            will be encoded as metadata into the file and this accommodation can be made automatically.

    Returns
    -------
    isotope_scoring_results : ``dict(...)``
        results dictionary with metadata and scoring information
    """
    # store metadata
    results = {
        'metadata': {
            'workflow': 'isotope_scoring_infusion',
            'lipidoz_version': VER,
            'oz_data_file': oz_data_file,
            'target_list_file': target_list_file, 
            'mz_tol': mz_tol,
            'd_label': d_label, 
            'd_label_in_nl': d_label_in_nl,
        },
        'targets': {},
    }
    # load the target list
    target_lipids, target_adducts = _load_target_list_infusion(target_list_file, ignore_preferred_ionization)
    # load the data 
    oz_data = MZA(oz_data_file, cache_scan_data=True, mza_version=mza_version)
    # main data processing
    n = len(target_lipids)
    i = 1
    for tlipid, tadduct in zip(target_lipids, target_adducts):
        # check for a stop event
        if early_stop_event is not None and early_stop_event.is_set():
            _debug_handler(debug_flag, debug_cb, msg='early stop event has been set')
            oz_data.close()
            return None
        if d_label is not None:
            # adjust formula by replacing protons with deuterium label
            tlipid.formula['H'] -= d_label
            tlipid.formula['D'] = d_label
        # check for SPLASH lipids
        remove_d = d_label if ((d_label is not None) and d_label_in_nl) else None
        adduct_formula = ms_adduct_formula(tlipid.formula, tadduct)
        # run the analysis for individual lipid species
        msg = '\n' + str(tlipid) + ' ' + tadduct + '\n=========================='
        _debug_handler(debug_flag, debug_cb, msg=msg)
        lipid_result = score_db_pos_isotope_dist_polyunsat_infusion(oz_data, adduct_formula, 
                                                                    tlipid.fa_carbon_chains, tlipid.fa_unsat_chains, 
                                                                    mz_tol, remove_d=remove_d, 
                                                                    debug_flag=debug_flag, debug_cb=debug_cb)
        # add individual result to full results
        rt_str = 'inf'
        if str(tlipid) in results['targets']:
            if tadduct in results['targets'][str(tlipid)]:
                results['targets'][str(tlipid)][tadduct][rt_str] = lipid_result
            else:
                results['targets'][str(tlipid)][tadduct] = {rt_str: lipid_result}
        else:
            results['targets'][str(tlipid)] = {tadduct: {rt_str: lipid_result}}
        # call progress callback function if provided
        if progress_cb is not None:
            progress_cb(str(tlipid), tadduct, i, n)
            i += 1
    # clean up
    oz_data.close()
    return results


def save_isotope_scoring_results(isotope_scoring_results, results_file_name, infusion=False):
    """
    save the results of the isotope scoring workflow (complete with metadata) to file in pickle format

    Parameters
    ----------
    isotope_scoring_results : ``dict(...)``
        results dictionary with metadata and scoring information
    results_file_name : ``str``
        filename and path to save the results file under, should have .loz file ending (maintains 
        compatibility with ``lipidoz_gui``)
    infusion : ``bool``, default=False
        indicates whether this is infusion data, if so the file extension should be .lozi
    """
    ext_should_be = '.lozi' if infusion else '.loz'
    ext = os.path.splitext(results_file_name)[-1]
    if ext != ext_should_be:
        msg = 'save_isotope_scoring_results: results file should have {} extension (was: {})'
        raise ValueError(msg.format(ext_should_be, ext))
    with open(results_file_name, 'wb') as pf:
        pickle.dump(isotope_scoring_results, pf)


def _write_metadata_to_sheet(metadata, workbook):
    sheet = workbook.get_worksheet_by_name('Metadata')
    fmt1 = workbook.add_format({'bold': True, 'align': 'right', 'font_size': 12})
    fmt2 = workbook.add_format({'align': 'left', 'font_size': 12})
    row = 0
    maxchar = 0
    for k, v in metadata.items():
        sheet.write(row, 0, k, fmt1)
        sheet.write(row, 1, v, fmt2)
        maxchar = max(len(str(v)), maxchar)
        row += 1    
    sheet.set_column(0, 0, 20)
    sheet.set_column(1, 1, maxchar)


def _write_targets_header_to_sheet(workbook):
    """
    21 columns
    """
    sheet = workbook.get_worksheet_by_name('Targets')
    headers = [
        'lipid', 'adduct', 'db_idx', 'db_pos', 
        'pre_mz', 'pre_rt', 'pre_xic_rt', 'pre_xic_ht', 'pre_xic_fwhm', 'pre_mz_cos_dist', 
        'ald_mz', 'ald_rt', 'ald_xic_rt', 'ald_xic_ht', 'ald_xic_fwhm', 'ald_mz_cos_dist', 'ald_rt_cos_dist', 
        'crg_mz', 'crg_rt', 'crg_xic_rt', 'crg_xic_ht', 'crg_xic_fwhm', 'crg_mz_cos_dist', 'crg_rt_cos_dist',
    ]
    fmt = workbook.add_format({'bold': True, 'font_size': 12})
    column_widths = [
        18, 8, 8, 8
    ] + [12 for _ in range(20)]
    for i in range(24):
        sheet.write(0, i, headers[i], fmt)
        sheet.set_column(i, i, column_widths[i])


def _add_line_to_targets_sheet(workbook, row, line_data):
    """
    24 columns
    """
    sheet = workbook.get_worksheet_by_name('Targets')
    formats = [workbook.add_format() for _ in range(24)]
    # .4f
    for i in [4, 9, 10, 15, 16, 17, 22, 23]:
        formats[i].set_num_format('0.0000')
    # .2f
    for i in [5, 6, 8, 11, 12, 14, 18, 19, 21]:
        formats[i].set_num_format('0.00')
    # d
    for i in [2, 3]:
        formats[i].set_num_format('0')
    # e
    for i in [7, 13, 19]:
        formats[i].set_num_format('0.00E+00')
    # apply conditional formatting - mz/rt cos dist
    for i in [9, 15, 16, 22, 23]:
        rule = {
            'type': '3_color_scale',
            'min_type': 'num',
            'mid_type': 'num',
            'max_type': 'num',
            'min_value': 0,
            'mid_value': 0.25,
            'max_value': 0.5,
            'min_color': '#AAAAFF',
            'mid_color': '#CCCCCC',
            'max_color': '#FFAAAA',
        }
        sheet.conditional_format(row, i, row, i, rule)
    # write the data
    for i in range(24):
        sheet.write(row, i, line_data[i], formats[i])


def write_isotope_scoring_report_xlsx(isotope_scoring_results, xlsx_file):
    """
    writes results of the isotope scoring workflow to an excel spreadsheet

    Parameters
    ----------
    isotope_scoring_results : ``dict(...)``
        results dictionary from isotope scoring workflow
    xlsx_file : ``str``
        filename to save report under
    """
    ext = os.path.splitext(xlsx_file)[-1]
    if ext != '.xlsx':
        msg = 'write_isotope_scoring_report_xlsx: report file should have .xlsx extension (was: "{}")'
        raise ValueError(msg.format(ext))
    # init workbook
    workbook = xlsxwriter.Workbook(xlsx_file, {'in_memory': True})
    metadata_sheet = workbook.add_worksheet('Metadata')
    targets_sheet = workbook.add_worksheet('Targets')
    # write metadata to the first workbook sheet
    _write_metadata_to_sheet(isotope_scoring_results['metadata'], workbook)
    # write lipid target data to second workbook sheet
    _write_targets_header_to_sheet(workbook)
    row = 1
    for lipid in isotope_scoring_results['targets']:
        for adduct in isotope_scoring_results['targets'][lipid]:
            for rt in isotope_scoring_results['targets'][lipid][adduct]:
                # get precursor info
                lipid_species_result = isotope_scoring_results['targets'][lipid][adduct][rt]
                if lipid_species_result is not None:
                    # make sure that the precursor was found, skip lipid completely if not
                    pre_data = [
                        lipid_species_result['precursor']['target_mz'],
                        lipid_species_result['precursor']['target_rt'],
                        lipid_species_result['precursor']['xic_peak_rt'],
                        lipid_species_result['precursor']['xic_peak_ht'],
                        lipid_species_result['precursor']['xic_peak_fwhm'],
                        lipid_species_result['precursor']['mz_cos_dist'],
                    ]
                    # iterate over db_idx, db_pos
                    for db_idx in lipid_species_result['fragments']:
                        for db_pos in lipid_species_result['fragments'][db_idx]:
                            # get the fragment data
                            fragment_data = []
                            for fragment in ['aldehyde', 'criegee']:
                                if lipid_species_result['fragments'][db_idx][db_pos][fragment]:
                                    fragment_data += [
                                        lipid_species_result['fragments'][db_idx][db_pos][fragment]['target_mz'],
                                        lipid_species_result['fragments'][db_idx][db_pos][fragment]['target_rt'],
                                        lipid_species_result['fragments'][db_idx][db_pos][fragment]['xic_peak_rt'],
                                        lipid_species_result['fragments'][db_idx][db_pos][fragment]['xic_peak_ht'],
                                        lipid_species_result['fragments'][db_idx][db_pos][fragment]['xic_peak_fwhm'],
                                        lipid_species_result['fragments'][db_idx][db_pos][fragment]['mz_cos_dist'],
                                        lipid_species_result['fragments'][db_idx][db_pos][fragment]['rt_cos_dist'],
                                    ]
                                else:
                                    fragment_data += [None, None, None, None, None, None, None]
                            fragment_found = False
                            for fd in fragment_data:
                                if fd is not None:
                                    fragment_found = True
                            if fragment_found:
                                # assemble the complete line data, write to workbook
                                line_data = [
                                    lipid, adduct, db_idx, db_pos
                                ] + pre_data + fragment_data
                                _add_line_to_targets_sheet(workbook, row, line_data)
                                row += 1
    # finish
    workbook.close()


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
    # load the data 
    oz_data = MZA(oz_data_file, cache_scan_data=True, mza_version=mza_version)
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

