"""
lipidoz/_util.py

Dylan Ross (dylan.ross@pnnl.gov)

    internal module with general utility functions
"""


import io

from matplotlib import pyplot as plt, image as mpimage

from mzapy.isotopes import monoiso_mass


def _polyunsat_ald_crg_formula(precursor_formula, db_position, db_idx):
    """
    produces the molecular formulas for aldehyde and criegee fragments corresponding to a polyunsaturated
    lipid species with double bond at the specified position (numbered from the end of the FA, i.e., n-X)
    db_idx corresponds to the index of the double bond, numbered from the end (e.g, a db_idx of 2 would
    indicate that the aldehyde and criegee fragments result from cleavage of the second double bound as
    counted from the end of the FA)

    Parameters
    ----------
    precursor_formula : ``dict(int:str)``
        molecular formula as a dictionary mapping elements (str) to their counts (int)
    db_position : ``int``
        double bond position, numbered from the end of the FA
    db_idx : ``int``
        index of double bond, numbered from the end of the FA

    Returns
    -------
    aldehyde : ``dict(int:str)``
        aldehyde fragment molecular formula as a dictionary mapping elements (str) to their counts (int)
    criegee : ``dict(int:str)``
        criegee fragment molecular formula as a dictionary mapping elements (str) to their counts (int)
    """
    aldehyde = precursor_formula.copy()
    aldehyde['C'] -= db_position
    aldehyde['H'] -= db_position * 2 - 2 * (db_idx - 1)
    aldehyde['O'] += 1
    criegee = aldehyde.copy()
    criegee['O'] += 1
    return aldehyde, criegee


def _monounsat_ald_crg_formula(precursor_formula, db_position):
    """
    produces the molecular formulas for aldehyde and criegee fragments corresponding to a monounsaturated
    lipid species with double bond at the specified position (numbered from the end of the FA, i.e., n-X)

    Parameters
    ----------
    precursor_formula : ``dict(int:str)``
        molecular formula as a dictionary mapping elements (str) to their counts (int)
    db_position : ``int``
        double bond position, numbered from the end of the FA

    Returns
    -------
    aldehyde : ``dict(int:str)``
        aldehyde fragment molecular formula as a dictionary mapping elements (str) to their counts (int)
    criegee : ``dict(int:str)``
        criegee fragment molecular formula as a dictionary mapping elements (str) to their counts (int)
    """
    # the results are the same as polyunsaturated version where db_idx == 1
    return _polyunsat_ald_crg_formula(precursor_formula, db_position, 1)


def _monounsat_ald_crg_mass(precursor_mass, db_position):
    """
    computes the monoisotopic mass corresponding to aldehyde and criegee fragment ions for an aribitrary
    lipid precursor mass

    Paramters
    ---------
    precursor_mass : ``float``
        mass of arbitrary precursor lipid
    db_position : int
        double bond position, numbered from the end of the FA

    Returns
    -------
    aldehyde : ``float``
        aldehyde fragment mass
    criegee : ``float``
        criegee fragment mass
    """ 
    # compute neutral loss formulas for ald and crg
    nl_form_ald, nl_form_crg = _monounsat_ald_crg_formula({'C': 0, 'H': 0, 'O': 0}, db_position)
    # compute netral losses and add those to the precursor
    return precursor_mass + monoiso_mass(nl_form_ald), precursor_mass + monoiso_mass(nl_form_crg) 


def _polyunsat_ald_crg_mass(precursor_mass, db_position, db_idx):
    """
    computes the monoisotopic mass corresponding to aldehyde and criegee fragment ions for an aribitrary
    lipid precursor mass

    Paramters
    ---------
    precursor_mass : ``float``
        mass of arbitrary precursor lipid
    db_position : ``int``
        double bond position, numbered from the end of the FA
    db_idx : ``int``
        index of double bond, numbered from the end of the FA

    Returns
    -------
    aldehyde : ``float``
        aldehyde fragment mass
    criegee : ``float``
        criegee fragment mass
    """ 
    # compute neutral loss formulas for ald and crg
    nl_form_ald, nl_form_crg = _polyunsat_ald_crg_formula({'C': 0, 'H': 0, 'O': 0}, db_position, db_idx)
    # compute netral losses and add those to the precursor
    return precursor_mass + monoiso_mass(nl_form_ald), precursor_mass + monoiso_mass(nl_form_crg)


def _calc_dbp_bounds(fa_carbon, fa_unsat, db_idx):
    """
    calculates the lower and upper bounds of double bond position given the number of carbons and unsaturations in 
    a fatty acid and the index of the double bond

    Parameters
    ----------
    fa_carbon : ``int``
        fatty acid carbon count
    fa_unsat : ``int``
        fatty acid unsaturation count
    db_idx : ``int``
        double bond index

    Returns
    -------
    db_pos_min : ``int``
        minimum double bond position
    db_pos_max : ``int``
        maximum double bond position
    """
    return 2 * db_idx - 1, fa_carbon - (2 * (fa_unsat - db_idx + 1))
    

def _debug_handler(debug_flag, debug_cb, msg=None, img=None):
    """
    deal with different debugging states automatically 
    debug_flag:
        None - do nothing
        'text' - produces text debugging messages only
        'full' - produces text debugging messages and displays plots (works with jupyter notebooks)
        'textcb' - produces text debugging messages but instead of printing it calls the debug_cb callback
                   with the message as an argument
    
    Parameters
    ----------
    debug_flag : ``str``
        specifies how to dispatch the message and/or plot, None to do nothing
    debug_cb : ``func``
        callback function that takes the debugging message as an argument, can be None if
        debug_flag is not set to 'textcb'
    msg : ``str``, optional
        debugging message (automatically prepended with "DEBUG: ")
    img : ``bytes``, optional
        plot image data to display
    """
    if debug_flag is not None:
        msg = 'DEBUG: ' + msg if msg is not None else None
        if debug_flag in ['text', 'full'] and msg is not None:
            print(msg)
        if debug_flag == 'full' and img is not None:
            with io.BytesIO(img) as buf:
                mimg = mpimage.imread(buf, format='png')
            plt.imshow(mimg)
            plt.show()
            plt.close()
        if debug_flag == 'textcb' and msg is not None:
            if debug_cb is not None:
                debug_cb(msg)
            else:
                ve = '_debug_handler: debug_flag was set to "textcb" but no debug_cb was provided'
                raise ValueError(ve)

