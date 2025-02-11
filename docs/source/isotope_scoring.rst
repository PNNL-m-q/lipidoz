``lipidoz.isotope_scoring``
=======================================
This module contains a function for performing lipid double bond position determination from OzID data using 
a method that is based on comparing observed and theoretical isotope distributions for putatative Oz fragments.


.. _individual-scoring-result-desc:

Structure of ``score_db_pos_isotope_dist_polyunsat`` Results
-------------------------------------------------------------
The :func:`lipidoz.isotope_scoring.score_db_pos_isotope_dist_polyunsat` function returns a dictionary containing all of the output 
information from a double bond position determination analysis for a single lipid species. The dictionary 
is organized into a few top-level sections: one for information related to the precursor (dictionary key 
``'precursor'``), and the others with information related to putative OzID fragments for each acyl chain 
(dictionary key ``'fragments'``). The fragments sections is further divided in a heirarchical 
fashion, first by double bond index, second by double bond position, with information related to the 
corresponding putative OzID fragments stored in sections ``'aldehyde'`` and ``'criegee'`` under that. For 
example, ``results['fragments'][2][9]['aldehyde']`` would contain information about the putative aldehyde OzID 
fragment corresponding to a double bond at position 9 with index of 2 (*i.e.*, the second double bond ordered 
from the end of the acyl chain). 

.. code-block:: python3
    :caption: Example ``score_db_pos_isotope_dist_polyunsat`` results dictionary

    results = {
        'precursor': {  # data associated with the lipid precursor
            'target_mz': 853.4567,
            'target_rt': 23.70,
            'xic_peak_rt': 23.71,  # XIC fitted peak parameter
            'xic_peak_ht': 1.23e5,  # XIC fitted peak parameter
            'xic_peak_fwhm': 1.4,  # XIC fitted peak parameter
            'mz_ppm': 11.621,  # isotope distribution scoring component
            'abun_percent': 4.15,  # isotope distribution scoring component
            'mz_cos_dist': 0.030,  # isotope distribution scoring component
            'isotope_dist_img': b'...',  # image in .png format stored as bytes
            'xic_fit_img': b'...',  # image in .png format stored as bytes
            'saturation_corrected': True
        },
        'fragments': {
            1: {
                # ... results for other double bond indices omitted
            },
            2: {
                # ... results for other double bond positions omitted
                8: {
                    'aldehyde': None,  # aldehyde fragment not found for this db position
                    'criegee': None,  # criegee fragment not found for this db position
                },
                9: {
                    'aldehyde': {
                        'target_mz': 625.424762,
                        'target_rt': 23.7,
                        'xic_peak_rt': 23.717556783082685,  # XIC fitted peak parameter
                        'xic_peak_ht': 12661.977204378885,  # XIC fitted peak parameter
                        'xic_peak_fwhm': 0.17212908575224203,  # XIC fitted peak parameter
                        'mz_ppm': 5.196747765566774,  # isotope distribution scoring component
                        'abun_percent': 1.3707328733140074,  # isotope distribution scoring component
                        'mz_cos_dist': 0.025067410332574203,  # isotope distribution scoring component
                        'rt_cos_dist': 0.13721195277045561,  # retention time agreement with precursor
                        'isotope_dist_img': b'...',  # image in .png format stored as bytes
                        'xic_fit_img': b'...',  # image in .png format stored as bytes
                        'saturation_corrected': False,
                    },
                    'criegee': {
                        'target_mz': 641.419677,
                        'target_rt': 23.7,
                        'xic_peak_rt': 23.72089709916982,  # XIC fitted peak parameter
                        'xic_peak_ht': 25546.76018712644,  # XIC fitted peak parameter
                        'xic_peak_fwhm': 0.17084010175727868,  # XIC fitted peak parameter
                        'mz_ppm': 7.422567626157515,  # isotope distribution scoring component
                        'abun_percent': 0.7105029034931888,  # isotope distribution scoring component
                        'mz_cos_dist': 0.03100184607806855,  # isotope distribution scoring component
                        'rt_cos_dist': 0.08396320539908453,  # retention time agreement with precursor
                        'isotope_dist_img': b'...',  # image in .png format stored as bytes
                        'xic_fit_img': b'...',  # image in .png format stored as bytes
                        'saturation_corrected': False,
                    },
                },
                10: {
                    'aldehyde': None,  # aldehyde fragment not found for this db position
                    'criegee': {
                        'target_mz': 627.404027,
                        'target_rt': 23.7,
                        'xic_peak_rt': 23.200000000000003,  # XIC fitted peak parameter
                        'xic_peak_ht': 1320.5005362849708,  # XIC fitted peak parameter
                        'xic_peak_fwhm': 0.21637538953545316,  # XIC fitted peak parameter
                        'mz_ppm': 31.88232636309292,  # isotope distribution scoring component
                        'abun_percent': 2.9104719246527786,  # isotope distribution scoring component
                        'mz_cos_dist': 0.20026290761290033,  # isotope distribution scoring component
                        'rt_cos_dist': 0.9794397320037685,  # retention time agreement with precursor
                        'isotope_dist_img': b'...',  # image in .png format stored as bytes
                        'xic_fit_img': b'...',  # image in .png format stored as bytes
                        'saturation_corrected': False,
                    },
                },
                # ... results for other double bond positions omitted 
            },
            3: {
                # ... results for other double bond indices omitted
            },
        },
    }


.. note::
    
    Results from :func:`lipidoz.isotope_scoring.score_db_pos_isotope_dist_polyunsat_infusion` are the same as for 
    :func:`lipidoz.isotope_scoring.score_db_pos_isotope_dist_polyunsat`, except all information related to retention 
    time, XIC, *etc.* are omitted. Likewise, results from :func:`lipidoz.isotope_scoring.score_db_pos_isotope_dist_polyunsat`
    and :func:`lipidoz.isotope_scoring.score_db_pos_isotope_dist_targeted` are in exactly the same format, the targeted
    variant just has results for fewer double bond positions.


Module Reference
---------------------------------------
.. autofunction:: lipidoz.isotope_scoring.score_db_pos_isotope_dist_polyunsat

.. autofunction:: lipidoz.isotope_scoring.score_db_pos_isotope_dist_targeted

.. autofunction:: lipidoz.isotope_scoring.score_db_pos_isotope_dist_polyunsat_infusion
