``lipidoz.workflows``
=======================================
This module defines the functional components for standard high-level OzID data processing workflows. The functions
fall broadly into two categories: those related to isotope distribution analysis and those related to the machine 
learning-based double bond determination.


.. _iso-scoring-target-list:

Isotope Scoring Target List Format
-----------------------------------
The isotope scoring workflow expects a target list in *.csv* format with 3 columns: lipid name, MS adduct, and 
target retention time. A single header row from the *.csv* file is always ignored. Lines starting with *#* are 
treated as comments and ignored.

.. code-block:: none
    :caption: Example target list for isotope scoring

    lipid,adduct,retention_time
    PE(17:0_18:1),[M-H]-,23.70
    PE(17:0_20:3),[M-H]-,22.99
    PE(17:0_22:4),[M-H]-,23.46
    #CE(18:1),[M-H]-,12.34  <- this line is commented out so it will be ignored
    PG(17:0_18:1),[M-H]-,23.70
    PG(17:0_20:3),[M-H]-,22.99
    PG(17:0_22:4),[M-H]-,23.46

.. note::
    
    Target list format for :func:`lipidoz.workflows.run_isotope_scoring_workflow_infusion` is the same, 
    but excluding the retention time column, and target list format for :func:`lipidoz.workflows.run_isotope_scoring_workflow_targeted`
    is likewise the same except for the inclusion of additional columns for targeted DB indices and positions.
    See examples below.

.. code-block:: none
    :caption: Example target list for isotope scoring (infusion)

    lipid,adduct
    PE(17:0_18:1),[M-H]-
    PE(17:0_20:3),[M-H]-
    PE(17:0_22:4),[M-H]-
    #CE(18:1),[M-H]-  <- this line is commented out so it will be ignored
    PG(17:0_18:1),[M-H]-
    PG(17:0_20:3),[M-H]-
    PG(17:0_22:4),[M-H]-


.. code-block:: none
    :caption: Example target list for isotope scoring (targeted)

    lipid,adduct,retention_time,db_idx,db_pos
    PE(17:0_18:1),[M-H]-,23.70,1,9
    PE(17:0_20:3),[M-H]-,22.99,1/2/3,6/9/12
    PE(17:0_22:4),[M-H]-,23.46,1/2/3/4,3/6/9/12
    #CE(18:1),[M-H]-,12.34,1,9  <- this line is commented out so it will be ignored
    # note that multiple target DB indices/positions can be included in one line 
    # and they are separated by /
    PG(17:0_18:1),[M-H]-,23.70,1,9
    PG(17:0_20:3),[M-H]-,22.99,1/2/3,6/9/12
    PG(17:0_22:4),[M-H]-,23.46,1/2/3/4,3/6/9/12


.. _lipidoz-results-desc:

Structure of *LipidOz* Results
------------------------------------------------------------
*LipidOz* now has multiple workflows for analyzing OzID data in different ways (*e.g.* 
isotope distribution analysis, machine-learning, hybrid approach), each of which produces
its own set of results in the form of extracted/processed data and metadata. The sections
below detail the structure of those individual results sets. In order to easily organize
the different results, an overarching datastructure, termed `lipidoz_results` is defined 
which is simply a dictionary with sections for storing the results from each of the different 
individual workflows. The structure of the `lipidoz_results` is as follows:


.. code-block:: python3
    :caption: Layout of ``lipidoz_results`` dictionary

    lipidoz_results = {
        # normal/infusion/targeted variants all get packed into this one
        'isotope_scoring_results': {...isotope_scoring_results...},  
        'preml_data': {...preml_data...},
        'ml_data': np.array(...),
        # when DL inference is run, put the predictions
        # and probabilities into arrays
        # and store the name of the parameters file used
        # to run the inference
        'ml_pred_lbls': np.array(...),
        'ml_pred_probs': np.array(...),
        'ml_params_file': 'resnet18_SPLA-ULSP-BTLE_params.pt'
    }


.. _iso-scoring-workflow-results:

Structure of ``run_isotope_scoring_workflow`` Results
------------------------------------------------------------
The ``run_isotope_scoring_workflow`` function returns a dictionary containing information from double bond 
determination analyses performed for a set of lipid species defined in a target list. The results are organized 
into two top-level sections: ``'metadata'`` and ``'targets'``. The ``'metadata'`` section contains metadata 
about the analysis including information like input files and tolerances used for data extraction. The ``
'targets'`` section contains the analysis results organized in a heirarchical fashion, first by lipid, then by 
MS adduct, finally by target retention time. The results for individual lipid species (defined by a combination of 
lipid and MS adduct) are stored underneath these sub-sections.


.. note::

    See :ref:`individual-scoring-result-desc` for details regarding the organization of the result sections for 
    individual lipid species.


.. note:: 
    
    Results from :func:`lipidoz.workflows.run_isotope_scoring_workflow_targeted` are the same as for 
    :func:`lipidoz.workflows.run_isotope_scoring_workflow`, except the metadata "workflow" entry will
    be set to "isotope_scoring_targeted"


.. code-block:: python3
    :caption: Example ``run_isotope_scoring_workflow`` results dictionary

    isotope_scoring_results = {
        'metadata': {
            'workflow': 'isotope_scoring',
            'lipidoz_version': 0.4.20,
            'oz_data_file': 'data/ozid_data_file.mza',
            'target_list_file': 'a_target_list.csv',
            'rt_tol': 0.25,
            'rt_peak_win': 1.5,
            'mz_tol': 0.05,
            'd_label': None,
            'd_label_in_nl': None,
        },
        'targets': {
            'PC(16:1_16:0)': {
                '[M+H]+': {
                    '21.05min': {
                        'precursor': {
                            'target_mz': 789.0123,
                            'target_rt': 23.45,
                            'xic_peak_rt': 23.45,
                            'xic_peak_ht': 1e5,
                            'xic_peak_fwhm': 0.15,
                            'mz_ppm': 10.1,
                            'abun_percent': 5.5,
                            'mz_cos_dist': 0.15,
                            'isotope_dist_img': ...,
                            'xic_fit_img': ...,
                            'saturation_corrected': False
                        },
                        'fragments': {
                            1: {
                                9: {
                                    'aldehyde': {
                                        'target_mz': 234.5678,
                                        'target_rt': 23.45,
                                        'xic_peak_rt': 23.45,
                                        'xic_peak_ht': 1e4,
                                        'xic_peak_fwhm': 0.25,
                                        'mz_ppm': 10.1,
                                        'abun_percent': 5.5,
                                        'mz_cos_dist': 0.15,
                                        'rt_cos_dist': 0.25,
                                        'isotope_dist_img': ...,
                                        'xic_fit_img': ...,
                                        'saturation_corrected': False,
                                    },
                                    # if the fragment was not found the section is set to None
                                    'criegee': None  
                                },
                                # more db positions ...
                            },
                            # more db indices ...
                        }
                    },
                    # more retention times ...
                },
                # more adducts ...
            },
            # more targets ...
        },
    }


.. _iso-scoring-workflow-results-inf:

Structure of ``run_isotope_scoring_workflow_infusion`` Results
----------------------------------------------------------------
The results from the infusion variant of the isotope scoring workflow are very similar
to those from the normal version, except any component having to do with retention time
is omitted.


.. code-block:: python3
    :caption: Example ``run_isotope_scoring_workflow_infusion`` results dictionary

    isotope_scoring_results = {
        'metadata': {
            'workflow': 'isotope_scoring_infusion',
            'lipidoz_version': 0.4.20,
            'oz_data_file': 'data/infusion_ozid_data_file.mza',
            'target_list_file': 'a_target_list.csv',
            'mz_tol': 0.05,
            'd_label': None,
            'd_label_in_nl': None,
        },
        'targets': {
            'PC(16:1_16:0)': {
                '[M+H]+': {
                    'infusion': {  # instead of a retention time the label here is just "infusion"
                        'precursor': {
                            'target_mz': 789.0123,
                            'mz_ppm': 10.1,
                            'abun_percent': 5.5,
                            'mz_cos_dist': 0.15,
                            'isotope_dist_img': ...,
                        },
                        'fragments': {
                            1: {
                                9: {
                                    'aldehyde': {
                                        'target_mz': 234.5678,
                                        'mz_ppm': 10.1,
                                        'abun_percent': 5.5,
                                        'mz_cos_dist': 0.15,
                                        'isotope_dist_img': ...,
                                    },
                                    # if the fragment was not found the section is set to None
                                    'criegee': None  
                                },
                                # more db positions ...
                            },
                            # more db indices ...
                        }
                    }
                },
                # more adducts ...
            },
            # more targets ...
        },
    }


.. _preml-dataset-desc:

Structure of ``collect_preml_dataset`` dataset
---------------------------------------------------
The :func:`lipidoz.workflows.collect_preml_dataset` function returns a dictionary containing minimally 
processed RTMZ data for a set of 
lipid species defined in a target list. The dataset contains extracted data for lipid precursor and aldehyde/criegee
OzID fragments for different double bond locations. The dataset is organized 
into two top-level sections: ``'metadata'`` and ``'targets'``. The ``'metadata'`` section contains metadata 
about the analysis including information like input files and tolerances used for data extraction. 
The ``'targets'`` section contains the data for individual lipid species, defined by the lipid, MS adduct, target
retention time, double bond index, and double bond position.

.. code-block:: python3
    :caption: Example ``collect_preml_data`` dataset

    pre_ml_dataset = {
        'metadata': {
            'workflow': 'pre_ml',
            'lipidoz_version': '0.4.20',
            'oz_data_file': '../../_data/Ultimate-Splash_NEG_O3_Run-1.mza',
            'target_list_file': 'test_target_list.csv', 
            'rt_tol': 0.2, 
            'd_label': 5, 
            'd_label_in_nl': False,
        },
        'targets': {
            'PE(18:1_17:0)|[M+H]+|23.70min|1|1': {  # <lipid>|<adduct>|<target_rt>|<db_idx>|<db_pos>
                'pre_data': # <raw RTMZ arrays for precursor>
                'ald_data': # <raw RTMZ arrays for aldehyde OzID fragment>
                'crg_data': # <raw RTMZ arrays for criegee OzID fragment>
                'pre_mz': mz,  # precursor m/z
                'ald_mz': ald_mz,  # aldehyde OzID fragment m/z
                'crg_mz': crg_mz,  # criegee OzID fragment m/z
                'rt': 23.70,  # target retention time
            },
            # ... data for other targets omitted
        },
    } 


Module Reference
---------------------------------------

Isotope Distribution Analysis
+++++++++++++++++++++++++++++++++++++++
.. autofunction:: lipidoz.workflows.run_isotope_scoring_workflow

.. autofunction:: lipidoz.workflows.run_isotope_scoring_workflow_targeted

.. autofunction:: lipidoz.workflows.run_isotope_scoring_workflow_infusion

.. autofunction:: lipidoz.workflows.save_isotope_scoring_results

.. autofunction:: lipidoz.workflows.write_isotope_scoring_report_xlsx

Machine Learning
+++++++++++++++++++++++++++++++++++++++
.. autofunction:: lipidoz.workflows.collect_preml_dataset

.. note::

    See :ref:`preml-dataset-desc` for details regarding the organization of the pre-ml dataset.

.. autofunction:: lipidoz.workflows.convert_multi_preml_datasets_labeled

.. autofunction:: lipidoz.workflows.convert_multi_preml_datasets_unlabeled


Hybrid Workflow
+++++++++++++++++++++++++++++++++++++++

.. note::
    
    The hybrid workflow uses ML inference to prioritize targets for full
    analysis using the targeted variant of the isotope distribution analysis.
    See :ref:`lipidoz-results-desc` for details on how the data from these 
    different steps is organized in the `lipidoz_results` dictionary that is
    returned by this function. 

.. autofunction:: lipidoz.workflows.hybrid_deep_learning_and_isotope_scoring

