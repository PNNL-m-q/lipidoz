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
    but excluding the retention time column

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
    
    Results from :func:`lipidoz.workflows.run_isotope_scoring_workflow_infusion` are the same as for 
    :func:`lipidoz.workflows.run_isotope_scoring_workflow`, except all information related to retention 
    time are omitted. 


.. code-block:: python3
    :caption: Example ``run_isotope_scoring_workflow`` results dictionary

    results = {
        'metadata': {
            'workflow': 'isotope_scoring',
            'lipidoz_version': '0.4.20',
            'oz_data_file': '../../_data/Ultimate-Splash_NEG_O3_Run-1.mza',
            'target_list_file': 'test_target_list.csv',
            'rt_tol': 0.2,  # retention time tolerance used for extracting MS1 spectra
            'rt_peak_win': 1.5,  # size of retention time window to extract for XIC fitting
            'mz_tol': 0.01,  # m/z tolerance for XIC extraction
            'd_label': 5,  # number of deuteriums for labeled lipid standards
            'd_label_in_nl': False,  # the deuterium labels are not part of the neutral loss
        },
        'targets': {
            'PE(18:1_17:0)': {
                '[M-H]-': {
                    '20.00min': {
                            # ... individual lipid species results
                        }
                    },
                    '22.22min': {
                            # ... individual lipid species results
                        }
                    },
                },
            },
            'PE(20:3_17:0)': {
                '[M-H]-': {
                    '21.21min': {
                            # ... individual lipid species results
                        }
                    },
                },
                '[M-HCOO]-': {
                    # ... individual lipid species results
                },
            },
            'PE(22:4_17:0)': {
                '[M-H]-': {
                    '17.38min': {
                            # ... individual lipid species results
                        }
                    },
                },
            },
            # ... results for other target lipids omitted
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
