==============================
Isotope Distribution Analysis
==============================
Isotope distribution analysis enables determination of lipid double bond position through examination of 
the isotope distributions for putative OzID fragments from LC-OzID-MS data. This analysis is implemented
in :func:`lipidoz.workflows.run_isotope_scoring_workflow`.


Inputs
------------------------------
* **OzID data in .mza format**
* **Lipid target list in .csv format † with the following columns:**
    * *lipid* -- lipid name in standard abbreviated format
    * *adduct* -- MS ionization state
    * *retention time* -- observed retention time
* **Parameters controlling data extraction/interpretation:**
    * *RT tolerance* -- determines the size of the retention time window (around fitted retention time) that is used for extraction of precursor and OzID fragment MS1 spectra
    * *RT extraction window* -- determines the size of the window (around target retention time) used for extraction of LC chromatograms for retention time fitting
    * *m/z tolerance* -- m/z tolerance for extraction of LC chromatograms
    * *deuterium label* -- count of deuteriums in labeled lipids (*e.g.*, in SPLASH LipidoMIX standards)
    * *deuterium label in NL* -- indicates whether the deuteriums in labeled lipids are present in the neutral loss from OzID fragmentation (and therefore must be accounted for in the fragment formulae)

† *The target list is expected to have one header row, and rows beginning with "#" are treated as comments and skipped*


Output
------------------------------
* **LipidOz Analysis Results** -- Full analysis results including plots and scoring information for all lipid targets, organized in a hierarchical dictionary datastructure (see :ref:`iso-scoring-workflow-results` for details). These results can be saved in a binary format for viewing later (use function: :func:`lipidoz.workflows.save_isotope_scoring_results`) or exported in tabular format as an Excel spreadsheet (use function: :func:`lipidoz.workflows.write_isotope_scoring_report_xlsx`).


Examples
------------------------------

.. code-block:: python3
    :caption: example of isotope distribution analysis with labeled standards

    from lipidoz.workflows import (
        run_isotope_scoring_workflow, save_isotope_scoring_results, write_isotope_scoring_report
    )


    def main():

        # OzID data file, in .mza format
        oz_data_file = 'data/SPLASH_d7_pos_Oz.mza'

        # target list in .csv format
        target_file = 'splash_d7_pos.csv'

        # setup parameters controlling data extraction/interpretation
        rt_tol = 0.2
        rt_ext_win = 1.5 
        mz_tol = 0.01
        # these parameters are specific for describing deuterium labels
        d_label = 7
        d_label_in_nl = True

        # run the analysis
        results = run_isotope_scoring_workflow(oz_data_file, target_file, 
                                               rt_tol, rt_ext_win, mz_tol,
                                               d_label=d_label, d_label_in_nl=d_label_in_nl)

        # when analysis is complete, save and export the results
        save_isotope_scoring_results(results, 'splash_d7_pos.loz')
        write_isotope_scoring_report_xlsx(results, 'splash_d7_pos_iso_scoring.xlsx')


    if __name__ == '__main__':
        main()


.. code-block:: python3
    :caption: example of isotope distribution analysis using the progress callback

    from lipidoz.workflows import (
        run_isotope_scoring_workflow, save_isotope_scoring_results, write_isotope_scoring_report
    )

    def progcb(lipid_name, adduct, pos, tot):
        # simple callback function that just prints out info
        s = 'lipid: {} adduct: {} ({} of {})'
        print(s.format(lipid_name, adduct, pos, tot), flush=True)


    def main():

        # OzID data file, in .mza format
        oz_data_file = 'data/total_lipids_1234_pos_Oz.mza'

        # target list in .csv format
        target_file = 'total_lipids_1234_pos.csv'

        # setup parameters controlling data extraction/interpretation
        rt_tol = 0.2
        rt_ext_win = 1.5 
        mz_tol = 0.01

        # run the analysis
        # here we will use the progress callback to report on progress after each lipid target
        results = run_isotope_scoring_workflow(oz_data_file, target_file, 
                                               rt_tol, rt_ext_win, mz_tol,
                                               progress_cb=progcb)

        # when analysis is complete, save and export the results
        save_isotope_scoring_results(results, 'total_lipids_1234_pos.loz')
        write_isotope_scoring_report_xlsx(results, 'total_lipids_1234_pos_iso_scoring.xlsx')


    if __name__ == '__main__':
        main()



