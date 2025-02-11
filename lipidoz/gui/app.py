"""
lipidoz/gui/app.py

Dylan Ross (dylan.ross@pnnl.gov)

    Main application code for LipidOz GUI

    TODO (Dylan Ross): add support for analysis of infusion data 

    TODO (Dylan Ross): once lipidoz supports looking for second order OzID fragments, need to implement probably an
                       option to look for them, as well as some logic necessary for viewing them in the results window
"""


import traceback
import queue
import threading
import pickle
import gzip
import sys

from lipidoz.workflows import run_isotope_scoring_workflow
from lipidoz.gui._config import MZA_VERSION, PROC_DEBUG
from lipidoz.gui.setup_window import SetupWindow
from lipidoz.gui.processing_window import ProcessingWindow
from lipidoz.gui.results_window import ResultsWindow


class LOzApp():
    """
    """

    def __init__(self):
        """
        """
        self.mza_version = MZA_VERSION
        # setup, processing, results windows
        self.setup_win = None
        self.processing_win = None
        self.results_win = None
        # processing parameters
        self.processing_params = None
        # results
        self.results = None

    def run(self):
        """
        main method for lipidoz app
        """
        # gather inputs using the setup window
        self.setup_win = SetupWindow()
        
        if self.setup_win.process_data:
            # if parameters for data processing provided, run the processing and display the results
            self.processing_params = self.setup_win.processing_params
            self._process_and_display_results()
        elif self.setup_win.view_results:
            results_file = self.setup_win.results_file
            # load results from file and display them
            self._load_and_display_results(results_file)
        else:
            sys.exit()

    def _process_and_display_results(self):
        """
        do the processing then display the results
        """
        # process the data using the processing window
        self._start_processing()
        # display the results using the results window
        if self.results is not None:
            self._display_results()
        else:
            sys.exit()

    def _load_and_display_results(self, results_file):
        """
        load the results from file and display them

        Parameters
        ----------
        results_file : ``str``
            lipidoz isotope scoring results file
        """
        try:
            with gzip.open(results_file, 'rb') as gzf:
                self.results = pickle.load(gzf)
        except Exception as e:
            # unable to load results
            sys.exit()
        self._display_results()

    def _processing_worker(self, results_queue, error_queue, message_queue, complete_event, early_stop_event):
        """
        """
        # TODO (Dylan Ross): handle params as a dataclass, include targeted variant of workflow function
        oz_file, target_file, rt_tol, rt_ext_win, mz_tol, d_label, d_label_in_nl = self.processing_params
        def progress_message(lipid_name, adduct, current, total):
            msg = '{} {} processing complete ({} of {})\n'
            message_queue.put(msg.format(lipid_name, adduct, current, total))
        try:
            if PROC_DEBUG:  # print complete debugging info to processing window
                results = run_isotope_scoring_workflow(oz_file, target_file, rt_tol, rt_ext_win, 
                                                    mz_tol, d_label=d_label, d_label_in_nl=d_label_in_nl,
                                                    progress_cb=progress_message, 
                                                    info_cb=lambda msg: message_queue.put(msg + '\n'),
                                                    ignore_preferred_ionization=True, 
                                                    mza_version=self.mza_version, early_stop_event=early_stop_event,
                                                    debug_flag='textcb', 
                                                    debug_cb=lambda msg: message_queue.put(msg + '\n'))
            else:
                results = run_isotope_scoring_workflow(oz_file, target_file, rt_tol, rt_ext_win, 
                                                    mz_tol, d_label=d_label, d_label_in_nl=d_label_in_nl,
                                                    progress_cb=progress_message, 
                                                    info_cb=lambda msg: message_queue.put(msg + '\n'),
                                                    ignore_preferred_ionization=True, 
                                                    mza_version=self.mza_version, early_stop_event=early_stop_event)
            # finish up
            complete_event.set()
            results_queue.put(results)
        except Exception as e:
            error_queue.put(traceback.format_exc(limit=6))
            error_queue.put(e)
            results = None

    def _start_processing(self):
        """
        starts up the isotope scoring workflow
        """
        # queue for storing the results
        results_queue = queue.Queue()
        # setup the processing window
        self.processing_win = ProcessingWindow()
        # setup and start the worker thread
        worker_thread = threading.Thread(target=self._processing_worker, 
                                         args=(results_queue, 
                                               self.processing_win.error_queue, self.processing_win.message_queue, 
                                               self.processing_win.complete_event, self.processing_win.early_stop_event))
        worker_thread.start()
        # start the processing window main loop
        self.processing_win.run()
        # get the results
        worker_thread.join()
        self.results = results_queue.get()
        # destroy the root window of the processing pane
        self.processing_win.win.destroy()
        # if anything happened except for successful processing and pressing the view results button, exit
        if not (self.processing_win.processing_complete and self.processing_win.view_results):
            sys.exit()
        
    def _display_results(self):
        """
        display results with the results window
        """
        self.results_win = ResultsWindow(self.results["isotope_scoring_results"])
        # check for back to setup flag
        if self.results_win.back_to_setup:
            # start from the beginning
            self.run()
