""" 
lipidoz_gui/processing_window.py

Dylan Ross (dylan.ross@pnnl.gov)

    processing window component of GUI for LipidOz utility

    TODO (Dylan Ross): Refactor the ProcessingWindow to be more like SetupWindow and ResultsWindow, with the individual
                       components of the UI separated into their own functions
"""

from tkinter import *
from tkinter import ttk, messagebox
import queue
import threading

from lipidoz import __version__ as LOZ_VER
from lipidoz.gui._config import ICO_PATH


class ProcessingWindow():
    """
    """

    def __init__(self):
        """
        """
        # flags indicating whether processing has completed and other statuses
        self.processing_complete = False
        self.cancelled = False
        self.view_results = False
        # attributes for threading-related stuff
        self.message_queue = queue.Queue()
        self.error_queue = queue.Queue()
        self.complete_event = threading.Event()
        self.early_stop_event = threading.Event()
        # init window
        self.win = Tk()
        self.win.title('LipidOz ({}) | Processing'.format(LOZ_VER))
        # set the window icon
        self.win.iconbitmap(ICO_PATH)
        self.win.minsize(530, 300)
        # set up window layout cols and rows
        self.win.rowconfigure(0, weight=1)
        self.win.columnconfigure(0, weight=1)
        # init the main frame that will contain all other widgets
        self.main_frm = ttk.Frame(self.win, padding=(10, 10, 10, 10))
        self.main_frm.grid(row=0, column=0, sticky=(N, S, E, W))
        # set up the main frame's grid structure
        self.main_frm.rowconfigure(0, weight=1)
        self.main_frm.rowconfigure(1)
        self.main_frm.columnconfigure(0, minsize=250, weight=1)
        # set up progress frame
        self.prog_frm = ttk.Frame(self.main_frm)
        self.prog_frm.grid(row=0, column=0, sticky=(N, S, E, W))
        self.prog_frm.rowconfigure(0, weight=1)
        self.prog_frm.columnconfigure(0, weight=1)
        self.prog_frm.columnconfigure(1)
        self.prog_frm_scroll = ttk.Scrollbar(self.prog_frm, orient='vertical')
        self.prog_frm_scroll.grid(row=0, column=1, sticky=(N, S))
        self.prog_frm_lb = Listbox(self.prog_frm, yscrollcommand=self.prog_frm_scroll.set, activestyle='none')
        self.prog_frm_lb.grid(row=0, column=0, sticky=(N, S, E, W))
        self.prog_frm_scroll.config(command=self.prog_frm_lb.yview)
        # set up the bottom frame with buttons
        self.bott_frm = ttk.Frame(self.main_frm, padding=(5, 5, 5, 5))
        self.bott_frm.columnconfigure(0, weight=1)
        self.bott_frm.columnconfigure(1, weight=1)
        self.bott_frm.grid(row=1, column=0, sticky=(W, E))
        self.bott_frm_viewres_btn = ttk.Button(self.bott_frm, 
                                               text='View Results', 
                                               command=self._viewres_btn_cb,
                                               state=DISABLED)
        self.bott_frm_viewres_btn.grid(row=0, column=1, sticky=(E,))
        self.bott_frm_cancel_btn = ttk.Button(self.bott_frm, text='Cancel', command=self._cancel_btn_cb)
        self.bott_frm_cancel_btn.grid(row=0, column=0, sticky=(W,))
        
    def _update_prog_frm(self):
        """
        """
        while self.message_queue.qsize() > 0:
            self._print_to_prog_frm(self.message_queue.get())
        if self.complete_event.is_set():
            self._processing_complete()
        if self.error_queue.qsize() > 0:
            err = self.error_queue.get()
            messagebox.showerror(title='Error while processing data', message=err)
            self.win.quit()
        else:
            # only recursively run 'after' if not complete
            self.win.after(500, self._update_prog_frm)

    def run(self):
        """ 
        start up the window main loop 
        """
        self.win.after(0, self._update_prog_frm)
        self.win.mainloop()

    def _cancel_btn_cb(self):
        """
        callback for the "Cancel" button, sets a flag and exits
        """
        self.cancelled = True
        self._close()

    def _viewres_btn_cb(self):
        """
        callback for the "View Results" button, sets a flag and exits
        """
        self.view_results = True
        self._close()

    def _print_to_prog_frm(self, line):
        """
        add line to the progress frame
        """
        self.prog_frm_lb.config(state=NORMAL)
        self.prog_frm_lb.insert(END, line)
        self.prog_frm_lb.yview(END)  # auto-scroll as lines are added
        self.prog_frm_lb.config(state=DISABLED)

    def _processing_complete(self):
        """
        sets the view results button to be enabled so user can proceed to viewing results
        """
        self.processing_complete = True
        self.bott_frm_viewres_btn.config(state=NORMAL)

    def _close(self):
        """
        handle all the possible conditions before closing the window
        """
        # if processing is not complete, set the early stop event and wait before closing
        if not self.complete_event.is_set():
            self.early_stop_event.set()
            msg = ('Early termination requested, application will exit once processing has been stopped gracefully. '
                    'This may take several seconds.')
            messagebox.showwarning(title='Stopping Processing...', message=msg)
        self.win.quit()

