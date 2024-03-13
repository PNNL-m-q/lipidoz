""" 
lipidoz_gui/setup_window.py

Dylan Ross (dylan.ross@pnnl.gov)

    Setup window component of GUI for LipidOz utility

"""


from re import sub
from tkinter import *
from tkinter import ttk, filedialog, messagebox
from pathlib import Path
import os
import signal
import subprocess

from lipidoz import __version__ as LOZ_VER
from lipidoz.gui._config import DEFAULT_RT_TOL, DEFAULT_RT_EXT_WIN, DEFAULT_MZ_TOL, ICO_PATH, PLATFORM_IS_WINDOWS

        
class SetupWindow():
    """
    """

    def __init__(self):
        """
        """
        # status flags
        self.view_results = False
        self.process_data = False
        # input values (data processing)
        self.processing_params = None
        # input values (load results)
        self.results_file = None
        # init window and configure
        self.win = Tk()
        self._configure_window()
        # init and configure the main frame that will contain all other widgets
        self._setup_main_frame()
        # start up the window main loop
        self.win.mainloop()

    def _configure_window(self):
        # configure window
        self.win.title('LipidOz ({}) | Setup'.format(LOZ_VER))
        # set the window icon
        self.win.iconbitmap(ICO_PATH)
        # set up window layout cols and rows
        self.win.rowconfigure(0, weight=1)
        self.win.columnconfigure(0, weight=1)

    def _setup_main_frame(self):
        """
        main frame has 3 rows and 1 column
            row 0 -> upper frame
            row 1 -> horizontal separator
            row 2 -> lower frame
        """
        self.main_frm = ttk.Frame(self.win, padding=(10, 10, 10, 10))
        self.main_frm.grid(row=0, column=0, sticky=(N, S, E, W))
        # set up the main frame's grid structure
        self.main_frm.rowconfigure(0)
        self.main_frm.rowconfigure(1)
        self.main_frm.rowconfigure(2)
        self.main_frm.columnconfigure(0, minsize=350, weight=1)
        # add horizontal separator
        self.main_frm_hsep = ttk.Separator(self.main_frm, orient=HORIZONTAL)
        self.main_frm_hsep.grid(row=1, column=0, sticky=(E, W))
        # set up the each of the sub-frames
        self._setup_upper_frame()
        self._setup_lower_frame()

    def _setup_upper_frame(self):
        """
        upper frame has 9 rows and 4 columns
            row 0, column 0 -> "Process OzID Data"
            row 1, column 0 -> "OzID data file"
            row 1, column 1 -> input box
            row 1, column 2 -> [...] button -> file select dialog (.mza)
            row 1, column 3 -> [...] button -> file select dialog (convert .d)
            row 2, column 0 -> "target list"  TODO (Dylan Ross): add checkbox for targeted variant
            row 2, column 1 -> input box
            row 2, column 2 -> [...] button -> file select dialog
            row 3, column 0 -> "RT tolerance"
            row 3, column 1 -> input box
            row 3, column 2 -> "min"
            row 4, column 0 -> "RT extraction window"
            row 4, column 1 -> input box
            row 4, column 2 -> "min"
            row 5, column 0 -> "m/z tolerance"
            row 5, column 1 -> input box
            row 5, column 2 -> "Da"
            row 6, column 0 -> "deuterium label"
            row 6, column 1 -> input box
            row 7, column 0 -> "deuterium label in neutral loss"
            row 7, column 1 -> checkbox
            row 8, column 2 -> [Process Data] button
        """
        self.upper_frm = ttk.Frame(self.main_frm, padding=(10, 10, 10, 10))
        self.upper_frm.grid(row=0, column=0, sticky=(N, S, E, W))
        # configure grid
        self.upper_frm.rowconfigure(0)
        self.upper_frm.rowconfigure(1)
        self.upper_frm.rowconfigure(2)
        self.upper_frm.rowconfigure(3)
        self.upper_frm.rowconfigure(4)
        self.upper_frm.rowconfigure(5)
        self.upper_frm.rowconfigure(6)
        self.upper_frm.rowconfigure(7)
        self.upper_frm.rowconfigure(8)
        self.upper_frm.columnconfigure(0)
        self.upper_frm.columnconfigure(1)
        self.upper_frm.columnconfigure(2)
        self.upper_frm.columnconfigure(3)
        # add all of the labels
        self._add_labels_to_upper_frm()
        # add the file select [...] buttons
        self.upper_frm_ozfilesel_mza_btn = ttk.Button(self.upper_frm, text='.mza', command=self._ozfilesel_mza_cb)
        self.upper_frm_ozfilesel_mza_btn.grid(row=1, column=2, sticky=(W, E))
        self.upper_frm_ozfilesel_d_btn = ttk.Button(self.upper_frm, text='.d', command=self._ozfilesel_d_cb)
        self.upper_frm_ozfilesel_d_btn.grid(row=1, column=3, sticky=(W,))
        # disable .d file selection if not running on windows
        if not PLATFORM_IS_WINDOWS:
            self.upper_frm_ozfilesel_d_btn["state"] = "disabled"
        self.upper_frm_tarlstsel_btn = ttk.Button(self.upper_frm, text='...', command=self._tarlstsel_cb)
        self.upper_frm_tarlstsel_btn.grid(row=2, column=2, sticky=(W,))
        # add the entries
        self._add_entries_to_upper_frame()
        # add d label in nl checkbutton
        self.upper_frm_dlblnl_cbn_var = IntVar()
        self.upper_frm_dlblnl_cbn = ttk.Checkbutton(self.upper_frm, variable=self.upper_frm_dlblnl_cbn_var)
        self.upper_frm_dlblnl_cbn.grid(row=7, column=1, sticky=(W,))
        # add process data button
        self.upper_frm_procdata_btn = ttk.Button(self.upper_frm, text='Process Data', command=self._procdata_btn_cb)
        self.upper_frm_procdata_btn.grid(row=8, column=2, sticky=(E,))

    def _ozfilesel_mza_cb(self):
        """
        callback for when ozid data file is selected (.mza)
        """
        # offer a dialog that accepts .mza
        oz_file = filedialog.askopenfilename(title='Select OzID data file', 
                                             filetypes=[('MZA', '.mza')])
        self.upper_frm_ent_vars[0].set(oz_file)

    def _ozfilesel_d_cb(self):
        """
        callback for when ozid data file is selected (agilent .d)
        """
        # determine whether agilent .d conversion is supported (mza_converter is available)
        # offer a filedialog that accepts mza and agilent .d if so, otherwise just offer a dialog that accepts .mza
        oz_file_d = filedialog.askdirectory(title='Select OzID data file for conversion', 
                                          mustexist=True)
        # try file conversion
        oz_file = self._convert_agilent_d(oz_file_d)
        # set the file selection to the newly converted .mza file
        self.upper_frm_ent_vars[0].set(oz_file)

    def _convert_agilent_d(self, oz_file_d):
        """ 
        try to convert agilent .d into .mza with mza.exe
        """
        # do not try to convert if the selected directory does not end in .d
        if os.path.splitext(oz_file_d)[-1] != '.d':
            msg = 'Oz file name "{}" not valid. Must end with ".d"'
            messagebox.showerror(title='Cannot convert Agilent .d', message=msg.format(oz_file_d))
            return ''
        # make sure mza.exe is in the same dir as the data
        mza_converter = os.path.join(Path(oz_file_d).parent.absolute(), 'mza.exe')
        if not os.path.isfile(mza_converter):
            msg = 'cannot find MZA converter (mza.exe), must be in same directory as data'
            messagebox.showerror(title='Cannot convert Agilent .d', message=msg)
            return ''
        # before converting, check that the output file does not already exist
        output = os.path.splitext(oz_file_d)[0] + '.mza'
        if os.path.isfile(output):
            msg = 'Oz file "{}" seems to already have been converted to MZA, using "{}" as input'
            messagebox.showwarning(title='File already converted', message=msg.format(oz_file_d, output))
            return output
        # put up a window with output from mza.exe output
        conv_win = _MzaConversionWindow(self.win, oz_file_d, mza_converter)
        # get the name of the converted oz file (will be '' if conversion didnt work)
        return conv_win.oz_file

    def _tarlstsel_cb(self):
        """
        callback for when target list file is selected
        """
        # TODO (Dylan Ross): check the state of targeted checkbox and store that so that the targeted
        #                    variant of the workflow function gets used in the processing window
        target_list = filedialog.askopenfilename(title='Select target list file', 
                                             filetypes=[('Comma separated values', '.csv')])
        self.upper_frm_ent_vars[1].set(target_list)

    def _procdata_btn_cb(self):
        """
        callback for the "Process Data" button, sets a flag and exits
        """
        # check that the data file and target list values are set before proceeding

        self.process_data = True
        # store the entry values
        self.processing_params = (self.upper_frm_ent_vars[0].get(),
                                  self.upper_frm_ent_vars[1].get(),
                                  self.upper_frm_ent_vars[2].get(),
                                  self.upper_frm_ent_vars[3].get(),
                                  self.upper_frm_ent_vars[4].get(),
                                  self.upper_frm_ent_vars[5].get(),
                                  bool(self.upper_frm_dlblnl_cbn_var.get()))
        self._close()

    def _add_labels_to_upper_frm(self):
        """
        """
        self.upper_frm_lbls = [None, None, None, None, None, None, None, None]
        labels = ['Process OzID Data', 
                  'OzID data file', 'target list', 'RT tolerance', 'RT extraction window', 
                  'm/z tolerance', 'deuterium label', 'deuterium label in NL']
        for i, label in enumerate(labels):
            if i == 0:
                self.upper_frm_lbls[i] = ttk.Label(self.upper_frm, text=label, 
                                                   relief=GROOVE, borderwidth=2, padding=(5, 5, 5, 5))
            else:
                self.upper_frm_lbls[i] = ttk.Label(self.upper_frm, text=label)
            self.upper_frm_lbls[i].grid(row=i, column=0, sticky=(E,))

    def _add_entries_to_upper_frame(self):
        """
        """
        self.upper_frm_ents = [None, None, None, None, None, None]
        ent_widths = [40, 40, 5, 5, 5, 5]
        stickies = [(W, E), (W, E), (W,), (W,), (W,), (W,)]
        self.upper_frm_ent_vars = [StringVar(), StringVar(), DoubleVar(), DoubleVar(), DoubleVar(), IntVar()]
        for i in range(6):
            self.upper_frm_ents[i] = ttk.Entry(self.upper_frm, justify=RIGHT, 
                                               width=ent_widths[i], textvariable=self.upper_frm_ent_vars[i])
            self.upper_frm_ents[i].grid(row=i + 1, column=1, sticky=stickies[i])
        # set a few default values
        self.upper_frm_ent_vars[2].set(DEFAULT_RT_TOL)
        self.upper_frm_ent_vars[3].set(DEFAULT_RT_EXT_WIN)
        self.upper_frm_ent_vars[4].set(DEFAULT_MZ_TOL)

    def _setup_lower_frame(self):
        """
        lower frame has 2 rows and 3 columns
            row 0, column 0 -> "Load Previous Results"
            row 1, column 0 -> "isotope scoring results"
            row 1, column 1 -> input box
            row 1, column 2 -> [...] button -> file select dialog
            row 2, column 2 -> [View Results] button
        """
        self.lower_frm = ttk.Frame(self.main_frm, padding=(10, 10, 10, 10))
        self.lower_frm.grid(row=2, column=0, sticky=(N, S, E, W))
        # configure grid
        self.lower_frm.rowconfigure(0)
        self.lower_frm.rowconfigure(1)
        self.lower_frm.rowconfigure(2)
        self.lower_frm.columnconfigure(0)
        self.lower_frm.columnconfigure(1)
        self.lower_frm.columnconfigure(2)
        # add label
        self.lower_frm_lbl_1 = ttk.Label(self.lower_frm, text='View Previous Results', 
                                         relief=GROOVE, borderwidth=2, padding=(5, 5, 5, 5))
        self.lower_frm_lbl_1.grid(row=0, column=0)
        self.lower_frm_lbl_2 = ttk.Label(self.lower_frm, text='results file')
        self.lower_frm_lbl_2.grid(row=1, column=0, sticky=(E,))
        # add input box
        self.lower_frm_ent_var = StringVar()
        self.lower_frm_ent = ttk.Entry(self.lower_frm, justify=LEFT, width=40, textvariable=self.lower_frm_ent_var)
        self.lower_frm_ent.grid(row=1, column=1, sticky=(W, E))
        # add file select button
        self.lower_frm_resfilesel_btn = ttk.Button(self.lower_frm, text='...', command=self._resfilesel_cb)
        self.lower_frm_resfilesel_btn.grid(row=1, column=2, sticky=(W,))
        # add view results button
        self.lower_frm_viewews_btn = ttk.Button(self.lower_frm, 
                                                text='View Results', 
                                                command=self._viewres_btn_cb)
        self.lower_frm_viewews_btn.grid(row=2, column=2, sticky=(E,))

    def _viewres_btn_cb(self):
        """
        callback for the "View Results" button, sets a flag and exits
        """
        self.view_results = True
        self.results_file = self.lower_frm_ent_var.get()
        self._close()

    def _resfilesel_cb(self):
        """
        callback for when results file is selected
        """
        res_file = filedialog.askopenfilename(title='Select results file', 
                                              filetypes=[('isotope scoring results', '.loz')])
        self.lower_frm_ent_var.set(res_file)

    def _close(self):
        """
        handle all the possible conditions before closing the window
        """
        self.win.quit()
        self.win.destroy()
  

class _MzaConversionWindow():
    """ simple window that displays output from mza conversion """

    def __init__(self, win, oz_file_d, mza_converter):
        """  """
        # this will be set if the file conversion is successful
        self.oz_file = ''
        # store the input file and MZA converter
        self.oz_file_d = oz_file_d
        self.mza_converter = mza_converter
        # store the conversion command
        self.cmd = [self.mza_converter, '-file', self.oz_file_d, '-intensityThreshold', '20']
        self.result = None
        # init window
        self.win = Toplevel(win)
        self.win.title('LipidOz | MZA Conversion')
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
        self.main_frm.rowconfigure(2)
        self.main_frm.columnconfigure(0, minsize=250, weight=1)
        # set up progress frame
        self.prog_frm = ttk.Frame(self.main_frm)
        self.prog_frm.grid(row=0, column=0, sticky=(N, S, E, W))
        self.prog_frm.rowconfigure(0, weight=1)
        self.prog_frm.columnconfigure(0, weight=1)
        self.prog_frm.columnconfigure(1)
        self.prog_frm_scroll = ttk.Scrollbar(self.prog_frm, orient='vertical')
        self.prog_frm_scroll.grid(row=0, column=1, sticky=(N, S))
        self.prog_frm_scroll2 = ttk.Scrollbar(self.prog_frm, orient='horizontal')
        self.prog_frm_scroll2.grid(row=1, column=0, sticky=(E, W))
        self.prog_frm_lb = Listbox(self.prog_frm, 
                                   xscrollcommand=self.prog_frm_scroll2.set, 
                                   yscrollcommand=self.prog_frm_scroll.set, 
                                   activestyle='none')
        self.prog_frm_lb.grid(row=0, column=0, sticky=(N, S, E, W))
        self.prog_frm_scroll.config(command=self.prog_frm_lb.yview)
        self.prog_frm_scroll2.config(command=self.prog_frm_lb.xview)
        # set up the bottom frame with buttons
        self.bott_frm = ttk.Frame(self.main_frm, padding=(5, 5, 5, 5))
        self.bott_frm.columnconfigure(0, weight=1)
        self.bott_frm.columnconfigure(1, weight=1)
        self.bott_frm.grid(row=2, column=0, sticky=(W, E))
        self.bott_frm_ok_btn = ttk.Button(self.bott_frm, 
                                          text='Ok', 
                                          command=self._ok_btn_cb,
                                          state=DISABLED)
        self.bott_frm_ok_btn.grid(row=0, column=1, sticky=(E,))
        self.bott_frm_cancel_btn = ttk.Button(self.bott_frm, text='Cancel', command=self._cancel_btn_cb)
        self.bott_frm_cancel_btn.grid(row=0, column=0, sticky=(W,))
        self._run()

    def _run(self):
        """ 
        start up the window main loop 
        """
        # print some initial text to the window before trying to run the command
        self._print_to_prog_frm('Converting {} to MZA with mza.exe'.format(self.oz_file_d))
        self._print_to_prog_frm('------------------------------')
        self._print_to_prog_frm('COMMAND:')
        self._print_to_prog_frm(' '.join(self.cmd))
        self._print_to_prog_frm('------------------------------')
        self._do_conversion()
        self.win.mainloop()

    def _do_conversion(self):
        """ """
        # try running the conversion
        self.result = subprocess.Popen(self.cmd, 
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       universal_newlines=True)
        # check progress every second
        self.win.after(1000, self._check_conversion_complete)

    def _check_conversion_complete(self):
        if self.result.poll() is not None:
            # complete
            self._conversion_complete()
        else:
            # check again in a second
            self.win.after(500, self._check_conversion_complete)

    def _conversion_complete(self):
        """ """
        stdout, stderr = self.result.communicate()
        self._print_to_prog_frm('RETURN CODE:')
        self._print_to_prog_frm(str(self.result.returncode))
        self._print_to_prog_frm('------------------------------')
        self._print_to_prog_frm('STDOUT:')
        self._print_to_prog_frm(stdout)
        self._print_to_prog_frm('------------------------------')
        self._print_to_prog_frm('STDERR:')
        self._print_to_prog_frm(stderr)
        self._print_to_prog_frm('------------------------------')
        self.bott_frm_ok_btn.config(state=NORMAL)
        output = os.path.splitext(self.oz_file_d)[0] + '.mza'
        if self.result.returncode == 0 and os.path.isfile(output):
            self._print_to_prog_frm('>>> SUCCESS <<<')
            self.oz_file = output
        else:
            self._print_to_prog_frm('>>> FAILED <<<')
            messagebox.showerror(message='MZA conversion failed')

    def _ok_btn_cb(self):
        """ callback for ok button """
        self._close()

    def _cancel_btn_cb(self):
        """ callback for cancel button """
        messagebox.showwarning(message='Cancelled MZA conversion')
        if self.result is not None:
            # stop processing if its happening
            self.result.send_signal(signal.SIGTERM)
        self._close()

    def _print_to_prog_frm(self, line):
        """
        add line to the progress frame
        """
        # auto wrap text to 50 characters
        if len(line) <= 50:
            lines = [line]
        else:
            lines = []
            words = line.split(' ')
            line_acc = ''
            for word in words:
                if len(line_acc) + len(word) <= 50:
                    line_acc += word + ' '
                else:
                    lines.append(line_acc + '\n')
                    line_acc = '    ' + word + ' '
            lines.append(line_acc)
        self.prog_frm_lb.config(state=NORMAL)
        # print all of the lines
        for line in lines:
            self.prog_frm_lb.insert(END, line)
        self.prog_frm_lb.config(state=DISABLED)

    def _close(self):
        """
        handle all the possible conditions before closing the window
        """
        self.win.quit()
        self.win.destroy()
        
