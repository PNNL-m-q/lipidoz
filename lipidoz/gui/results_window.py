""" 
lipidoz_gui/results_window.py

Dylan Ross (dylan.ross@pnnl.gov)

    results window component of GUI for LipidOz utility

"""


from tkinter import *
from tkinter import ttk, filedialog, messagebox
import io
import os
import traceback

from PIL import Image, ImageTk
import numpy as np
from matplotlib.pyplot import cm

from lipidoz import __version__ as LOZ_VER
from lipidoz.workflows import save_isotope_scoring_results, write_isotope_scoring_report_xlsx
from lipidoz.gui._config import PLOT_HEIGHT, XIC_PLOT_WIDTH, MS1_PLOT_WIDTH, ICO_PATH, SAT_CORR_COLOR


class ResultsWindow():
    """
    """

    def __init__(self, isotope_scoring_results):
        """
        setup the window then start the mainloop
        """
        # status flag
        self.back_to_setup = False
        # store the results
        self.results = isotope_scoring_results
        #self._remove_nones_from_results()
        # init window and configure
        self.win = Tk()
        self._configure_window()
        # init and configure the main frame that will contain all other widgets
        self._setup_main_frame()
        # start up the window main loop
        self.win.mainloop()

    def _configure_window(self):
        # configure window
        self.win.title('LipidOz ({}) | Results'.format(LOZ_VER))
        # set the window icon
        self.win.iconbitmap(ICO_PATH)
        # set up window layout cols and rows
        self.win.rowconfigure(0, weight=1)
        self.win.columnconfigure(0, weight=1)

    def _setup_main_frame(self):
        """
        main frame has 3 rows and 1 column
            row 0 -> metadata frame
            row 1 -> content frame
            row 2 -> button frame
        """
        self.main_frm = ttk.Frame(self.win, padding=(10, 10, 10, 10))
        self.main_frm.grid(row=0, column=0, sticky=(N, S, E, W))
        # set up the main frame's grid structure
        self.main_frm.rowconfigure(0)  # label above the main content pane
        self.main_frm.rowconfigure(1, weight=1)  # main content pane
        self.main_frm.rowconfigure(2)  # bottom container for buttons
        self.main_frm.columnconfigure(0, minsize=400, weight=1)
        # set up the each of the sub-frames
        self._setup_metadata_frame()
        self._setup_content_frame()
        self._setup_button_frame()
        # callback for pressing backspace or delete key on selected item in treeview menu
        self.main_frm.bind_all('<BackSpace>', self._tree_item_delete_cb)
        self.main_frm.bind_all('<Delete>', self._tree_item_delete_cb)

    def _setup_metadata_frame(self):
        """
        create and populate the top frame position with experiment metadata
        """
        # NOT USED FOR NOW...

    def _setup_content_frame(self):
        """
        content frame has 1 row and 2 columns
            column 0 -> browseable menu of lipid species
            column 1 -> plots of lipid OzID results 
        """
        self.content_frm = ttk.Frame(self.main_frm)
        self.content_frm.grid(row=1, column=0, sticky=(N, S, E, W))
        # content frame is split side by side
        # left panel contains the treeview menu
        # right panel displays plots and other info
        self.content_frm.rowconfigure(0, weight=1)
        self.content_frm.columnconfigure(0, weight=1)
        self.content_frm.columnconfigure(1, weight=1)
        # setup left treeview menu frame
        self._setup_menu_frame()
        # setup the right pane of the content frame (plots frame)
        self._setup_plots_frame()

    def _setup_button_frame(self):
        """
        setup frame that contains Save Results and Export Results buttons
        """
        self.butt_frm = ttk.Frame(self.main_frm, padding=(5, 5, 5, 5))
        self.butt_frm.grid(row=2, column=0, sticky=(W, E))
        self.butt_frm.columnconfigure(0)
        self.butt_frm.columnconfigure(1, weight=1)
        self.butt_frm.columnconfigure(2)
        self.butt_frm_back_btn = ttk.Button(self.butt_frm, text='Back to Setup', command=self._back_btn_cb)
        self.butt_frm_back_btn.grid(row=0, column=0, sticky=(W,))
        self.butt_frm_saveres_btn = ttk.Button(self.butt_frm, text='Save Results', command=self._saveres_btn_cb)
        self.butt_frm_saveres_btn.grid(row=0, column=1, sticky=(E,))
        self.butt_frm_expres_btn = ttk.Button(self.butt_frm, text='Export Results', command=self._expres_btn_cb)
        self.butt_frm_expres_btn.grid(row=0, column=2, sticky=(E,))

    def _setup_menu_frame(self):
        """
        create and configure the menu frame which lives inside the content frame and holds the
        treeview display of lipid species
        """
        self.menu_frm = ttk.Frame(self.content_frm)
        self.menu_frm.grid(row=0, column=0, sticky=(N, S, E, W))
        self.menu_frm.rowconfigure(0, weight=1)
        self.menu_frm.columnconfigure(0, weight=1)
        self.menu_frm.columnconfigure(1)
        self.menu_frm_scroll = ttk.Scrollbar(self.menu_frm, orient='vertical')
        self.menu_frm_scroll.grid(row=0, column=1, sticky=(N, S))
        self.menu_frm_tv = ttk.Treeview(self.menu_frm, selectmode='browse', padding=(5, 5, 5, 5), show='tree')
        self.menu_frm_tv.configure(yscroll=self.menu_frm_scroll.set)
        self.menu_frm_tv.grid(row=0, column=0, sticky=(N, S, E, W))
        self.menu_frm_scroll.config(command=self.menu_frm_tv.yview)
        # callback for selecting an item in the treeview menu
        self.menu_frm_tv.bind('<<TreeviewSelect>>', self._tree_select_cb)
        # add all the lipid data to the treeview menu
        self._add_lipid_data_to_tree()

    def _setup_plots_frame(self):
        """
        create and configure the plots frame which lives inside the content frame and holds the
        plots of lipid OzID results

        plots frame has 8 rows and 2 columns
            row 0, column 0 -> "precursor xic"
            row 0, column 1 -> "precursor isotope distribution"
            row 1, column 0 -> precursor chromatogram fit plot
            row 1, column 1 -> precursor isotope distribution plot
            row 2, column 0 -> "aldehyde xic"
            row 2, column 1 -> "aldehyde isotope distribution"
            row 3, column 0 -> aldehyde chromatogram fit plot
            row 3, column 1 -> aldehyde isotope distribution plot
            row 4, column 0 -> "criegee xic"
            row 4, column 1 -> "criegee isotope distribution"
            row 5, column 0 -> criegee chromatogram fit plot
            row 5, column 1 -> criegee isotope distribution plot
            row 6, column 0 -> "scores"
            row 7, column 0 -> scores
        """
        self.plots_frm = ttk.Frame(self.content_frm)
        self.plots_frm.grid(row=0, column=1, sticky=(N, S, E, W))
        # set figure dimensions
        self.ht = PLOT_HEIGHT
        self.xwt, self.iwt = XIC_PLOT_WIDTH, MS1_PLOT_WIDTH
        # set up the grid
        self.plots_frm.columnconfigure(0, weight=1, minsize=self.xwt)
        self.plots_frm.columnconfigure(1, weight=1, minsize=self.iwt)
        self.plots_frm.rowconfigure(0, weight=0)
        self.plots_frm.rowconfigure(1, weight=1, minsize=self.ht)
        self.plots_frm.rowconfigure(2, weight=0)
        self.plots_frm.rowconfigure(3, weight=1, minsize=self.ht)
        self.plots_frm.rowconfigure(4, weight=0)
        self.plots_frm.rowconfigure(5, weight=1, minsize=self.ht)
        self.plots_frm.rowconfigure(6, weight=0)
        self.plots_frm.rowconfigure(7, weight=0)
        # precursor section
        self.plots_frm_pre_lbl1 = ttk.Label(self.plots_frm, text='Precursor XIC')
        self.plots_frm_pre_lbl1.grid(row=0, column=0, sticky=(W,))
        self.plots_frm_pre_lbl2 = ttk.Label(self.plots_frm, text='Precursor Isotope Distribution')
        self.plots_frm_pre_lbl2.grid(row=0, column=1, sticky=(W,))
        self.plots_frm_pre_xic = Canvas(self.plots_frm, width=self.xwt, height=self.ht)
        self.plots_frm_pre_xic.grid(row=1, column=0, sticky=(N, S, E, W), padx=10, pady=10)
        self.plots_frm_pre_xic.bind('<Double-Button-1>', self._pre_xic_dblclk)
        self.plots_frm_pre_isodist = Canvas(self.plots_frm, width=self.iwt, height=self.ht)
        self.plots_frm_pre_isodist.grid(row=1, column=1, sticky=(N, S, E, W), padx=10, pady=10)
        self.plots_frm_pre_isodist.bind('<Double-Button-1>', self._pre_isodist_dblclk)
        # aldehyde section
        self.plots_frm_ald_lbl1 = ttk.Label(self.plots_frm, text='Aldehyde XIC')
        self.plots_frm_ald_lbl1.grid(row=2, column=0, sticky=(W,))
        self.plots_frm_ald_lbl2 = ttk.Label(self.plots_frm, text='Aldehyde Isotope Distribution')
        self.plots_frm_ald_lbl2.grid(row=2, column=1, sticky=(W,))
        self.plots_frm_ald_xic = Canvas(self.plots_frm, width=self.xwt, height=self.ht)
        self.plots_frm_ald_xic.grid(row=3, column=0, sticky=(N, S, E, W), padx=10, pady=10)
        self.plots_frm_ald_xic.bind('<Double-Button-1>', self._ald_xic_dblclk)
        self.plots_frm_ald_isodist = Canvas(self.plots_frm, width=self.iwt, height=self.ht)
        self.plots_frm_ald_isodist.grid(row=3, column=1, sticky=(N, S, E, W), padx=10, pady=10)
        self.plots_frm_ald_isodist.bind('<Double-Button-1>', self._ald_isodist_dblclk)
        # criegee section
        self.plots_frm_crg_lbl1 = ttk.Label(self.plots_frm, text='Criegee XIC')
        self.plots_frm_crg_lbl1.grid(row=4, column=0, sticky=(W,))
        self.plots_frm_crg_lbl2 = ttk.Label(self.plots_frm, text='Criegee Isotope Distribution')
        self.plots_frm_crg_lbl2.grid(row=4, column=1, sticky=(W,))
        self.plots_frm_crg_xic = Canvas(self.plots_frm, width=self.xwt, height=self.ht)
        self.plots_frm_crg_xic.grid(row=5, column=0, sticky=(N, S, E, W), padx=10, pady=10)
        self.plots_frm_crg_xic.bind('<Double-Button-1>', self._crg_xic_dblclk)
        self.plots_frm_crg_isodist = Canvas(self.plots_frm, width=self.iwt, height=self.ht)
        self.plots_frm_crg_isodist.grid(row=5, column=1, sticky=(N, S, E, W), padx=10, pady=10)
        self.plots_frm_crg_isodist.bind('<Double-Button-1>', self._crg_isodist_dblclk)
        # clear out all of the plots
        self._clear_plots_frame()
        # scores section
        #self.plots_frm_scores_lbl = ttk.Label(self.plots_frm, text='Scores')
        #self.plots_frm_scores_lbl.grid(row=6, column=0, sticky=(W,))
        self._setup_scores_frm()

    def _setup_scores_frm(self):
        """
        """
        self.scores_frm = ttk.Frame(self.plots_frm)
        self.scores_frm.grid(row=7, column=0, sticky=(N, S, E, W), columnspan=2)
        # setup grid
        self.scores_frm.columnconfigure(0, weight=0, minsize=75)
        self.scores_frm.columnconfigure(1, weight=0, minsize=75)
        self.scores_frm.columnconfigure(2, weight=0, minsize=75)
        self.scores_frm.columnconfigure(3, weight=1, minsize=150)
        self.scores_frm.rowconfigure(0, weight=0)
        self.scores_frm.rowconfigure(1, weight=0)
        self.scores_frm.rowconfigure(2, weight=0)
        # add  labels
        self.scores_frm_ald_lbl = ttk.Label(self.scores_frm, text='Aldehyde')
        self.scores_frm_ald_lbl.grid(row=0, column=1)
        self.scores_frm_crg_lbl = ttk.Label(self.scores_frm, text='Criegee')
        self.scores_frm_crg_lbl.grid(row=0, column=2)
        self.scores_frm_mzcos_lbl = ttk.Label(self.scores_frm, text='m/z cosine distance')
        self.scores_frm_mzcos_lbl.grid(row=1, column=0, sticky=(E,))
        self.scores_frm_rtcos_lbl = ttk.Label(self.scores_frm, text='RT cosine distance')
        self.scores_frm_rtcos_lbl.grid(row=2, column=0, sticky=(E,))
        self.scores_frm_mzcos_ald_lbl = ttk.Label(self.scores_frm)
        self.scores_frm_mzcos_ald_lbl.grid(row=1, column=1)
        self.scores_frm_mzcos_crg_lbl = ttk.Label(self.scores_frm)
        self.scores_frm_mzcos_crg_lbl.grid(row=1, column=2)
        self.scores_frm_rtcos_ald_lbl = ttk.Label(self.scores_frm)
        self.scores_frm_rtcos_ald_lbl.grid(row=2, column=1)
        self.scores_frm_rtcos_crg_lbl = ttk.Label(self.scores_frm)
        self.scores_frm_rtcos_crg_lbl.grid(row=2, column=2)
        # add a legend explaining how to interpret plots
        # TODO: needs updating
        # legend_text = ('XICs: grey trace = raw data, red dashed line and shaded area = target RT with RT '
        #                'extraction window,\n        blue crossing lines = fitted peak(s), blue shaded region = '
        #                'extraction window\nMS1 spectra: grey trace = raw data, blue crossing lines = '
        #                'fitted isotopic peaks,\n                red crossing lines = theoretical isotopic peaks')
        legend_text = ""
        self.scores_frm_legend_lbl = ttk.Label(self.scores_frm, text=legend_text, relief=GROOVE)
        self.scores_frm_legend_lbl.grid(row=0, column=3, rowspan=3, sticky=(N, S, E, W))

    def _add_lipid_data_to_tree(self):
        """
        iterate over the lipids in the results file and add them into the treeview
        """
        # clear out menu if it has already been filled
        for i in self.menu_frm_tv.get_children():
            self.menu_frm_tv.delete(i)
        # map iid to tuples with all the indexing information for all leaf nodes in the menu
        self._tree_leaf_nodes = {}
        # map sets of square images to idbidx
        self._tree_square_imgs = {}
        # build the tree
        for lipid in self.results['targets']:
            ilipid = self.menu_frm_tv.insert('', END, text=f"{lipid}", open=False)
            ilipid_children = 0
            for adduct in self.results['targets'][lipid]:
                iadduct = self.menu_frm_tv.insert("", END, text=f"{adduct}", open=False)
                iadduct_children = 0
                ilipid_children += 1
                self.menu_frm_tv.move(iadduct, ilipid, ilipid_children)
                for rt in self.results['targets'][lipid][adduct]:
                    if self.results['targets'][lipid][adduct][rt] is not None:
                        irt = self.menu_frm_tv.insert('', END, text=f"{rt}", open=False)
                        irt_children = 0
                        iadduct_children += 1
                        self.menu_frm_tv.move(irt, iadduct, iadduct_children)
                        for db_idx in self.results['targets'][lipid][adduct][rt]['fragments']:
                            idbidx = self.menu_frm_tv.insert("", END, text=f"{db_idx=}", open=False)
                            idbidx_children = 0
                            irt_children += 1
                            self.menu_frm_tv.move(idbidx, irt, irt_children)
                            # sort by ascending composite score
                            db_posns = [_ for _ in self.results['targets'][lipid][adduct][rt]['fragments'][db_idx].keys()]
                            comp_scores = [self._get_composite_score(self.results['targets'][lipid][adduct][rt]['fragments'][db_idx][_]) for _ in db_posns]
                            self._gen_colored_square_imgs(comp_scores, idbidx)
                            for sort_i in np.argsort(comp_scores):
                                db_pos = db_posns[sort_i]
                                fragment_data = self.results['targets'][lipid][adduct][rt]['fragments'][db_idx][db_pos]
                                if fragment_data['aldehyde'] is not None or fragment_data['criegee'] is not None:
                                    # only add entries when there is actual data to look at
                                    idbpos = self.menu_frm_tv.insert("", END, text=f"{db_pos=:<2d}",
                                                                     open=False,
                                                                     image=self._tree_square_imgs[idbidx][sort_i])
                                    idbidx_children += 1
                                    self.menu_frm_tv.move(idbpos, idbidx, idbidx_children)
                                    # add an entry to self._tree_leaf_nodes
                                    idx_info = (lipid, adduct, rt, db_idx, db_pos)
                                    self._tree_leaf_nodes[idbpos] = idx_info

    def _gen_colored_square_imgs(self, comp_scores, idbidx):
        """ """
        cmap = cm.get_cmap('coolwarm')
        # map composite scores from 0 -> 0.265 * 2 into the range 0.2 -> 0.8, reversed
        def get_color(comp_score):
            remapped = comp_score / (0.265 * 2)
            rgb = cmap(remapped)[:3]
            return tuple([int(255 * _) for _ in rgb])
        cs = [get_color(score) for score in comp_scores]
        self._tree_square_imgs[idbidx] = [ImageTk.PhotoImage(Image.new('RGB', (10, 10), c)) for c in cs]

    def _tree_select_cb(self, event):
        """
        callback that gets called whenever the selection on the tree menu changes
        """
        if len(sel := self.menu_frm_tv.selection()) == 1:
            sel_idx = sel[0]   
            if sel_idx in self._tree_leaf_nodes:
                # an individual DB position result is selected
                sel_idx_info = self._tree_leaf_nodes[sel_idx]
                # update the plots page with the results from the selection
                self._populate_plots_frame(sel_idx_info)
                self._update_scores(sel_idx_info)
            else:
                # an intermediate value is selected
                # clear the display
                self._clear_plots_frame()
                self._clear_scores()

    def _tv_get_sib_or_parent_idx(self, idx):
        """ get index of a sibling node or parent node if there are no siblings """
        if (sib_idx := self.menu_frm_tv.prev(idx)) != "":
            return sib_idx
        if (sib_idx := self.menu_frm_tv.next(idx)) != "":
            return sib_idx
        return self.menu_frm_tv.parent(idx)

    def _tv_focus_and_select(self, idx):
        """ move focus and selection to specified index in TreeView """
        self.menu_frm_tv.focus(idx)
        self.menu_frm_tv.selection_set(idx)
    
    def _tv_move_up_and_detach_old_level(self, sel_idx):
        """ move focus one level up and detach the old lower level """
        old_idx = sel_idx
        sel_idx = self.menu_frm_tv.parent(old_idx)
        self._tv_focus_and_select(sel_idx)
        # now detach the old level
        self.menu_frm_tv.detach(old_idx)
        # return the new selection index (one level up)
        return sel_idx

    def _tree_item_delete_cb(self, event):
        """
        callback that gets called when delete key pressed
        """
        # TODO: Delete other nodes besides leaf nodes?
        if len(sel := self.menu_frm_tv.selection()) == 1:
            sel_idx = sel[0]
            if sel_idx in self._tree_leaf_nodes:
                # the selected node is a leaf node
                sel_idx_info = self._tree_leaf_nodes[sel_idx]
                lipid, adduct, rt, db_idx, db_pos = sel_idx_info
                # move focus to another sibling and select it
                old_idx = sel_idx
                sel_idx = self._tv_get_sib_or_parent_idx(old_idx)
                self._tv_focus_and_select(sel_idx)
                self.menu_frm_tv.detach(old_idx)
                # remove the selected entry from results
                self.results['targets'][lipid][adduct][rt]['fragments'][db_idx].pop(db_pos)
                # run checks to see if we need to prune back 
                # prune to db_idx level
                if self.results["targets"][lipid][adduct][rt]["fragments"][db_idx] == {}:
                    self.results["targets"][lipid][adduct][rt]["fragments"].pop(db_idx)
                    # after the previous leaf node was removed above, the sel_idx should have
                    # been moved to the container db_idx index. Move focus one level up and 
                    # delete the previous level
                    sel_idx = self._tv_move_up_and_detach_old_level(sel_idx)
                # prune to RT level
                if self.results["targets"][lipid][adduct][rt]["fragments"] == {}:
                    self.results["targets"][lipid][adduct].pop(rt)
                    # after the previous level was removed, the sel_idx should have 
                    # been moved to the container RT level. Move focus one level up and 
                    # delete the previous level
                    sel_idx = self._tv_move_up_and_detach_old_level(sel_idx)
                # prune to adduct level
                if self.results["targets"][lipid][adduct] == {}:
                    self.results["targets"][lipid].pop(adduct)
                    # after the previous level was removed, the sel_idx should have 
                    # been moved to the container adduct level. Move focus one level up and 
                    # delete the previous level
                    sel_idx = self._tv_move_up_and_detach_old_level(sel_idx)
                # prune to lipid level
                if self.results["targets"][lipid] == {}:
                    self.results["targets"].pop(lipid)
                    # after the previous level was removed, the sel_idx should have 
                    # been moved to the container lipid level. Move focus one level up and 
                    # delete the previous level
                    sel_idx = self._tv_move_up_and_detach_old_level(sel_idx)
            # else: the selected node was an intermediate level

    def _populate_plots_frame(self, idx_info):
        """
        add plots to the plot frame from the results corresponding to indexing info
        """
        # NOTE: It is possible for ald_data to not be None but one (or both?) of the actual XIC or 
        #       isotope distribution images to be None. In such cases, simply checking whether 
        #       ald_data is None causes an attempt to load the None as an image which leads to weird
        #       behavior where the images do not clear fully when switching between targets. The
        #       same applies for crg_data. The solution is two part: (1) always clear out all of 
        #       the plots at the beginning of this method and (2) add a second layer of checking 
        #       within each of the blocks that deal with populating plots from the ald_data or
        #       crg_data.
        # clear out any existing plots
        self._clear_plots_frame()
        # indexing information
        lipid, adduct, rt, db_idx, db_pos = idx_info
        # TODO: Perform some indexing here at the top for the precursor section and the 
        #       selected fragment section to avoid these long repetitive multi-part 
        #       indexing operations to grab different pieces of information to display.
        # fetch precursor image data
        resample = Image.BICUBIC
        img_data = self.results['targets'][lipid][adduct][rt]['precursor']['xic_fit_img']
        self.pre_xic_img = Image.open(io.BytesIO(img_data))
        self.plots_frm_pre_xic_img = ImageTk.PhotoImage(self.pre_xic_img.resize((self.xwt, self.ht), resample))
        self.plots_frm_pre_xic.create_image(0, 0, anchor=NW, image=self.plots_frm_pre_xic_img)
        # if saturation correction was performed in RT fitting add a little green square to the plot
        if self.results['targets'][lipid][adduct][rt]['precursor']['saturation_corrected']:
            # make an image to annotate precursor if saturation was corrected
            self.pre_sat_corr_img = ImageTk.PhotoImage(Image.new('RGB', (20, 20), SAT_CORR_COLOR))
            self.plots_frm_pre_xic.create_image(0, 0, image=self.pre_sat_corr_img)
        img_data = self.results['targets'][lipid][adduct][rt]['precursor']['isotope_dist_img']
        self.pre_isodist_img = Image.open(io.BytesIO(img_data))
        self.plots_frm_pre_isodist_img = ImageTk.PhotoImage(self.pre_isodist_img.resize((self.iwt, self.ht), resample))
        self.plots_frm_pre_isodist.create_image(0, 0, anchor=NW, image=self.plots_frm_pre_isodist_img)
        # fetch aldehyde image data (if present)
        if (ald_data := self.results['targets'][lipid][adduct][rt]['fragments'][db_idx][db_pos]['aldehyde']) is not None:
            # XIC
            if (xic_img := ald_data["xic_fit_img"]) is not None:
                self.ald_xic_img = Image.open(io.BytesIO(xic_img))
                self.plots_frm_ald_xic_img = ImageTk.PhotoImage(self.ald_xic_img.resize((self.xwt, self.ht), resample))
                self.plots_frm_ald_xic.create_image(0, 0, anchor=NW, image=self.plots_frm_ald_xic_img)
                if ald_data['saturation_corrected']:
                    # make an image to annotate precursor if saturation was corrected
                    self.ald_sat_corr_img = ImageTk.PhotoImage(Image.new('RGB', (20, 20), SAT_CORR_COLOR))
                    self.plots_frm_ald_xic.create_image(0, 0, image=self.ald_sat_corr_img)
            # isotope distribution
            if (iso_img := ald_data["isotope_dist_img"]) is not None:
                self.ald_isodist_img = Image.open(io.BytesIO(iso_img))
                self.plots_frm_ald_isodist_img = ImageTk.PhotoImage(self.ald_isodist_img.resize((self.iwt, self.ht), resample))
                self.plots_frm_ald_isodist.create_image(0, 0, anchor=NW, image=self.plots_frm_ald_isodist_img)
        # fetch criegee image data (if present)
        if (crg_data := self.results['targets'][lipid][adduct][rt]['fragments'][db_idx][db_pos]['criegee']) is not None:
            # XIC
            if (xic_img := crg_data["xic_fit_img"]) is not None:
                self.crg_xic_img = Image.open(io.BytesIO(xic_img))
                self.plots_frm_crg_xic_img = ImageTk.PhotoImage(self.crg_xic_img.resize((self.xwt, self.ht), resample))
                self.plots_frm_crg_xic.create_image(0, 0, anchor=NW, image=self.plots_frm_crg_xic_img)
                if crg_data['saturation_corrected']:
                    # make an image to annotate precursor if saturation was corrected
                    self.crg_sat_corr_img = ImageTk.PhotoImage(Image.new('RGB', (20, 20), SAT_CORR_COLOR))
                    self.plots_frm_crg_xic.create_image(0, 0, image=self.crg_sat_corr_img)
            # isotope distribution    
            if (iso_img := crg_data["isotope_dist_img"]) is not None:
                self.crg_isodist_img = Image.open(io.BytesIO(iso_img))
                self.plots_frm_crg_isodist_img = ImageTk.PhotoImage(self.crg_isodist_img.resize((self.iwt, self.ht), resample))
                self.plots_frm_crg_isodist.create_image(0, 0, anchor=NW, image=self.plots_frm_crg_isodist_img)

    def _clear_plots_frame(self):
        """
        clear all of the plots from the plots frame
        """
        # set all of the Image references to None (full-sized images)
        self.pre_xic_img = None
        self.pre_isodist_img = None
        self.ald_xic_img = None
        self.ald_isodist_img = None
        self.crg_xic_img = None
        self.crg_isodist_img = None
        # set all of the PhotoImage references to None
        self.plots_frm_pre_xic_img = None
        self.plots_frm_pre_isodist_img = None
        self.plots_frm_ald_xic_img = None
        self.plots_frm_ald_isodist_img = None
        self.plots_frm_crg_xic_img = None
        self.plots_frm_crg_isodist_img = None
        self.pre_sat_corr_img = None
        self.ald_sat_corr_img = None
        self.crg_sat_corr_img = None

    def _get_composite_score(self, target_data):
        """
        composite score is the average of m/z and rt cosine distances for ald and crg oz fragments, 
        missing values are set to 1
        """
        # fetch aldehyde scores (if present)
        ald_data = target_data['aldehyde'] if 'aldehyde' in target_data else None
        scores = [1, 1, 1, 1]
        if ald_data is not None:
            scores[0] = ald_data['mz_cos_dist']
            scores[1] = ald_data['rt_cos_dist']
        # fetch criegee scores (if present)
        crg_data = target_data['criegee'] if 'criegee' in target_data else None
        if crg_data is not None:
            scores[2] = crg_data['mz_cos_dist']
            scores[3] = crg_data['rt_cos_dist']
        return sum(scores) / 4.

    def _update_scores(self, idx_info):
        """
        """
        # indexing information
        lipid, adduct, rt, db_idx, db_pos = idx_info
        # fetch aldehyde image data (if present)
        ald_data = self.results['targets'][lipid][adduct][rt]['fragments'][db_idx][db_pos]['aldehyde']
        if ald_data is not None:
            self.scores_frm_mzcos_ald_lbl.configure(text='{:.4f}'.format(ald_data['mz_cos_dist']))
            self.scores_frm_rtcos_ald_lbl.configure(text='{:.4f}'.format(ald_data['rt_cos_dist']))
        else:
            self.scores_frm_mzcos_ald_lbl.configure(text='')
            self.scores_frm_rtcos_ald_lbl.configure(text='')
        # fetch criegee image data (if present)
        crg_data = self.results['targets'][lipid][adduct][rt]['fragments'][db_idx][db_pos]['criegee']
        if crg_data is not None:
            self.scores_frm_mzcos_crg_lbl.configure(text='{:.4f}'.format(crg_data['mz_cos_dist']))
            self.scores_frm_rtcos_crg_lbl.configure(text='{:.4f}'.format(crg_data['rt_cos_dist']))
        else:
            self.scores_frm_mzcos_crg_lbl.configure(text='')
            self.scores_frm_rtcos_crg_lbl.configure(text='')

    def _clear_scores(self):
        """
        """
        self.scores_frm_mzcos_ald_lbl.configure(text='')
        self.scores_frm_rtcos_ald_lbl.configure(text='')
        self.scores_frm_mzcos_crg_lbl.configure(text='')
        self.scores_frm_rtcos_crg_lbl.configure(text='')

    def _back_btn_cb(self):
        """
        callback for back to setup button, set status flag and exit
        """
        self.back_to_setup = True
        self._close()

    def _check_ext(self, fname, ext):
        """
        checks a file name's extension, and if it does not match ext, add ext to the filename and return that. 
        if the name already matches just return the filename
        """
        return fname + '.' + ext if os.path.splitext(fname)[-1] != ext else fname

    def _saveres_btn_cb(self):
        """
        """
        # the results are propted to be saved in the same directory as the data file
        oz_data_dir = os.path.abspath(os.path.split(self.results['metadata']['oz_data_file'])[0])
        out_f = filedialog.asksaveasfilename(title='Save LipidOz isotope scoring results',
                                             initialdir=oz_data_dir, 
                                             filetypes=[('isotope scoring results', '.loz')])
        if out_f is not None and out_f != '':
            try:
                # check for lack of file extension and fix if needed
                out_f = self._check_ext(out_f, 'loz')
                save_isotope_scoring_results(self.results, out_f)
            except Exception as e:
                # do a message box with the error
                msg = str(traceback.format_exc(limit=3)) + '\n' + str(e)
                messagebox.showerror(title='Unable to save results', message=msg)

    def _expres_btn_cb(self):
        """
        """
        # the exported results are prompted to save in the current working directory
        out_f = filedialog.asksaveasfilename(title='Export istope scoring results to spreadsheet',
                                             initialdir=os.getcwd(),
                                             filetypes=[('Excel Spreadsheet', '.xlsx')])
        if out_f is not None and out_f != '':
            try:
                # check for lack of file extension and fix if needed
                out_f = self._check_ext(out_f, 'xlsx')
                write_isotope_scoring_report_xlsx(self.results, out_f)
            except Exception as e:
                # do a message box with the error
                msg = str(traceback.format_exc(limit=3)) + '\n' + str(e)
                messagebox.showerror(title='Unable to export results', message=msg)

    def _pre_xic_dblclk(self, event):
        """
        """
        if self.pre_xic_img is not None:
            self._view_full_size_plot(self.pre_xic_img)

    def _pre_isodist_dblclk(self, event):
        """
        """
        if self.pre_isodist_img is not None:
            self._view_full_size_plot(self.pre_isodist_img)

    def _ald_xic_dblclk(self, event):
        """
        """
        if self.ald_xic_img is not None:
            self._view_full_size_plot(self.ald_xic_img)

    def _ald_isodist_dblclk(self, event):
        """
        """
        if self.ald_isodist_img is not None:
            self._view_full_size_plot(self.ald_isodist_img)

    def _crg_xic_dblclk(self, event):
        """
        """
        if self.crg_xic_img is not None:
            self._view_full_size_plot(self.crg_xic_img)

    def _crg_isodist_dblclk(self, event):
        """
        """
        if self.crg_isodist_img is not None:
            self._view_full_size_plot(self.crg_isodist_img)

    def _view_full_size_plot(self, img):
        pimg = ImageTk.PhotoImage(img)
        win = Toplevel(self.win)
        lbl = Label(win, image=pimg).pack()
        win.mainloop()

    def _close(self):
        """
        handle all the possible conditions before closing the window
        """
        self.win.quit()
        self.win.destroy()

