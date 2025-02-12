"""
lipidoz/_util.py

Dylan Ross (dylan.ross@pnnl.gov)

    internal module with general utility functions
"""


from typing import Optional, Tuple, Union, Dict
from sqlite3 import connect
from multiprocessing import JoinableQueue, Process
import queue
import os
import io

from matplotlib import pyplot as plt, image as mpimage
import lzf
import numpy as np
from mzapy import MZA
from mzapy.isotopes import monoiso_mass


# type annotation that covers both MZA and UIMF readers
type OzData = Union[MZA, CustomUReader]

# molecular formula: dict mapping atom to count
type Formula = Dict[str, int]


def _polyunsat_ald_crg_formula(precursor_formula: Formula, 
                               db_position: int, 
                               db_idx: int
                               ) -> Tuple[Formula, Formula] :
    """
    produces the molecular formulas for aldehyde and criegee fragments corresponding to a polyunsaturated
    lipid species with double bond at the specified position (numbered from the end of the FA, i.e., n-X)
    db_idx corresponds to the index of the double bond, numbered from the end (e.g, a db_idx of 2 would
    indicate that the aldehyde and criegee fragments result from cleavage of the second double bound as
    counted from the end of the FA)

    Parameters
    ----------
    precursor_formula
        molecular formula as a dictionary mapping elements (str) to their counts (int)
    db_position 
        double bond position, numbered from the end of the FA
    db_idx 
        index of double bond, numbered from the end of the FA

    Returns
    -------
    aldehyde, criegee
        aldehyde and criegee fragment molecular formulas as a dicts mapping elements (str) to their counts (int)
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

    Parameters
    ----------
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

    Parameters
    ----------
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


def new_lipidoz_results():
    """
    create a new lipidoz results dictionary with empty fields 
    """
    lipidoz_results = {
        "isotope_scoring_results": None,
        "preml_data": None,
        "ml_data": None,
        "ml_pred_lbls": None,
        "ml_pred_probs": None,
        "ml_params_file": None
    }
    return lipidoz_results


# TODO (Dylan Ross): This really needs to be pared down to a single class with only the functionality
#                    that is absolutely required for OzID data analysis. I copied this in here from 
#                    other existing code for a different purpose to quickly get things up and running 
#                    but there is a lot that is unecessary. 


# TODO: This needs proper attribution of the original UIMFReader code with license and all.

class _UReader():
    """
    Minimal object for reading data from UIMF files
    """

    def __init__(self, path):
        """
        _UReader initialized with path to UIMF file. Should be closed when done

        TODO (Dylan Ross): accomodate manually modifying the mass calibration parameters

        Parameters
        ----------
        path : ``str``
            path to UIMF file
        workers : ``int``
            number of worker processes for deconding and decompressing spectra
        """
        self._con = connect(path)
        self._cur = self._con.cursor()
        # store number of m/z bins
        self._n_bins = self._get_n_bins()
        # store number of frames
        self._n_frames = self._get_n_frames()
        # store the max number of scans per frame (Prescan_TOFPulses)
        self._max_frame_n_scans = self._get_max_frame_n_scans()
        # store the average duration of all IM frames (in milliseconds)
        self._avg_frame_duration = self._get_avg_frame_duration()
        # store array of all m/zs
        self._all_mzs = self.mzbin_to_mz(np.arange(self._n_bins))        

    def close(self):
        """ clean up tasks """
        # close DB connection
        self._con.close()

    def _get_n_bins(self):
        """ fetches NBins value from Global_Params table """
        return int(self._cur.execute("SELECT ParamValue FROM Global_Params WHERE ParamName='Bins'").fetchone()[0])
    
    def _get_n_frames(self):
        """ fetches N_frames value from Global_Parameters table """
        return int(self._cur.execute('SELECT NumFrames FROM Global_Parameters').fetchone()[0])

    def _get_max_frame_n_scans(self):
        """ fetch the max n_scans for all frames (Prescan_TOFPulses from Global_Parameters) """
        return int(self._cur.execute('SELECT Prescan_TOFPulses FROM Global_Parameters').fetchone()[0])

    def _get_avg_frame_duration(self):
        """ computes the average frame duration across all IM frames in milliseconds """
        durations = [_[0] for _ in self._cur.execute('SELECT Duration FROM Frame_Parameters').fetchall()]
        return np.mean(durations) * 1000  # multiply by 1000 for milliseconds 

    def scan_to_ms(self, scan):
        """ 
        converts scan index into scan time in milliseconds, 

        Parameters
        ----------
        scan : ``int``
            scan index, operation is broadcastable to numpy.ndarray(int) for multiple scans if desired

        Returns
        -------
        scan_time : ``float``
            scan time in milliseconds
        """
        return  self._avg_frame_duration * scan / self._max_frame_n_scans

    @staticmethod
    def _decompress(blob):
        """
        Use LZF to decompress the encoded intensity array

        Parameters
        ----------
        blob : ``bytes``
            compressed and encoded intensity array, directly from SQLite3 query

        Returns
        -------
        enc_array : ``numpy.array(int32)``
            encoded intensity array (as 32 bit ints)
        """
        l = len(blob)
        if l == 0:
            return np.array([], dtype=np.int32)
        else:
            return np.frombuffer(lzf.decompress(blob, 100 * l), dtype=np.int32)

    @staticmethod
    def _decode(enc_array):
        """
        Decodes the encoded intensity array (output from _decompress), produces array of mzbins and intensities

        TODO (Dylan Ross): What type of encoding is being used here? -> RLZE!

        Parameters
        ----------
        enc_array : ``numpy.array(int32)``
            encoded intensity array (as 32 bit ints)

        Returns
        -------
        mzb_i : ``list(tuple(int, int))``
            pairs of m/z bin and intensity values
        """
        mzb_i = []
        #prev_value = 0
        bin_idx = 0
        for value in enc_array:
            if value < 0:
                # negative values increment bin index
                bin_idx += -int(value)
            #elif value == 0 and prev_value < -1e10:
            #    # ?
            #    pass
            else:
                mzb_i.append((bin_idx, value))
                bin_idx += 1
                #prev_value = value
        return mzb_i
    
    @staticmethod
    def _decode_for_atd(enc_array, mzbin_min, mzbin_max):
        i_sum = 0
        bin_idx = 0
        for value in enc_array:
            if value < 0:
                # negative values increment bin index
                bin_idx += -int(value)
            else:
                if bin_idx >= mzbin_min:
                    if bin_idx <= mzbin_max:
                        i_sum += value
                        bin_idx += 1
                    else:
                        break
        return i_sum

    def _accum_spectra(self, qry):
        """
        backend method for accumulating spectra, takes a query then does the extraction and processing 
        and returns the intensities array. Enables using slightly different queries without repeating the 
        logic for the actual data extraction and processing.

        Parameters
        ----------
        qry : ``str``
            query for selecting spectra (by frame and/or scan)
        
        Returns
        -------
        mzs : ``numpy.array(float)``
            m/z values for summed mass spectrum
        intensities : ``numpy.array(int)``
            summed mass spectrum, indices are m/z bins
        """
        intensities = np.zeros(self._n_bins)
        for _, blob in self._cur.execute(qry).fetchall():
            for mzb, i in _UReader._decode(_UReader._decompress(blob)):
                intensities[mzb] += i
        return self._all_mzs, intensities

    def accum_spectra_allframes(self, 
                                scans=None, scan_min=None, scan_max=None, skip_frame_1=False):
        """ 
        sum together spectra from all frames and specified scans, returns array of summed intensities by m/z bins
        scans can be specified as a list or as min/max values for a range

        Parameters
        ----------
        scans : ``list(int)``, optional
            list of scans to accumulate spectra from
        scan_min, scan_max : ``int``, optional
            alternative to scans, specify min/max scan values to include a range of scans
            ignored if scans is provided or if only min or max is provided
        skip_frame_1 : ``bool``, default=False
            skip the first frame when accumulating spectra from multiple frames

        Returns
        -------
        intensities : ``numpy.array(int)``
            summed mass spectrum, indices are m/z bins
        """
        qry = "SELECT ScanNum, Intensities FROM Frame_Scans"
        if scans is not None:
            s = ','.join([str(_) for _ in scans])
            qry += " WHERE ScanNum IN ({})".format(s)
        elif scan_min is not None and scan_max is not None:
            qry += " WHERE ScanNum>={} AND ScanNum<={}".format(scan_min, scan_max)
        if skip_frame_1:
            where = "WHERE " if "WHERE" not in qry else "AND "
            qry += " " + where + "FrameNum>1"
        return self._accum_spectra(qry)
    
    def extract_atd_allframes(self,
                              target_mz, mz_ppm,
                              skip_frame_1=False):
        """
        Extract an arrival time distribution for a target m/z +/- ppm
        """
        # compute the m/z tolerance from ppm
        tol = target_mz * mz_ppm / 1e6
        # figure out the min/max mz bins to include
        mzbin_min = self.mz_to_mzbin(target_mz - tol)
        mzbin_max = self.mz_to_mzbin(target_mz + tol)
        # set up the query to get the data
        qry = "SELECT ScanNum, Intensities FROM Frame_Scans" 
        if skip_frame_1:
            qry += " WHERE FrameNum>1"
        # create the ATD arrays to accumulate into
        #atd_at = self.scan_to_ms(np.arange(self._max_frame_n_scans))
        atd_i = np.zeros(self._max_frame_n_scans)
        for scan, blob in self._cur.execute(qry):
            atd_i[scan] += _UReader._decode_for_atd(_UReader._decompress(blob), mzbin_min, mzbin_max)
        return self.scan_to_ms(np.arange(self._max_frame_n_scans)), atd_i        
        
    def _get_mz_cal_params(self, frames):
        """ get m/z calibration parameters for a specified range of frames, average them together """
        qry = "SELECT ParamValue FROM V_Frame_Params WHERE ParamName='{}'"
        slope_qry = qry.format('CalibrationSlope')
        intercept_qry = qry.format('CalibrationIntercept')
        if frames is not None:
            f = ','.join([str(_) for _ in frames])
            slope_qry += " AND FrameNum IN ({})".format(f)
            intercept_qry += " AND FrameNum IN ({})".format(f)
        slopes = [float(_[0]) for _ in self._cur.execute(slope_qry).fetchall()]
        intercepts = [float(_[0]) for _ in self._cur.execute(intercept_qry).fetchall()]
        return np.mean(slopes), np.mean(intercepts)
    
    def _get_mz_cal_params_all_frames(self):
        """ get m/z calibration parameters for all frames, average them together """
        # set frames to None to include all frames
        return self._get_mz_cal_params(None)

    def _get_bin_width(self):
        """ fetches BinWidth value from GlobalParams table """
        qry = "SELECT ParamValue FROM Global_Params WHERE ParamName='BinWidth'"
        return float(self._cur.execute(qry).fetchone()[0])
    
    def mzbin_to_mz(self, mzbin):
        """
        convert mzbin to m/z using TOF calibration parameters

        Parameters
        ----------
        mzbin : ``int``
            mzbin

        Returns
        -------
        mz : ``float``
            m/z
        """
        bin_width = self._get_bin_width()
        cal_slope, cal_intercept = self._get_mz_cal_params_all_frames()
        return (cal_slope / 1e4 * (mzbin * bin_width * 10. - cal_intercept * 1e4))**2.
    
    def mz_to_mzbin(self, mz):
        # this is really terrible, there is a much better way to do this to be sure
        return np.argmin(np.abs(self._all_mzs - mz))


class CustomUReader(_UReader):
    """ custom UIMF reader to deal with these particular converted SLIM data files """

    def __init__(self, path):
        super().__init__(path)
        # this holds the arrays used for indexing out XICs and MS1 spectra
        # set by self.accum_frame_spectra_allscans method
        self.__frame_spectra_allscans = None

    def _get_n_bins(self):
        """ OVERRIDE: fetches NBins value from Global_Parameters table """
        return int(self._cur.execute("SELECT Bins FROM Global_Parameters").fetchone()[0])

    def _get_bin_width(self):
        """ OVERRIDE: fetches BinWidth value from GlobalParameters table """
        return float(self._cur.execute("SELECT BinWidth FROM Global_Parameters").fetchone()[0])

    def _get_mz_cal_params(self, frames):
        """ OVERRIDE: get m/z calibration parameters for a specified range of frames, average them together """
        qry = "SELECT CalibrationSlope, CalibrationIntercept FROM Frame_Parameters"
        if frames is not None:
            qry += f"WHERE FrameNum in ({','.join([str(_) for _ in frames])})"
        slopes, intercepts = np.array(self._cur.execute(qry).fetchall()).T
        return np.mean(slopes), np.mean(intercepts)
    
    def _get_max_frame_n_scans(self):
        """ OVERRIDE fetch the max n_scans for all frames (from Frame_Parameters) """
        return int(self._cur.execute('SELECT MAX(Scans) FROM Frame_Parameters').fetchone()[0])

    def _accum_frame_spectra_allscans(self, skip_frame_1):
        """ 
        sum together spectra from all scans (effectively ignoring arrival time dimension)
        separately for each frame and keep track of frame RT. Returns an array of 
        retention times for all frames, the full m/z array and a transposed list of 
        mass spectrum intensities that can be indexed by slicing the full m/z array. 

        .. note::

            Run this once to initialize the arrays, then it should be very fast to index out
            XICs using m/z selection. It is also possible to index with RT range and extract
            the m/z dimension to get MS1 spectrum. 

        .. note::

            This returns a dict with the arrays described in the Returns section below
            packed into it (keys are just the array variable names as strings)

        Parameters
        ----------
        skip_frame_1 : ``bool``
            skip the first frame when accumulating spectra

        Returns 
        -------
        frame_rts : ``numpy.ndarray(float)``
            Retention times, shape = (n_frames,)
        spectra : ``numpy.array(int)``
            2D array of summed mass spectrum intensities, shape = (n_mz_bins, n_frames)
            transposed from a list of mass spectrum intensities such that indexing by 
            the first dimension (m/z bins) produces the corresponding XIC
        """
        qry1 = "SELECT FrameNum, StartTime FROM Frame_Parameters"
        if skip_frame_1:
            qry1 += " WHERE FrameNum>1"
        qry2 = "SELECT ScanNum, Intensities FROM Frame_Scans WHERE FrameNum={}"
        frames, rts, spectra = [], [], []
        for frame, rt in self._cur.execute(qry1).fetchall():
            #print(f"{frame=}")
            frames.append(frame)
            rts.append(rt)
            # just take the intensities from the summed spectra
            spectra.append(self._accum_spectra(qry2.format(frame))[1])
        rts = np.array(rts)
        spectra = np.array(spectra).T
        return {"rts": rts, "spectra": spectra}
    
    def accum_frame_spectra_allscans(self, skip_frame_1):
        """
        Use self._accume_frame_spectra_allscans method to generate arrays of data with 
        IM dimension effectively collapsed that can be used to quickly index and select
        out XICs and MS1 spectra. 
        
        These arrays are stored in a dict in self._frame_spectra_allscans
        """
        if self.__frame_spectra_allscans is None:
            self.__frame_spectra_allscans = self._accum_frame_spectra_allscans(skip_frame_1)

    def collect_xic_arrays_by_mz(self, 
                                 mz_min: float, 
                                 mz_max: float, 
                                 rt_bounds: Optional[Tuple[float, float]] = None, 
                                 mslvl: int = 1, 
                                 verbose: bool = False) -> None :
        """
        """
        # this method requires the data from self.__frame_spectra_allscans
        if self.__frame_spectra_allscans is None:
            msg = ("CustomUReader: collect_xic_arrays_by_mz: "
                   "self.__frame_spectra_allscans is not set, call "
                   "self.accum_frame_spectra_allscans before calling this method")
            raise RuntimeError(msg)
        assert mslvl == 1, "MS levels other than MS1 not implemented"
        mz_idx = (self._all_mzs >= mz_min) & (self._all_mzs <= mz_max)
        xic_i = self.__frame_spectra_allscans["spectra"][mz_idx].sum(axis=0)
        xic_rt = self.__frame_spectra_allscans["rts"]
        if rt_bounds is not None:
            rt_min, rt_max = rt_bounds
            rt_idx = (self.__frame_spectra_allscans["rts"] >= rt_min) & (self.__frame_spectra_allscans["rts"] <= rt_max)
            xic_i = xic_i[rt_idx]
            xic_rt = xic_rt[rt_idx]
        return xic_rt, xic_i        

    def collect_ms1_arrays_by_rt(self, 
                                 rt_min: float, 
                                 rt_max: float, 
                                 mz_bounds: Optional[Tuple[float, float]] = None) -> None :
        """
        """
        # this method requires the data from self.__frame_spectra_allscans
        if self.__frame_spectra_allscans is None:
            msg = ("CustomUReader: collect_ms1_arrays_by_rt: "
                   "self.__frame_spectra_allscans is not set, call "
                   "self.accum_frame_spectra_allscans before calling this method")
            raise RuntimeError(msg)
        rt_idx = (self.__frame_spectra_allscans["rts"] >= rt_min) & (self.__frame_spectra_allscans["rts"] <= rt_max)
        ms1_mz = self._all_mzs
        ms1_i = self.__frame_spectra_allscans["spectra"].T[rt_idx].sum(axis=0)
        if mz_bounds is not None:
            mz_min, mz_max = mz_bounds
            mz_idx = (self._all_mzs >= mz_min) & (self._all_mzs <= mz_max)
            ms1_mz = ms1_mz[mz_idx]
            ms1_i = ms1_i[mz_idx]
        return ms1_mz, ms1_i

