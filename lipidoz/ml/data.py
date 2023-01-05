"""
lipidoz/ml/data.py

Dylan Ross (dylan.ross@pnnl.gov)

    module with utilities for handling ML data 
"""


import pickle

from scipy import interpolate
import numpy as np
from torch.utils.data import Dataset, DataLoader
from torchvision import transforms
from sklearn.model_selection import ShuffleSplit

from mzapy.peaks import lerp_2d

from lipidoz._pyliquid import parse_lipid_name
from lipidoz._util import _debug_handler


def load_preml_data(preml_file):
    """ 
    loads a pre-ml dataset (produced by :func:`lipidoz.workflows.collect_pre_ml_dataset`) from specified file 

    Parameters
    ----------
    preml_file : ``str``
        path to pre-ml dataset file

    Returns
    -------
    preml_data : ``dict(...)``
        pre-ml dataset
    """
    with open(preml_file, 'rb') as pf:
        preml_data = pickle.load(pf)
    return preml_data


def load_ml_targets(ml_target_file, rt_corr_func=None):
    """
    loads list of annotated lipids from .csv formatted target list with columns:

        * lipid -- lipid name in standard format (``str``) 
        * adduct -- MS adduct (``str``)
        * retention time -- target retention time (``float``)
        * true_dbidx -- known double bond index (``int``)
        * true_dbopos -- known double bond position(s), separated by - if multiple at this index, gets unpacked into 
          a ``list(int)`` with all double bond positions at that index

    A list of lipid targets is returned, each defined by the information above, but grouped by lipid/adduct/retention
    time. Each lipid target contains:
    
        * lipid (``str``)
        * adduct (``str``)
        * retention time (rounded to 2 decimal places, ``float``)
        * annotations (``dict(int:list(int))``), a dict mapping db indices to lists of db positions

    Parameters
    ----------
    ml_target_file : ``str``
        path to target list .csv
    rt_corr_func : ``func``, optional
        apply rt correction using provided function

    Returns
    -------
    targets : ``dict(...)``
        dict mapping (lipid, adduct, retention time) to annotations (known db indices and positions)
    """
    target_data = {}
    with open(ml_target_file, 'r') as f:
        next(f)
        for line in f.readlines():
            lipid, adduct, retention_time, true_dbidx, true_dbpos = line.strip().split(',')
            # lipid name may need to be shuffled around a bit for consistency (specifically, chain order with _ sep.)
            lipid = str(parse_lipid_name(lipid))
            retention_time = rt_corr_func(float(retention_time)) if rt_corr_func is not None else float(retention_time)
            retention_time = round(retention_time, 2)
            true_dbidx = int(true_dbidx)
            true_dbpos = [int(_) for _ in true_dbpos.split('-')] if '-' in true_dbpos else [int(true_dbpos)]
            if (lipid, adduct, retention_time) in target_data:
                target_data[(lipid, adduct, retention_time)][true_dbidx] = true_dbpos
            else:
                target_data[(lipid, adduct, retention_time)] = {true_dbidx: true_dbpos}
    return target_data


def split_true_and_false_preml_data(preml_data, targets):
    """
    Splits a pre-ml dataset into true/false annotated examples based on a list of target lipids with annotated double
    bond positions and indices. The annotated double bond position and indices determine which entries are put into
    the true split for a given lipid target, the rest of the entries for that target are put into the false split. 
    
    Parameters
    ----------
    preml_data : ``dict(...)``
        a pre-ml dataset produced by :func:`lipidoz.workflows.collect_pre_ml_dataset`
    targets : ``dict(...)``
        dict mapping (lipid, adduct, retention time) to annotations (known db indices and positions), 
        produced by :func:`lipidoz.ml.data.load_ml_targets`

    Returns
    -------
    true_preml_data : ``dict(str:dict(...))``
    false_preml_data : ``dict(str:dict(...))``
        new pre-ml datasets with entries in 'targets' split by annotation (T/F) 
    """
    # set up empty split datasets
    true_preml_data = {}
    false_preml_data = {}
    # copy over metadata from the original pre-ml dataset
    true_preml_data['metadata'] = preml_data['metadata']
    false_preml_data['metadata'] = preml_data['metadata']
    # iterate through lipid targets and sort into true/false datasets 
    true_preml_data['targets'], false_preml_data['targets'] = {}, {}
    for target, annotation in targets.items(): 
        target_prefix = '{}|{}|{:.2f}min|'.format(*target)
        for k, v in preml_data['targets'].items():
            if target_prefix in k:
                db_idx, db_pos = [int(_) for _ in k.split('|')[-2:]]
                if (db_idx in annotation) and (db_pos in annotation[db_idx]):
                    true_preml_data['targets'][k] = v
                else:
                    false_preml_data['targets'][k] = v
    return true_preml_data, false_preml_data


def _get_binned_data_from_rtmz_arrays(rtmz, rt_bounds, mz_bounds, rt_bins=24, mz_bins=400, method='lerp2d'):
    """
    takes RTMZ array data ([retention time], [m/z], [intensity]) and returns a 2D array of data from specified 
    retention time and m/z bounds, binned at specified dimensions

    Default binning dimensions correspond to default sizes used for machine learning data

    Different methods can be used to covert the sparse scan data into continuous binned data:

        * "accumulate" - iterate over all of the rt bins and sum together all spectra within that range, will leave
                         gaps in the data when the rt binning is more closely spaced than the original scan sampling 
                         in the RT dimension
        * "lerp2d" - uses 2D linear interpolation across all datapoints to produce binned data, can cause some 
                     artifacting from gaps in the original sparse scan data

    Parameters
    ----------
    rtmz : ``tuple(numpy.ndarray(float))``
        arrays of retention time, m/z, and intensity components of RTMZ data
    rt_bounds : ``tuple(float)``
        lower, upper retention time bounds to bin data within
    mz_bounds : ``tuple(float)``
        lower, upper m/z bounds to bin data within
    rt_bins : ``int``, default=24
        number of bins to use for retention time dimension
    mz_bins : ``int``, default=400
        number of bins to use for m/z dimension
    method : ``str``, default='lerp2d'
        method to use for converting sparse scan data to continuous binned data, must be "accumulate" or "lerp2d"

    Returns
    -------
    binned_data : ``numpy.ndarray(float)``
        binned data array with shape (rt_bins, mz_bins)
    """
    rt_min, rt_max = rt_bounds
    rt_range = rt_max - rt_min
    rt_density = rt_bins / rt_range
    rt_binned = np.linspace(rt_min, rt_max, rt_bins)
    rt_bin_spacing = rt_binned[1] - rt_binned[0]
    mz_min, mz_max = mz_bounds
    mz_range = mz_max - mz_min
    mz_density = mz_bins / mz_range
    mz_binned = np.linspace(mz_min, mz_max, mz_bins)
    if method == 'accumulate':
        binned_data = []
        for rtb in range(rt_bins):
            # collect m/z, intensity data that falls within the bounds of this bin
            # rt values in rt_binned represent the bin centers so select values with RT in range: 
            # bin center +/- (spacing/2)
            mz, i = [], []
            for rt_, mz_, i_ in zip(*rtmz):
                if abs(rt_ - rt_binned[rtb]) <= rt_bin_spacing / 2.:
                    mz.append(mz_)
                    i.append(i_)
            if len(mz) < 2:
                # no point trying to do interpolation and stuff if there are not enough scans at this RT bin
                binned_data.append(np.zeros(mz_binned.shape))
            else:
                mz, i = np.array(mz), np.array(i)
                # sort before interpolation
                idx = np.argsort(mz)
                mz, i = mz[idx], i[idx]
                # interpolate mass spectrum, add to binned array
                bd = interpolate.interp1d(mz, i, kind='linear', fill_value=0, bounds_error=False)(mz_binned)
                binned_data.append(bd)
        return np.array(binned_data)
    elif method == 'lerp2d':
        # linear interpolation in 2D
        return lerp_2d(*rtmz, rt_min, rt_max, rt_density, mz_min, mz_max, mz_density, )[2].T
    else:
        msg = '_get_binned_data_from_rtmz_arrays: method must be "accumulate" or "lerp2d", was: "{}"'.format(method)
        raise ValueError(msg)


def _get_aug_rt_bounds(rt):
    """ returns a list of rt bounds with window size and target rt shifted for data augmentation """
    aug_rt_bounds = [
        [rt - 0.3, rt + 0.3],
        [rt - 0.25, rt + 0.25],
        [rt - 0.05 - 0.2, rt - 0.05 + 0.2],
        [rt - 0.025 - 0.2, rt - 0.025 + 0.2],
        [rt + 0.05 - 0.2, rt + 0.05 + 0.2],
        [rt + 0.025 - 0.2, rt + 0.025 + 0.2],
        [rt - 0.175, rt + 0.175],
        [rt - 0.025 - 0.175, rt - 0.025 + 0.175],
        [rt + 0.025 - 0.175, rt + 0.025 + 0.175],
        [rt - 0.15, rt + 0.15,],
        [rt - 0.225, rt + 0.225],
    ]
    return aug_rt_bounds


def preml_to_ml_data(preml_data, rt_sampling_augment=False, normalize_intensity=True,
                     debug_flag=None, debug_cb=None):
    """
    Takes pre-ml dataset and performs standard binning and grouping of precursor/fragment RTMZ profiles into arrays

    Parameters
    ----------
    preml_data : ``dict(...)``
        pre-ml dataset (produced by :func:`lipidoz.workflows.collect_pre_ml_dataset`)
    rt_sampling_augment : ``bool``, default=False
        re-sample RT dimension from RTMZ data multiple times in order to augment training examples (~10x)
    normalize_intensity : ``bool``, default=True
        normalize the intensities in each 2D RTMZ array so that they are in the range 0->1 
    debug_flag : ``str``, optional
        specifies how to dispatch the message and/or plot, None to do nothing
    debug_cb : ``func``, optional
        callback function that takes the debugging message as an argument, can be None if
        debug_flag is not set to 'textcb'

    Returns
    -------
    ml_data : ``numpy.ndarray``
        array of binned data for ML with shape: (targets, 3, 24, 400)
    """
    binned_data = []
    i = 0
    n = len(preml_data['targets'])
    for v in preml_data['targets'].values():
        rt_bounds = [v['rt'] - 0.2, v['rt'] + 0.2]
        pre = _get_binned_data_from_rtmz_arrays(v['pre_data'], rt_bounds, [v['pre_mz'] - 1.5, v['pre_mz'] + 2.5])
        ald = _get_binned_data_from_rtmz_arrays(v['ald_data'], rt_bounds, [v['ald_mz'] - 1.5, v['ald_mz'] + 2.5])
        crg = _get_binned_data_from_rtmz_arrays(v['crg_data'], rt_bounds, [v['crg_mz'] - 1.5, v['crg_mz'] + 2.5])
        if normalize_intensity:
            pre /= np.max(pre)
            ald /= np.max(ald)
            crg /= np.max(crg)
        binned_data.append([pre, ald, crg])
        if rt_sampling_augment:
            for rtb in _get_aug_rt_bounds(v['rt']):
                pre = _get_binned_data_from_rtmz_arrays(v['pre_data'], rtb, [v['pre_mz'] - 1.5, v['pre_mz'] + 2.5])
                ald = _get_binned_data_from_rtmz_arrays(v['ald_data'], rtb, [v['ald_mz'] - 1.5, v['ald_mz'] + 2.5])
                crg = _get_binned_data_from_rtmz_arrays(v['crg_data'], rtb, [v['crg_mz'] - 1.5, v['crg_mz'] + 2.5])
                if normalize_intensity:
                    pre /= np.max(pre)
                    ald /= np.max(ald)
                    crg /= np.max(crg)
                binned_data.append([pre, ald, crg])
        i += 1
        msg = 'finished binning {:3d} of {:3d} targets'.format(i, n)
        _debug_handler(debug_flag, debug_cb, msg=msg)
    return np.array(binned_data)


class _OzIDDataset(Dataset):
    """ internal class used for loading data for ML, subclassed from ``torch.utils.data.Dataset`` """

    def __init__(self, true_data, false_data, transform=None):
        """
        inits an instance of _OzIDDataset using True and False data

        Parameters
        ----------
        true_data : ``numpy.ndarray``
        false_data : ``numpy.ndarray``
            arrays of true, false binned data for ML with shapes (N, 3, 24, 400), where N is the number of 
            training examples in each set
        transform : ``?``, optional
            transform to apply to data when accessed (e.g. converting to FloatTensor)
        """
        self.true_data = true_data
        self.false_data = false_data
        self.transform = transform
        self.num_true = len(self.true_data)
        self.num_false = len(self.false_data)
        self.total_num = self.num_true + self.num_false
        self.classes = ["False", "True"]

    def __len__(self):
        return self.total_num

    def __getitem__(self, idx):
        if idx < self.num_true:
            mat = self.true_data[idx]
            label = 1
        else:
            mat = self.false_data[idx-self.num_true]
            label = 0
        # reshape for ToTensor
        # numpy.ndarray (rt_bins x mz_bins x 3) -> (3 x rt_bins x mz_bins)
        mat = mat.transpose(1, 2, 0)
        if self.transform:
            mat = self.transform(mat)
        # training example + corresponding label
        return mat, label


def _split_training_and_val_data(true_data, false_data, val_size, random_state):
    """ returns true_train, true_val, false_train, false_val """
    ssplit = ShuffleSplit(n_splits=1, test_size=val_size, random_state=random_state)
    for train_index, val_index in ssplit.split(true_data):
        true_train, true_val = true_data[train_index], true_data[val_index]
    for train_index, val_index in ssplit.split(false_data):
        false_train, false_val = false_data[train_index], false_data[val_index]
    return true_train, true_val, false_train, false_val


def get_dataloaders_for_ml(true_data, false_data, val_size=0.2, batch_size=128, shuffle=True, random_state=420):
    """
    Splits true/false data into training/validation sets, returns ``torch.utils.data.DataLoader`` instances for each
    along with corresponding dataset sizes. Includes the ``torchvision.transforms.ToTensor`` transform in datasets

    Parameters
    ----------
    true_data : ``numpy.ndarray``
    false_data : ``numpy.ndarray``
        arrays of true, false binned data for ML with shapes (N, 3, 24, 400), where N is the number of 
        training examples in each set
    val_size : ``float``, default=0.2
        proportion of dataset to split into validation set
    batch_size : ``int``, default=128
        ``batch_size`` parameter for dataloaders, given the proportion of True/False samples in the complete trianing
        data (~7.5%), a batch size of 128 should contain around 10 True examples on average
    shuffle : ``bool``, default=True
        ``shuffle`` parameter for dataloaders
    random_state : ``int``, default=420
        pRNG seed for deterministic splitting results

    Returns
    -------
    dataloaders : ``dict(str:torch.utils.data.DataLoader)``
        'train' and 'validate' Dataloaders (``torch.utils.data.DataLoader``)
    dataset_sizes : ``dict(str:int)``
        'train' and 'validate' dataset sizes
    """
    transform = transforms.Compose([transforms.ToTensor()])
    true_train, true_val, false_train, false_val = _split_training_and_val_data(true_data, false_data, 
                                                                                val_size, random_state)
    train_dset = _OzIDDataset(true_data=true_train, false_data=false_train, transform=transform)
    val_dset = _OzIDDataset(true_data=true_val, false_data=false_val, transform=transform)

    train_dataloader = DataLoader(train_dset, batch_size=batch_size, shuffle=shuffle)
    val_dataloader = DataLoader(val_dset, batch_size=batch_size, shuffle=shuffle)
    dataloaders = {'train': train_dataloader, 'val': val_dataloader}
    dataset_sizes = {'train': len(train_dset), 'val': len(val_dset)}
    return dataloaders, dataset_sizes

