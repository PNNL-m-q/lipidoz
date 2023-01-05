"""
lipidoz/ml/view.py

Dylan Ross (dylan.ross@pnnl.gov)

    module with utilities for viewing ML data 
"""


import numpy as np
from matplotlib import pyplot as plt, rcParams, colors as mcolors


def plot_preml_example(pre_rtmz, ald_rtmz, crg_rtmz, rt, mzs, 
                        rt_tol=0.2, mz_range=(1.5, 2.5), figname=None, rgb=True):
    """
    Produces a plot of a pre-machine learning example which consists of raw RTMZ array data (arrays of retention time, 
    m/z, and intensities) for a precursor and aldehyde/criegee Oz fragments. Default rt_tol and mz_range correspond to
    the normal values used for extracting and binning ML data.

    Parameters
    ----------
    pre_rtmz : ``tuple(numpy.ndarray(float))``
    ald_rtmz : ``tuple(numpy.ndarray(float))``
    crg_rtmz : ``tuple(numpy.ndarray(float))``
        raw RTMZ data arrays for precursor and aldehyde/criegee Oz fragments
    rt : ``float``
        central retention time value to display (same for precursor and Oz fragments)
    mzs : ``tuple(float)``
        monoisotopic masses for precursor, aldehyde, and criegee
    rt_tol : ``float``, default=0.2
        tolerance of rt values (relative to rt) to display
    mz_range : ``tuple(float)``, default=(1.5, 2.5)
        lower, upper range of m/z values to display (relative to monoisotopic mass), same for precursor and Oz fragments
    figname : ``str``, optional
        if provided, save the image to the specified file
    rgb : ``bool``, default=True
        if True, use Reds, Greens, Blues color scales for each panel, otherwise use viridis for all
    """
    if rgb:
        colors = ['Reds', 'Greens', 'Blues']
    else:
        colors = ['plasma', 'plasma', 'plasma']
    rcParams['font.size'] = 6
    fig, axs = plt.subplots(nrows=3, figsize=(3, 2))
    for ax in axs:
        ax.set_facecolor('#CCCCCC')
    axs[0].scatter(pre_rtmz[1], pre_rtmz[0], c=pre_rtmz[2], cmap=colors[0], marker='s', s=0.2)
    axs[1].scatter(ald_rtmz[1], ald_rtmz[0], c=ald_rtmz[2], cmap=colors[1], marker='s', s=0.2)
    axs[2].scatter(crg_rtmz[1], crg_rtmz[0], c=crg_rtmz[2], cmap=colors[2], marker='s', s=0.2)
    # set the display bounds
    mz_lower, mz_upper = mz_range
    for ax, mz in zip(axs, mzs):
        ax.axvline(mz, ls='-', c='k', lw=0.5, zorder=-1)
        ax.set_ylim((rt - rt_tol, rt + rt_tol))
        ax.set_xlim((mz - mz_lower, mz + mz_upper))
        ax.axhline(rt, ls='-', c='k', lw=0.5, zorder=-1)
    fig.tight_layout()
    if figname is not None:
        plt.savefig(figname, dpi=350, bbox_inches='tight')
    else:
        plt.show()
    plt.close()


def plot_ml_example(pre_binned, ald_binned, crg_binned, figname=None, rgb=True):
    """
    Produces a plot of a machine learning example which consists of binned RTMZ data for a precursor and 
    aldehyde/criegee Oz fragments

    Parameters
    ----------
    pre_binned : ``numpy.ndarray(float)``
    ald_binned : ``numpy.nparray(float)``
    crg_binned : ``numpy.ndarray(float)``
        binned RTMZ data for precursor and aldehyde/criegee Oz fragments
    figname : ``str``, optional
        if provided, save the image to the specified file
    rgb : ``bool``, default=True
        if True, use Reds, Greens, Blues color scales for each panel, otherwise use viridis for all
    """
    if rgb:
        colors = [
            mcolors.LinearSegmentedColormap.from_list("reds_", ['k', 'r']),
            mcolors.LinearSegmentedColormap.from_list("greens_", ['k', 'g']),
            mcolors.LinearSegmentedColormap.from_list("blues_", ['k', 'b'])
        ]
    else:
        colors = ['plasma', 'plasma', 'plasma']
    fig, axs = plt.subplots(nrows=4, figsize=(3.333, 2.5))
    axs[0].imshow(pre_binned, cmap=colors[0], aspect='auto')
    axs[1].imshow(ald_binned, cmap=colors[1], aspect='auto')
    axs[2].imshow(crg_binned, cmap=colors[2], aspect='auto')
    axs[3].imshow(np.array([pre_binned, ald_binned, crg_binned]).transpose(1, 2, 0), aspect='auto')
    for ax in axs:
        ax.set_xticks([])
        ax.set_yticks([])
    plt.tight_layout()
    if figname is not None:
        plt.savefig(figname, dpi=350, bbox_inches='tight')
    else:
        plt.show()
    plt.close()

