""" 
lipidoz_gui/_config.py

Dylan Ross (dylan.ross@pnnl.gov)

    internal module for configuration parameters
"""


import os


# mza version
MZA_VERSION = 'new'


# default values for data extraction parameters
DEFAULT_RT_TOL = 0.1
DEFAULT_RT_EXT_WIN = 1.5
DEFAULT_MZ_TOL = 0.01


# results window plot height/width (in inches)
PLOT_HEIGHT = 175
XIC_PLOT_WIDTH = int(1.20 * PLOT_HEIGHT)
MS1_PLOT_WIDTH = int(3.00 * PLOT_HEIGHT)


# color for saturation correction symbol on XICs
SAT_CORR_COLOR = (255, 0, 0)


# set the application icon path
ICO_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'lipidoz.ico')


# detect whether we are running on Windows
PLATFORM_IS_WINDOWS = os.name == 'nt'


# print complete debugging info to processing window
PROC_DEBUG = False
