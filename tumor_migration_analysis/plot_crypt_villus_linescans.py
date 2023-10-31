#!/opt/local/bin/python

"""

Plots linescans from average PIV data from uGut tissues.
Linescan data can be extracted easily using FIJI (with a wide linescan from crypt region to villus region and plot profile)

"""

import os
from copy import deepcopy

import numpy as np
from skimage.io._plugins import tifffile_plugin as tifffile
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import utility_functions as uf

def plot_linescan(ls_path,px_size=1,shift=0):

    ls_data = uf.get_dict_list(uf.read_file(ls_path,delim=","))
    print(ls_data)
    x_px = np.array([_['Distance_(pixels)'] for _ in ls_data[shift:]], dtype='float')
    x_um = x_px * px_size
    y = np.array([_['Gray_Value'] for _ in ls_data[shift:]], dtype='float')
    y_smooth = uf.smooth_moving_window(y, window_len=11, include_edges="On")

    fig,ax = plt.subplots(figsize=(10,5))
    ax.plot(x_um,y,'ro',ms=8,alpha=0.5)
    ax.plot(x_um,y_smooth,'r-',lw=2)
    ax.set_xlabel("Distance along Crypt-Villus Axis ($\mu$m)")
    ax.set_ylabel("Speed along Crypt-Villus Axis\n($\mu$m/min)")
    # plt.show()
    plt.tight_layout()
    plt.savefig(ls_path[:-4] + ".pdf")

def main():
    """Sets up the analysis for extracting PIV vectors.
    You should update the image path and pixel size here.
    You should not have to change anything in the rest of the script.

    """

    # #sets some initial parameters
    data_dir = "./migration"
    filename_list = ["JC8_WT_day3_1_crop _C_trimmed_edges_manual_mask_speed_along_c-v_axis_interp_linescan.csv",
                     "JC8_WT_day3_1_crop _C_trimmed_edges_manual_mask_speed_along_c-v_axis_linescan.csv"]
    px_size = 0.8817564589 * 25 / 1000 #mm/px (divided by half of interrogation window)

    for filename in filename_list:

        print(filename)
        ls_path = os.path.join(data_dir,filename)
        plot_linescan(ls_path, px_size=px_size,shift=1)

if __name__=="__main__":
    main()