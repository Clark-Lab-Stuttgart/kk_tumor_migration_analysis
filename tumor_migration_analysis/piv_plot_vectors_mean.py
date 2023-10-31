#!/opt/local/bin/python

"""

Plots PIV vector data (just vectors, no fame)

"""

import os
from copy import deepcopy

import numpy as np
from skimage.io._plugins import tifffile_plugin as tifffile
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import utility_functions as uf

def plot_vectors(stk_path,px_size=1,scale_factor=0.004,scale_length=0.1,vector_width=1.5):
    """Plots PIV vectors overlaid onto the orginal image stack
    The images are automatically written to a new directory.

    Parameters
    ----------
    stk_path : string
        the path to the image stack to be analyzed
    px_size : float
        the pixel size in um/px
    scale_factor : float
        a scaling to determine the vector size
    scale_length : int
        the length of the scale vector (in scaled units, usually um/min)
    vector_width : float
        the width of the PIV vectors for the quiver plots

    """

    # opens the image stack to get the aspect ratio
    stk = tifffile.imread(stk_path)
    width = stk[0].shape[1] * px_size
    height = stk[0].shape[0] * px_size
    ar = height / width

    # finds the PIV vector data
    data_dir = os.path.splitext(stk_path)[0] + "_piv_data"
    if not os.path.isdir(data_dir):
        raise FileNotFoundError("No PIV vector data found. Please run extraction script first.")

    # get unique basename list (from x coordinate data)
    basename_list = [_[:-6] for _ in os.listdir(data_dir) if '_x.dat' in _]
    basename_list = uf.natural_sort(basename_list)

    for tag in ["", "_interp"]: #plots both the raw and interpolated data

        print("Tag = ", tag)

        #makes new directory for plotting
        plot_dir = os.path.join(data_dir,"vectors_w%sum-min_scale"%(str(scale_length).replace(".","p")) + tag)
        uf.make_dir(plot_dir)

        x = np.array(uf.read_file(os.path.join(data_dir, basename_list[1] + "_x.dat")), dtype=float)
        y = np.array(uf.read_file(os.path.join(data_dir, basename_list[1] + "_y.dat")), dtype=float)

        U_stk = np.zeros(shape=(len(basename_list) - 1, x.shape[0], x.shape[1]))
        V_stk = deepcopy(U_stk)

        print(x.shape)
        print(U_stk.shape)

        #goes through each time frame
        for i, basename in enumerate(basename_list[1:]):

            #plots each frame
            print('Grabbing frame:', basename)
            U = np.array(uf.read_file(os.path.join(data_dir, basename + "_u%s.dat" % tag)),dtype=float)
            V = np.array(uf.read_file(os.path.join(data_dir, basename + "_v%s.dat" % tag)),dtype=float)

            U_stk[i] = U
            V_stk[i] = V


        print(U_stk)
        print(U_stk.shape)

        U_mean = np.nanmean(U_stk, axis=0)
        V_mean = np.nanmean(V_stk, axis=0)
        print(U_mean.shape)


        print(basename)
        basename = "_".join(basename.split("_")[:-1])
        print(basename)

        #saves mean data
        uf.save_data_array(U_mean,plot_dir + '/%s_mean_vectors_u%s.dat'%(basename,tag))
        uf.save_data_array(V_mean,plot_dir + '/%s_mean_vectors_v%s.dat'%(basename,tag))

        #sets up the plot and plots the image data underneath
        rcParams['axes.linewidth'] = 0
        fig, ax = plt.subplots(figsize=(6,6*ar))
        # ax.patch.set_facecolor('black')
        ax.imshow(stk[0],cmap='Greys_r',extent=(0,width,height,0))

        #plots vectors with color code according to angle
        phis = np.arctan2(V_mean, U_mean) * 180. / np.pi
        plt.quiver(x, y, U_mean, V_mean, phis.ravel(), cmap='rainbow', clim=(-180., 180.),
                   units='xy',scale=scale_factor,scale_units='x',width=vector_width)

        #makes an arrow for scale
        cm = LinearSegmentedColormap.from_list('cm', [(1,1,1),(1,1,1)])
        plt.quiver([width - width * 0.04], [height - height * 0.02], [scale_length], [0], [0], cmap=cm,
                   units='xy', scale=scale_factor, scale_units='x',width=vector_width)

        #finishes plot
        ax.set_xlim(0,width)
        ax.set_ylim(0,height)
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        fig.subplots_adjust(bottom=0,left=0,top=1,right=1)

        #saves plot as png
        plt.savefig(plot_dir + '/%s_mean_vectors%s.png'%(basename,tag), dpi=600)
        plt.close()

        ##now make a plot with the migration speed along the crypt-villus axis
        cv_angle = 26.867 ##in degrees ##this should be in a metadata file somewhere, fine to hard-code for now
        angular_diff = np.deg2rad(phis - cv_angle)
        vector_mags = np.linalg.norm((U_mean,V_mean), axis=0)
        v_projected = vector_mags * np.cos(angular_diff) * 60 #um/hr instead of um/min

        #makes plot
        fig, ax = plt.subplots(figsize=(8,8*ar))
        im = plt.imshow(v_projected,cmap='bwr',vmin=-4,vmax=4)
        cb = fig.colorbar(im, fraction=0.03, pad=0.04)
        cb.set_label('Speed along C-V Axis ($\mu$m/hr)', rotation=90)
        ax.set_xticks([])
        ax.set_yticks([])
        # plt.show()
        plt.tight_layout()
        plt.savefig(plot_dir + '/%s_speed_along_c-v_axis%s.png'%(basename,tag), dpi=600)
        plt.close()

        #saves speed data
        uf.save_data_array(v_projected,os.path.join(plot_dir, basename + "_speed_along_c-v_axis%s.dat" % tag))

        #make a 1D average along C-V axis
        #easiest to do this in FIJI and then run the plotting script ()


def main():
    """Sets up the analysis for extracting PIV vectors.
    You should update the image path and pixel size here.
    You should not have to change anything in the rest of the script.

    """

    # #sets some initial parameters
    stk_path = './JC8_WT_day3_1_crop _C_trimmed_edges_manual_mask_bf_overlay.tif'
    px_size = 0.8817564589 #um/px
    scale_factor = 0.002 #for scaling the PIV vectors on the plots
    scale_length = 0.1 #sets the length of the scale vector on the plots (in um/min)
    # #sets some initial parameters
    # stk_path = './sample_data/tumor_nuclei_small/tumor_nuclei_small.tif'
    # px_size = 0.91 #um/px
    # scale_factor = 0.004 #for scaling the PIV vectors on the plots
    # scale_length = 0.1 #sets the length of the scale vector on the plots (in um/min)

    plot_vectors(stk_path,px_size=px_size,scale_factor=scale_factor,scale_length=scale_length,vector_width=3)

if __name__=="__main__":
    main()