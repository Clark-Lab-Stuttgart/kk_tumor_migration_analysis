#!/Users/clark/anaconda3/bin/python

"""


"""

import os

import tifffile

from piv_extract_vectors import extract_vectors
from piv_plot_vectors import plot_vectors
from piv_analyze_vectors import analyze_vectors

def main():

    data_dir = './sample_data/PIV_karen/sample_full'
    # data_dir = './sample_data/tumor_nuclei_small'

    for filename in os.listdir(data_dir):
        # if 'runway_zoom.tif' in filename:
        if ".tif" in filename and not "mask" in filename:

            mask_path = os.path.join(data_dir, os.path.splitext(filename)[0] + "_mask.tif")
            mask = tifffile.imread(mask_path)

            stk_path = os.path.join(data_dir,filename)

            time_int = 30  # min
            px_size = 0.822  # um/px
            window_len = 50  # interrogation window length in um
            scale_factor = 0.004  # for scaling the PIV vectors on the plots
            scale_length = 0.1  # sets the length of the scale vector on the plots (in um/min)

            # time_int = 30  # min
            # px_size = 0.91  # um/px
            # window_len = 20  # interrogation window length in um
            # scale_factor = 0.004  # for scaling the PIV vectors on the plots
            # scale_length = 0.1  # sets the length of the scale vector on the plots (in um/min)

            extract_vectors(stk_path,time_int=time_int,px_size=px_size,window_length=window_len,mask=mask)
            plot_vectors(stk_path,px_size=px_size,scale_factor=scale_factor,scale_length=scale_length)
            analyze_vectors(stk_path)

if __name__ == "__main__":
    main()
