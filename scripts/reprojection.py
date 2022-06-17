#! /usr/bin/env python3

import glob
import os
import numpy as np
import argparse
import sys

from astropy.io import fits
from reproject import reproject_interp

def get_epoch_catalogues(epochs_file):
    """
    Read a file which contains a list of the catalogues to be read

    parameters
    ----------
    epochs_files : str
        A file which has a list of catalogues, one per line.

    returns
    -------
    files : list
        A list of filenames
    """
    files = list(map(str.strip, open(epochs_file).readlines()))
    return files


def fits_file_reprojection(epoch_files, mean_file, output_dir):

    # Get list of Epoch's and then mean image
    fits_files = epoch_files
    fits_files.append(mean_file[0])
  
    # Iterate through FITS files
    for f_file in fits_files:
        hdu = fits.open(f_file)[0]
        # Edit header
        new_header = hdu.header.copy()
        new_header['CTYPE1'] = 'RA---CAR'
        new_header['CTYPE2'] = 'DEC--CAR'

        if 'mean' not in f_file:
            new_header['CRPIX2'] -= new_header['CRVAL2'] / new_header['CDELT2']
            new_header['CRPIX1'] -= new_header['CRVAL1'] / new_header['CDELT1']
        else:
            new_header['CRPIX2'] -= new_header['CRVAL2'] / new_header['CD2_2']
            new_header['CRPIX1'] -= new_header['CRVAL1'] / new_header['CD1_1']

        new_header['CRVAL1'] = 0.0
        new_header['CRVAL2'] = 0.0
            

        new_image, footprint = reproject_interp(hdu, new_header) 

        # Change dtype to float32 for visualisation
        new_image = new_image.astype(np.float32)

        # Squeeze to remove extra dimensions
        new_image = np.squeeze(new_image)

        # Write out file 
        f_fileout_name = f_file.replace('.fits', '_reprojected.fits').split('/')[-1]
      
        fits.writeto(f'{output_dir}/{f_fileout_name}', new_image, new_header, overwrite=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Reproject epoch files for visualisation")
    group1.add_argument("--epochs", dest='epochs', type=str, default=None,
                        help="A file containing the list of epoch catalogues")
    group1.add_argument("--mean", dest='mean', type=str, default=None,
                        help="The mean image file")
    group1.add_argument("--reproj_dir", dest='reproject_img_dir', type=str, default=None,
                        help="The reprojected image output directory")                         
    args = parser.parse_args()

    if None in (args.epochs, args.mean):
        parser.print_help()
        sys.exit()

    files = get_epoch_catalogues(args.epochs)
    mean = get_epoch_catalogues(args.mean)

    fits_file_reprojection(files, mean, args.reproject_img_dir)