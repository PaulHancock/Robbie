#! /usr/bin/env python3

__author__ = 'Paul Hancock'
__date__ = '2022/05/30'

import argparse
from astropy.io import fits
import astropy.units as u
from radio_beam import Beam, Beams
import numpy as np
import AegeanTools.wcs_helpers as wcsh
from scipy import ndimage, signal, optimize
from astropy.convolution.kernels import Gaussian2DKernel

SIGMA_TO_FWHM = np.sqrt(8*np.log(2))

def get_common_beam(fnames):
    """
    """
    return cb


def find_gaussian_kernel(image_kern, target_kern):
    """
    """
    def diff(kern):
        return np.sum(np.abs(signal.convolve2d(kern.array, image_kern.array, mode='same')- target_kern.array))


    def func(x):
        a,b,pa = x
        B = Gaussian2DKernel(a/SIGMA_TO_FWHM, y_stddev=b/SIGMA_TO_FWHM, theta=pa, 
                            x_size=image_kern.data.shape[0], y_size=image_kern.data.shape[1])
        return diff(B)
    
    res = optimize.minimize(func, x0=[1,1,0])

    B_fit = Gaussian2DKernel(res.x[0]/SIGMA_TO_FWHM, y_stddev=res.x[1]/SIGMA_TO_FWHM, theta=res.x[2], 
                         x_size=image_kern.data.shape[0], y_size=image_kern.data.shape[1])
    return B_fit



def convolve_to(kern, image):
    return image



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Convolve images")
    group1.add_argument("--in", dest='files', type=str, default=None, nargs='+',
                        help="List of input images")
    group1.add_argument("--suffix", dest='suffix', type=str, default="convolved",
                        help="Output suffix for filenames.")
    results = parser.parse_args()

