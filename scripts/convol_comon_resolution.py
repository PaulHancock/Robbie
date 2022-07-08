#! /usr/bin/env python3

__author__ = 'Paul Hancock'
__date__ = '2022/07/08'

import argparse
import os
import sys

import AegeanTools.wcs_helpers as wcsh
import numpy as np
from astropy.convolution import Gaussian2DKernel
from astropy.convolution.kernels import Gaussian2DKernel
from astropy.io import fits
from radio_beam import Beam, Beams
from scipy import ndimage, optimize, signal

SIGMA_TO_FWHM = np.sqrt(8*np.log(2))


def get_common_beam(fnames):
    """
    Read headers from a list of fits files and determine the smallest
    common beam to which they could be convolved

    parameters
    ----------
    fnames : list([str, ...])
        List of filenames

    return
    ------
    cb : `radio_beam.Beam`
        The smallest common beam
    """
    headers = [fits.getheader(f) for f in fnames]
    beams = [Beam.from_fits_header(h) for h in headers]
    beams=Beams(beams=beams)
    cb = beams.common_beam()
    return cb


def beam_to_kernel(beam,
                   ra, dec,
                   wcshelper):
    """
    Convert a beam object into a Gaussian kernel.
    Since the sky->pix mapping may not be affine the sky location is required for this conversion.

    parameters
    ----------
    beam : `radio_beam.Beam`
        The beam object
    
    ra, dec : float, float
        The sky location
    
    wcshelper: `AegeanTools.wcs_helpers.wcshelper`
        An augmented WCS object.
    
    returns
    -------
    kern : `astropy.convolution.Gaussian2DKernel`
        The kernel
    """
    _, _, a, b, pa = wcshelper.sky2pix_ellipse((ra,dec), beam.major.value, beam.minor.value, beam.pa.value)
    kern = Gaussian2DKernel(a/SIGMA_TO_FWHM, y_stddev=b/SIGMA_TO_FWHM, theta=pa)
    return kern


def find_gaussian_kernel(image_kern, target_kern):
    """
    Find the kernel B which when convolved with image_kern will produce target_kern
    This is solving A*B = C for B, where A and C are known and * is the 2d convolution.

    parameters
    ----------
    image_kern, target_kern : `astropy.convolution.Gaussian2DKernel`
        The kernels A and C
    
    return
    ------
    B : `astropy.convolution.Gaussian2DKernel`
    """
    def diff(kern):
        return np.sum(np.abs(signal.convolve2d(kern.array, image_kern.array, mode='same')- target_kern.array))


    def func(x):
        a,b,pa = x
        B = Gaussian2DKernel(a/SIGMA_TO_FWHM, y_stddev=b/SIGMA_TO_FWHM, theta=pa, 
                            x_size=image_kern.shape[0], y_size=image_kern.shape[1])
        return diff(B)
    
    res = optimize.minimize(func, x0=[1,1,0])

    B_fit = Gaussian2DKernel(res.x[0]/SIGMA_TO_FWHM, y_stddev=res.x[1]/SIGMA_TO_FWHM, theta=res.x[2], 
                         x_size=image_kern.shape[0], y_size=image_kern.shape[1])
    return B_fit


def convolve_to(psf, infile, outfile):
    """
    Convolved an image to the target psf.
    Save the resulting image.

    parameters
    ----------
    psf : `radio_beam.Beam`
        The target psf (must be same/larger than that of the input file)

    infile : str
        Path to read input image

    outfile : str
        Path to save output
    """
    hdu = fits.open(infile)

    # determine the current psf
    wcshelper = wcsh.WCSHelper.from_header(hdu[0].header)
    beam = Beam.from_fits_header(hdu[0].header)
    im_kern = beam_to_kernel(beam, hdu[0].header['CRVAL1'], hdu[0].header['CRVAL2'], wcshelper)
    t_kern = beam_to_kernel(psf, hdu[0].header['CRVAL1'], hdu[0].header['CRVAL2'], wcshelper)

    # find the kernel that will get us from im_kern -> t_kern
    fit_kern = find_gaussian_kernel(im_kern, t_kern)

    # convolve and update the image data
    data = signal.convolve2d(hdu[0].data, fit_kern, mode='same')
    hdu[0].data = data

    # update the header info
    header = hdu[0].header
    header['BMAJ'] = psf.major.value
    header['BMIN'] = psf.minor.value
    header['BPA'] = psf.pa.value

    # save the file
    hdu.writeto(outfile, overwrite=True)
    print("Wrote {0}".format(outfile))
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Convolve images")
    group1.add_argument("--in", dest='files', type=str, default=None, nargs='+',
                        help="List of input images")
    group1.add_argument("--suffix", dest='suffix', type=str, default="convolved",
                        help="Output suffix for filenames.")
    results = parser.parse_args()

    if results.files is None or len(results.files) <2:
        parser.print_help()
        sys.exit(0)

    psf = get_common_beam(results.files)
    for f in results.files:
        base, ext = os.path.splitext(f)
        outfile = "{0}_{1}{2}".format(base, results.suffix, ext)
        convolve_to(psf, f, outfile)

