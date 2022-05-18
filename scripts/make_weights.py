#! /usr/bin/env python

from astropy.io import fits
import os
import argparse

__author__ = 'Paul Hancock'
__date__ = '2022-05-18'

def make_weight(fname, val, outfile=None):
    hdu = fits.open(fname)
    hdu[0].data *= 0
    hdu[0].data += val
    
    if outfile is None:
        outfile = os.path.basename(fname)
        outfile = outfile.replace('.fits','.weight.fits')

    hdu.writeto(outfile, overwrite=True)
    print(f"Wrote {outfile}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("file", type=str,
                        help="Filename (fits image)")
    parser.add_argument("-w", dest="weight",type=float, default=1,
                        help="Weight")
    parser.add_argument('-o','--out', dest="outfile", type=str, default=None,
                        help="outfile name")
    results = parser.parse_args()
    make_weight(results.file, results.weight, results.outfile)