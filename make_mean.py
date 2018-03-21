#!/usr/bin/env python
from __future__ import print_function

from astropy.io import fits
import numpy as np
import sys
import os


def mean(cube, out):
    """
    Flatten an image cube into a median image.
    Write an output file

    Prameters
    ---------
    cube : str
        Name of the cube image

    out : str
        Output file name
    """
    print("reading {0}".format(cube))
    hdu = fits.open(cube)
    data = hdu[0].data

    print("mean calc")
    median = np.mean(data, axis=0)
    hdu[0].data = median

    hdu.writeto(out)
    print("wrote {0}".format(out))
    return


if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print("make_mean.py cube outfile")
        sys.exit(1)

    infile = sys.argv[1]
    outfile = sys.argv[2]

    mean(infile, out=outfile)
