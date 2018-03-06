#! /usr/bin/env python
from __future__ import print_function

from astropy.io import fits
import numpy as np
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
    cubes = ['cube_154MHz.fits', 'cube_185MHz.fits']
    for c in cubes:
        out = 'mean'+c[4:]
        if not os.path.exists(out):
            mean(c, out=out)
