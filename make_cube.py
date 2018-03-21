#!/usr/bin/env python
from __future__ import print_function

from astropy.io import fits
import numpy as np
import sys
import os

def stack(files, out):
    """
    Combine a list of fits files into a single cube and save
    the output to out.

    Parameters
    ----------
    files : list
        List of files
    
    out : str
        Filename to save
    """

    ref = fits.open(files[0])
    data = np.empty((len(files), ref[0].data.shape[0], ref[0].data.shape[1]),
                    dtype=np.float32)

    for i, f in enumerate(files):
        print('add {0}'.format(f))
        hdu = fits.open(f)
        data[i, :, :] = hdu[0].data
    
    ref[0].data = data
    ref.writeto(out, overwrite=True)
    print("wrote {0}".format(out))

if __name__ == "__main__":

    if len(sys.argv) <= 2:
        print("make_cube.py outfile.fits file1.fits file2.fits ...")
        sys.exit(1)

    outname = sys.argv[1]
    files = sys.argv[2:]

    if len(files) < 2:
        print("not enough files, need at least 2 to make a cube")
        print("given {0}".format(files))
        sys.exit(1)
    stack(files=files, out=outname)
