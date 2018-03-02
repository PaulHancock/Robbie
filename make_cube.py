#! /usr/bin/env python
from __future__ import print_function

from glob import glob
from astropy.io import fits
import numpy as np
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
    ref.writeto(out)
    print("wrote {0}".format(out))

if __name__ == "__main__":
    for f in [154, 185]:
        fname = 'K2_{0}MHz.dat'.format(f)
        files = [a.split()[0] for a in open(fname).readlines()]
        # use the warped files to make the cube!
        files = [a[:-5]+'_warped.fits' for a in files]
        outname='cube_{0}MHz.fits'.format(f)
        if not os.path.exists(outname):
            stack(files=files, out=outname) 
