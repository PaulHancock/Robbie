#!/usr/bin/env python
from __future__ import print_function

from astropy.io import fits
import numpy as np
import argparse
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
    data = np.empty((len(files), ref[0].data.shape[-2], ref[0].data.shape[-1]),
                    dtype=np.float32)

    for i, f in enumerate(files):
        print('add {0}'.format(f))
        hdu = fits.open(f)
        data[i, :, :] = hdu[0].data
    
    ref[0].data = data
    ref.writeto(out, overwrite=True)
    print("wrote {0}".format(out))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Combine images into a cube")
    group1.add_argument("--infile", dest='infile', type=str, default=None,
                        help="A list of fits images in a file. [optional]")
    group1.add_argument("--in", dest='files', type=str, default=None, nargs='+',
                        help="Explicit list of files to include.")
    group1.add_argument("--out", dest='outfile', type=str, default=None,
                        help="output filename")
    results = parser.parse_args()


    if (results.infile is None) and (results.files is None):
        parser.print_help()
        sys.exit()

    if results.infile is not None:
        files = [l.strip() for l in open(results.infile).readlines()]
    else:
        files = results.files

    if len(files) < 2:
        print("not enough files, need at least 2 to make a cube")
        print("given {0}".format(files))
        sys.exit(1)
    stack(files=files, out=results.outfile)
