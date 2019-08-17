#!/usr/bin/env python
from __future__ import print_function

from astropy.io import fits
import argparse
import numpy as np
import sys
import os


def mean_from_files(files, out):
    """
    Create a mean image from the input list of images.
    Write an output file

    parameters
    ----------
    files : [str, str, [str,...]]
        list of file names

    out : str
        Output file name
    """
    print("Reading {0}".format(files[0]))
    hdu = fits.open(files[0])
    data = hdu[0].data

    for f in files[1:]:
        print("Adding {0}".format(f))
        data += fits.getdata(f)
    data /= len(files)

    hdu[0].data = data
    hdu.writeto(out)
    print("Wrote {0}".format(out))
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Combine images into a mean image")
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
        print("not enough files, need at least 2 to make a mean image")
        print("given {0}".format(files))
        sys.exit(1)

    mean_from_files(files=files, out=results.outfile)
