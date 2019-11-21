#!/usr/bin/env python
from AegeanTools.regions import Region
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
import argparse
import scipy.ndimage.morphology as morph
from scipy import ndimage
import numpy as np
import os
import sys

__author__ = 'Paul Hancock'
__date__ = '2018-08-29'


def filter_cat(catalogue, image, outcat, region=None):
    im = fits.open(image)
    cat = Table.read(catalogue).to_pandas()
    wcs = WCS(im[0].header, naxis=2)
    # Convert ra/dec to pixel values (in x,y)
    pix = np.int32(wcs.all_world2pix(cat[['ra', 'dec']].values, 1)).T
    
    # constrain pixel values to be within the image
    pix[0] = np.clip(pix[0], a_min=0, a_max=im[0].data.shape[0]-1)
    pix[1] = np.clip(pix[1], a_min=0, a_max=im[0].data.shape[1]-1)

    struct1 = ndimage.generate_binary_structure(2, 1)
    dl = np.bitwise_not(morph.binary_dilation(np.bitwise_not(np.isfinite(im[0].data)),
                                              iterations=3,
                                              structure=struct1))
    # image mask is where we *havent* masked the image (note [y,x])
    mask = np.where(dl[pix[1], pix[0]])[0]
    tab = Table.from_pandas(cat.iloc[mask])
    if region is not None:
        # exclude regions of sky that are outside of the mask.
        reg = Region.load(region)
        ra = tab['ra']
        dec = tab['dec']
        mask = reg.sky_within(ra, dec, degin=True)
        tab = tab[mask]
    # don't write empty files
    if len(tab) > 0:
        tab.write(outcat,  overwrite=True)
        print("Wrote {0}".format(outcat))
    else:
        print("Empty table. No output written.")
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Filter a transients catalogue")
    group1.add_argument("--incat", dest='incat', type=str, default=None,
                        help="The input catalogue.")
    group1.add_argument("--image", dest='image', type=str, default=None,
                        help='The input image')
    group1.add_argument('--region',dest='region', type=str, default=None,
                         help='A region file to provide additional filtering')
    group1.add_argument("--outcat", dest='outcat', type=str, default=None,
                        help="The output catalogue")

    results = parser.parse_args()

    if None in [results.incat, results.image, results.outcat]:
        parser.print_help()
        sys.exit(1)

    filter_cat(catalogue=results.incat,
               image=results.image,
               outcat=results.outcat,
               region=results.region)
