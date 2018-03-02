#! python
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
import scipy.ndimage.morphology as morph
from scipy import ndimage
import numpy as np
import os
import sys

__author__ = 'Paul Hancock'
__date__ = ''


def filter_cat(catalogue, image, outcat):
    im = fits.open(image)
    cat = Table.read(catalogue).to_pandas()
    wcs = WCS(im[0].header, naxis=2)
    # Convert ra/dec to pixel values (in x,y)
    pix = np.int32(wcs.all_world2pix(cat[['ra', 'dec']].values, 1)).T
    struct1 = ndimage.generate_binary_structure(2, 1)
    dl = np.bitwise_not(morph.binary_dilation(np.bitwise_not(np.isfinite(im[0].data)),
                                              iterations=3,
                                              structure=struct1))
    # image mask is where we *havent* masked the image (note [y,x])
    mask = np.where(dl[pix[1], pix[0]])[0]
    Table.from_pandas(cat.iloc[mask]).write(outcat,  overwrite=True)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print "usage ${0} incat.fits image.fits outcat.fits".format(__file__)
        sys.exit(1)
    catalogue, image, outcat = sys.argv[-3:]
    if os.path.exists(outcat):
        os.remove(outcat)
    filter_cat(catalogue, image, outcat)
