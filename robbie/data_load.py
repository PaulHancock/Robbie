"""
Loads all the data required by robbi from the data directory.
"""

import os

datadir = os.path.join(os.path.dirname(__file__))

# Hard code the path of the reference catalogue will be
REF_CAT = os.path.join(datadir, 'GLEAM_ref_cat.fits')