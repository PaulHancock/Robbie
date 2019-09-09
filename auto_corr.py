#! /usr/bin/env python
from __future__ import print_function
__author__=['Tim White', 'Paul Hancock']
__date__ = '2019/09/06'

import argparse
from astropy.io import fits
import sqlite3
import numpy as np
import sys
import os

NPIX = 500


def autocorr(x):
    corr = np.correlate(x, x, mode='full')
    half = corr[corr.size/2:]
    normed = half/half[0]
    return normed


def get_pixels(n, cube):
    # TODO learn how to avoid nan pixels
    x = np.random.randint(cube.shape[1],size=n)
    y = np.random.randint(cube.shape[2],size=n)
    pix = cube[:,x,y]
    return pix


def make_acorr_stats(n,cube):
    pix = get_pixels(n,cube)
    acorr = np.array( [ autocorr(pix[:,i]) for i in range(n)])
    mean = np.nanmean(acorr,axis=0)
    std = np.nanstd(acorr,axis=0)
    return mean, std


def get_cube_endof(cube):
    nsamples = cube.shape[0]
    mean, std = make_acorr_stats(NPIX,cube)
    # detect the first element that has correlation consistent with zero
    # print(mean)
    # print(std)
    # print(mean-std)
    fzero = np.min(np.where(mean-std < 0))
    # print(fzero)
    ndof = nsamples - 1 - fzero
    return ndof


def get_db_endof(db):
    """

    parameters
    ----------
    db : str
       Database filename.

    return
    ------
    endof : int
       The effective number of degrees of freedom
    """
    conn = sqlite3.connect(db)
    cur = conn.cursor()
    cur.execute("SELECT count(*) FROM epochs")
    epochs = cur.fetchone()[0]
    cur.execute("SELECT uuid FROM stats ORDER BY RANDOM() LIMIT ?", (NPIX,))
    uuids = [a[0] for a in cur.fetchall()]

    fluxes = np.zeros(shape=(epochs, NPIX))

    for i, uuid in enumerate(uuids):
        cur.execute("""
        SELECT peak_flux - mean_peak_flux
        FROM sources JOIN stats ON sources.uuid = stats.uuid WHERE stats.uuid = ? """, (uuid,))
        fluxes[:, i] = zip(*cur.fetchall())[0]

    acorr = np.array([autocorr(fluxes[:,i]) for i in range(NPIX)])
    mean = np.nanmean(acorr, axis=0)
    std = np.nanstd(acorr, axis=0)
    fzero = np.min(np.where(mean-std < 0))
    ndof = epochs - 1 - fzero
    return ndof



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Combine images into a cube")
    group1.add_argument("--dbname", dest='db', type=str, default=None,
                        help="Database name. [optional]")
    group1.add_argument('--cuebname', dest='cube', type=str, default=None,
                        help='Image cube name. [optional]')
    results = parser.parse_args()

    if (results.db is None) and (results.cube is None):
        print("ERROR: Either dbname or cubename need to be specified.")
        parser.print_help()
        sys.exit(1)

    if results.cube:
        print("Reading cube from {0}".format(results.cube))
        cube = fits.open(results.cube)[0].data
        if len(cube.shape) != 3:
            print("{0} needs to have 3 axes, but it has {1}".format(results.cube, len(cube.shape)))
            sys.exit(1)
        print("Effective degrees of freedom: {0}".format(get_cube_endof(cube)))
    elif results.db:
        print("Reading data from {0}".format(results.db))
        ndof = get_db_endof(results.db)
        print("Effective degrees of freedom: {0}".format(ndof))
