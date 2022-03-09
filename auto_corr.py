#! /usr/bin/env python
from __future__ import print_function
__author__=['Tim White', 'Paul Hancock']
__date__ = '2019/11/21'

import argparse
from astropy.io import fits
from astropy.table import Table
import sqlite3
import numpy as np
import sys
import os

NPIX = 500


def autocorr(x):
    corr = np.correlate(x, x, mode='full')
    half = corr[corr.size//2:]
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
    cur.execute("SELECT DISTINCT uuid FROM sources ORDER BY RANDOM() LIMIT ?", (NPIX,))
    uuids = [a[0] for a in cur.fetchall()]

    fluxes = np.zeros(shape=(epochs, NPIX))

    for i, uuid in enumerate(uuids):
        cur.execute("""
        SELECT s1.peak_flux - s2.mean_peak_flux
        FROM sources as s1 JOIN (SELECT uuid, AVG(peak_flux) as mean_peak_flux FROM sources WHERE uuid=?) as s2 ON s1.uuid = s2.uuid WHERE s1.uuid = ? """, (uuid,uuid))
        fluxes[:, i] = list(zip(*cur.fetchall()))[0]

    acorr = np.array([autocorr(fluxes[:,i]) for i in range(NPIX)])
    mean = np.nanmean(acorr, axis=0)
    std = np.nanstd(acorr, axis=0)
    fzero = np.min(np.where(mean-std < 0))
    ndof = epochs - 1 - fzero
    return ndof


def get_table_endof(filename):
    """

    parameters
    ----------
    filename : str
       Table filename.

    return
    ------
    endof : int
       The effective number of degrees of freedom
    """
    tab = Table.read(filename)
    flux_cols = [a for a in tab.colnames if a.startswith('peak_flux')]

    # Choose N random rows without repetition
    nitems = min(NPIX, len(tab))
    idx = np.random.choice(range(len(tab)), nitems, replace=False)
    fluxes = tab[idx][flux_cols]

    # construct the autocorrelation and determine ndof
    # ensure that the light curves have zero mean
    acorr = np.array([autocorr(np.array(list(arr)) -  np.nanmean(list(arr)))
                      for arr in fluxes])
    mean = np.nanmean(acorr, axis=0)
    std = np.nanstd(acorr, axis=0)
    fzero = np.min(np.where(mean-std < 0))
    epochs = len(flux_cols)
    ndof = epochs - 1 - fzero
    return ndof


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Combine images into a cube")
    group1.add_argument("--dbname", dest='db', type=str, default=None,
                        help="Database name. [optional]")
    group1.add_argument('--cubename', dest='cube', type=str, default=None,
                        help='Image cube name. [optional]')
    group1.add_argument('--table', dest='table', type=str, default=None,
                        help='Flux table name. [optional]')
    results = parser.parse_args()

    ndof = 0 # default
    if results.cube:
        print("Reading cube from {0}".format(results.cube))
        cube = fits.open(results.cube)[0].data
        if len(cube.shape) != 3:
            print("{0} needs to have 3 axes, but it has {1}".format(results.cube, len(cube.shape)))
            sys.exit(1)
        ndof = get_cube_endof(cube)

    elif results.db:
        print("Reading data from {0}".format(results.db))
        ndof = get_db_endof(results.db)

    elif results.table:
        print("Reading data from {0}".format(results.table))
        ndof = get_table_endof(results.table)

    else:
        print("ERROR: One of dbname, cubename, or table need to be specified.")
        parser.print_help()
        sys.exit(1)

    print("Effective degrees of freedom: {0}".format(ndof))
