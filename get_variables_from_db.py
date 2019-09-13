#! /usr/bin/env python

from __future__ import print_function, division

import sqlite3
import astropy
import astropy.table
import numpy as np
import argparse
import os


__author__  = ["Paul Hancock"]
__date__ = '2019/09/13'


def get_variable_stats(cur):
    """
    Extract the stats table from the given database and return the corresponding table.
    Returns all the information from the stats table, as well as the average ra/dec from the sources table.

    parameters
    =========
    cur : sqlite.connection.cursor
        A cursor into the database of interest

    returns
    =======
    tab : astropy.table.Table
        The same data as an astropy table
    """
    # grab the data, doing the average ra/dec from the sources table
    cur.execute("""
    SELECT avg(ra) as ra, avg(dec) as dec, stats.uuid,
    stats.mean_peak_flux, stats.std_peak_flux, stats.chisq_peak_flux, stats.m, stats.md, stats.pval_peak_flux
    FROM  sources JOIN stats ON sources.uuid = stats.uuid
    GROUP by sources.uuid
    """)

    # convert the table data into numpy arrays and then astropy Table columns
    rows = cur.fetchall()
    ra, dec, uuid, mean_peak_flux, std_peak_flux, chisq_peak_flux, m, md, pval_peak_flux =  map(np.array, zip(*rows))
    ra = astropy.table.Column(data=ra, name='ra')
    dec = astropy.table.Column(data=dec, name='dec')
    uuid = astropy.table.Column(data=uuid, name='uuid')#, dtype=str) #need to be explicit with type here or the final fits file ends up with floats
    mean_peak_flux = astropy.table.Column(data=mean_peak_flux, name='mean_peak_flux')
    std_peak_flux = astropy.table.Column(data=std_peak_flux, name='std_peak_flux')
    chisq_peak_flux = astropy.table.Column(data=chisq_peak_flux, name='chisq_peak_flux')
    m = astropy.table.Column(data=m, name='m')
    md = astropy.table.Column(data=md, name='md')
    pval_peak_flux = astropy.table.Column(data=pval_peak_flux, name='pval_peak_flux')

    # smoosh the columns into a table
    tab = astropy.table.Table()
    tab.add_columns([uuid,ra,dec,mean_peak_flux,std_peak_flux,chisq_peak_flux, m,md,pval_peak_flux])
    return tab

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Export variability stats from db to a fits file")
    group1.add_argument("--name", dest='name', type=str, default=None,
                        help="Database filename.")
    group1.add_argument("--out", dest='out', type=str, default=None,
                        help="Output fits file")
    results = parser.parse_args()

    if None in (results.name, results.out):
        parser.print_help()
        sys.exit(1)

    conn = sqlite3.connect(results.name)
    c = conn.cursor()
    tab = get_variable_stats(c)

    if os.path.exists(results.out):
        os.remove(results.out)
    tab.write(results.out)
