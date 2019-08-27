#!/usr/bin/env python

from __future__ import print_function
__author__ = 'Paul Hancock'
__date__ = '2019/08/18'

import argparse
import sqlite3
import numpy as np
from scipy import stats
import sys



def make_table(cur):
    """
    Create a table to store the variability stats for each source
    Clear the table if it already exists

    parameters
    ----------
    cur : sqlite3.connection.cursor
        Cursor for database connection
    """
    cur.execute("""CREATE TABLE IF NOT EXISTS stats (
    uuid TEXT,
    mean_peak_flux NUMERIC,
    std_peak_flux NUMERIC,
    chisq_peak_flux NUMERIC,
    m NUMERIC,
    md NUMERIC,
    pval_peak_flux NUMERIC)""")
    # clear the table 
    cur.execute("DELETE FROM stats")
    return


def calc_stats(cur, ndof=None):
    """
    Compute the mean and std of each light curve

    parameters
    ----------
    cur : sqlite3.connection.cursor
        Cursor for the database connection
    """
    cur.execute("SELECT DISTINCT uuid FROM sources")
    sources = cur.fetchall()
    for s in sources:
        cur.execute("SELECT peak_flux FROM sources WHERE uuid=? ORDER BY epoch", s)
        fluxes = np.array([i[0] for i in cur.fetchall()])
        cur.execute("SELECT err_peak_flux FROM sources WHERE uuid=? ORDER BY epoch", s)
        err = np.array([i[0] for i in cur.fetchall()])
        # don't include fit errors in the stats calculation
        mask = np.where(err>0)
        npts = len(mask[0])
        
        if npts < 2:
            pval = 0.
            md = 0.
            mean = 0.
        else:
            # modulation index
            mean = np.mean(fluxes[mask])
            std = np.std(fluxes[mask])
            m = std/mean
            # chi squared
            chisq = np.sum((fluxes[mask] - mean)**2 / err[mask]**2)
            # pvalue
            if ndof is None:
                ndof = npts - 1
            pval = stats.chi2.sf(chisq, ndof)
            pval = max(pval, 1e-10)
            # debiased modulation index
            desc = np.sum((fluxes[mask] - mean)**2) - np.sum(err[mask]**2)
            #print(mean, desc, npts)
            md = 1./mean * np.sqrt(np.abs(desc)/npts)
            if desc < 0:
                md *= -1
        # add all to the stats table
        cur.execute("""INSERT INTO stats(uuid, mean_peak_flux, std_peak_flux, m, md, chisq_peak_flux, pval_peak_flux)
        VALUES (?,?,?,?,?,?,?)""",
                    (s[0], mean, std, m,  md, chisq, pval))
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Calculate variability stats")
    group1.add_argument("--name", dest='name', type=str, default=None,
                        help="Database filename.")
    group1.add_argument("--ndof", dest='ndof', type=float, default=None,
                        help="Effective number of degrees of freedom. Defualt: N=epochs-1")
    group1.add_argument("--correct", dest="correct", action='store_true', default=False,
                        help="Compute and remove the mean light curve (if there are more than 5 rows). Default: False.")
    results = parser.parse_args()

    if results.name is None:
        parser.print_help()
        sys.exit(1)

    conn = sqlite3.connect(results.name)
    c = conn.cursor()
    make_table(c)
    calc_stats(c, results.ndof)
    conn.commit()
    conn.close()
