#!/usr/bin/env python

from __future__ import print_function
__author__ = 'Paul Hancock'
__date__ = '2019/08/18'

import argparse
import sqlite3
import numpy as np
from scipy import stats
import sys
import astropy.table
from astropy.table import Table, Column
from join_catalogues import write_table
import multiprocessing as mp


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

    ndof : int or None
        Number of degrees of freedom for the light curves. None -> npts-1
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
            std = 0.
            m = 0.
            chisq = 0.
        else:
            # modulation index
            mean = np.mean(fluxes[mask])
            std = np.std(fluxes[mask])
            m = std/mean
            # chi squared
            chisq = np.sum((fluxes[mask] - mean)**2 / err[mask]**2)
            # pvalue
            if ndof is None:
                ndof = max(1,npts - 1)
            else:
                ndof = max(1,ndof)
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


def calc_stats_table(filename, ndof=None, start=0, stride=1):
    """
    Compute various stats for each light curve
    Start at a given row number and process only nrows from this table.

    parameters
    ----------
    filename : str
        The filename of the table to be read

    ndof : int or None
        Number of degrees of freedom for the light curves. None -> npts-1

    start : int
        Starting row (default=0)

    stride : int
        Process every Nth row of the table. Default =1

    return
    ------
    tab : `astropy.table.Table`
        A table of stats
    """
    # print("loading {0}".format(filename))
    tab = Table.read(filename)
    # print("done")
    start = max(start, 0)
    start = min(start, len(tab))
    tab = tab[start::stride]
    # print("Using {0} rows, with start {1} and stride {2}".format(len(tab), start, stride))

    flux_cols = [a for a in tab.colnames if a.startswith('peak_flux')]
    err_flux_cols = [a for a in tab.colnames if a.startswith('err_peak_flux')]
    src_stats = np.zeros(shape=(len(tab), 6))
    for i, row in enumerate(tab):
        fluxes = np.array(list(row[flux_cols]))
        err = np.array(list(row[err_flux_cols]))
        mask = np.where(err>0)
        npts = len(mask[0])
        
        if npts < 2:
            pval = 0.
            md = 0.
            mean = 0.
            std = 0.
            m = 0.
            chisq = 0.
        else:
            # modulation index
            mean = np.mean(fluxes[mask])
            std = np.std(fluxes[mask])
            m = std/mean
            # chi squared
            chisq = np.sum((fluxes[mask] - mean)**2 / err[mask]**2)
            # pvalue
            if ndof is None:
                ndof = max(1,npts - 1)
            else:
                ndof = max(1,ndof)
            pval = stats.chi2.sf(chisq, ndof)
            pval = max(pval, 1e-10)
            # debiased modulation index
            desc = np.sum((fluxes[mask] - mean)**2) - np.sum(err[mask]**2)
            #print(mean, desc, npts)
            md = 1./mean * np.sqrt(np.abs(desc)/npts)
            if desc < 0:
                md *= -1
        src_stats[i, :] = [mean, std, m,  md, chisq, pval]

    stats_tab = Table()
    stats_tab.add_column(tab['uuid'])
    stats_tab.add_column(Column(src_stats[:,0], name='mean_peak_flux'))
    stats_tab.add_column(Column(src_stats[:,1], name='std_peak_flux'))
    stats_tab.add_column(Column(src_stats[:,2], name='m'))
    stats_tab.add_column(Column(src_stats[:,3], name='md'))
    stats_tab.add_column(Column(src_stats[:,4], name='chisq_peak_flux'))
    stats_tab.add_column(Column(src_stats[:,5], name='pval_peak_flux'))
    return stats_tab


def calc_stats_table_parallel(filename, ndof=None, nprocs=1):
    """
    Compute various stats for each light curve using multiple cores if available

    parameters
    ----------
    filename : str
        The filename of the table to be read

    ndof : int or None
        Number of degrees of freedom for the light curves. None -> npts-1

    nprocs : int
        Number of processes to use simultaneously

    return
    ------
    tab : `astropy.table.Table`
        A table of stats, not necessarily in the same order as the input!
    """
    results = []

    def collect_result(result):
        results.append(result)

    pool = mp.Pool(nprocs)
    for i in range(nprocs):
        pool.apply_async(calc_stats_table,
                         args=[filename, ndof],
                         kwds={'start':i, 'stride':nprocs},
                         callback=collect_result)
    pool.close()
    pool.join()
    stats_tab = astropy.table.vstack(results)
    return stats_tab


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Calculate variability stats")
    group1.add_argument("--dbname", dest='db', type=str, default=None,
                        help="Database filename.")
    group1.add_argument("--table", dest='table', type=str, default=None,
                        help="Table filename. [requires --out]")
    group1.add_argument("--out", dest='out', type=str, default=None,
                        help="Output filename.")
    group1.add_argument("--ndof", dest='ndof', type=float, default=None,
                        help="Effective number of degrees of freedom. Defualt: N=epochs-1")
    group1.add_argument("--cores", dest='cores', type=int, default=None,
                        help="Number of cores to use: Default all")
    
    results = parser.parse_args()

    if results.cores is None:
        results.cores = mp.cpu_count()

    if results.db:
        conn = sqlite3.connect(results.db)
        c = conn.cursor()
        make_table(c)
        calc_stats(c, results.ndof)
        conn.commit()
        conn.close()
    elif results.table:
        if not results.out:
            print("ERROR: --table requires --out to be set")
            sys.exit(1)
        tab = calc_stats_table_parallel(results.table, results.ndof, nprocs=results.cores)
        write_table(tab, results.out)
    else:
        parser.print_help()
        sys.exit(1)

