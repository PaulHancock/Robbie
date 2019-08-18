#!/usr/bin/env python
__author__ = 'Paul Hancock'
__date__ = ''

import argparse
from astropy.table import Table
import sqlite3
import pandas as pd
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
        fluxes = np.array(cur.fetchall())
        cur.execute("SELECT err_peak_flux FROM sources WHERE uuid=? ORDER BY epoch", s)
        err = np.array(cur.fetchall())

        # don't include fit errors in the stats calculation
        mask = np.where(err>0)
        mean = np.mean(fluxes[mask])
        std = np.std(fluxes[mask])
        m = std/mean
        chisq = np.sum(( fluxes[mask] - mean)**2 / err[mask]**2)
        npts = len(mask)
        if npts <2:
            pval = 0
        else:
            if ndof is None:
                ndof = npts - 1
            pval = stats.chi2.sf(chisq, ndof)
            pval = max(pval, 1e-10)
        # debiased modulation index
        desc = np.sum((fluxes[mask] - mean)**2) - np.sum(err[mask]**2)
        md = 1./mean * np.sqrt(np.abs(desc)/len(mask))
        md = md*((desc > 0)*2 - 1)
        # add all to the stats table
        cur.execute("""INSERT INTO stats(uuid, mean_peak_flux, std_peak_flux, m, md, chisq_peak_flux, pval_peak_flux)
        VALUES (?,?,?,?,?,?,?)""",
                    (s[0], mean, std, m,  md, chisq, pval))
    return


def chisq(series, fluxes=[], errs=[]):
    f = series[fluxes].values
    e = series[errs].values
    m = e == e
    f = f[m]
    e = e[m]
    mean = np.mean(f)
    val = np.sum((f - mean) ** 2 / e**2)
    return val


def pval(series, fluxes=[], errs=[], ndof=None):
    chi = chisq(series, fluxes, errs)
    e = series[errs].values
    m = e == e
    npts = len(m)
    if npts < 2:
        return 0
    if ndof is None:
        ndof = npts - 1
    p = stats.chi2.sf(chi, ndof)
    return max(p, 1e-10)


def modulation_index(series, fluxes=[], errs=[]):
    f = series[fluxes]
    return np.std(f)/np.mean(f)


def debias_modulation_index(series, fluxes=[], errs=[]):
    f = series[fluxes].values
    e = series[errs].values
    m = e == e
    f = f[m]
    e = e[m]
    mean = np.mean(f)
    desc = np.sum((f - mean)**2) - np.sum(e**2)
    md = 1./mean * np.sqrt(np.abs(desc)/len(m))
    md = md*((desc > 0)*2 - 1)

    return md


def norm(series, fluxes=[]):
    f = series[fluxes].values
    mean = np.mean(f)
    return f/mean


def load_table(filename, correct=False):
    """
    Load the given table, and apply corrections if requested.

    parameters
    ==========
    filename : string
      The name of the file to read

    correct : bool
      Compute and remove the 'mean light curve' if correct=True, 
      and there are more than 5 rows in the table.

    return
    ======
    df : pandas.df
      A data frame of the given table.
    """
    tab = Table.read(filename)
    df = tab.to_pandas()

    ## The following corrections are not required unless your data has some
    ## real problems with instrumental variability.
    if correct:
        flux_cols = [n for n in df.columns if n.startswith('peak')]
        err_flux_cols = [n for n in df.columns if n.startswith('err_peak')]
        if len(df) > 5:
            mean_fluxes = df[flux_cols].apply(norm, axis=1, fluxes=flux_cols)
            mean_lc = mean_fluxes.median(axis=0)
            df[flux_cols] = df[flux_cols] / mean_lc
            df[err_flux_cols] = df[err_flux_cols].divide(mean_lc.values)
    return df


def add_stats(df, outfile=None, ndof=None):
    """

    Parameters
    ----------
    tab : pandas.Dataframe

    outfile: string

    ndof: float

    Returns
    -------

    """
    flux_cols = [n for n in df.columns if n.startswith('peak')]
    err_flux_cols = [n for n in df.columns if n.startswith('err_peak')]

    if ndof is None:
        ndof = len(flux_cols) - 1

    mean_flux = df[flux_cols].mean(axis=1)
    df['mean_peak_flux'] = pd.Series(mean_flux, index=df.index)

    std_flux = df[flux_cols].std(axis=1)
    df['std_peak_flux'] = pd.Series(std_flux, index=df.index)

    chi_flux = df.apply(chisq, axis=1, fluxes=flux_cols, errs=err_flux_cols)
    df['chisq_peak_flux'] = pd.Series(chi_flux, index=df.index)

    pval_flux = df.apply(pval, axis=1, fluxes=flux_cols, errs=err_flux_cols, ndof=ndof)
    df['pval_peak_flux'] = pd.Series(pval_flux, index=df.index)

    m = df.apply(modulation_index, axis=1, fluxes=flux_cols, errs=err_flux_cols)
    df['m'] = pd.Series(m, index=df.index)

    md = df.apply(debias_modulation_index, axis=1, fluxes=flux_cols, errs=err_flux_cols)
    df['md'] = pd.Series(md, index=df.index)


    tab2 = Table.from_pandas(df)
    fmt = outfile.split('.')[-1]
    if 'vot' in fmt:
        fmt = 'votable'
    tab2.write(outfile, overwrite=True, format=fmt)


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
