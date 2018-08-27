#!/usr/bin/env python
__author__ = 'Paul Hancock'
__date__ = ''

import argparse
from astropy.table import Table
import pandas as pd
import numpy as np
from scipy import stats
import sys

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


def load_corrected_table(filename):
    tab = Table.read(filename)
    df = tab.to_pandas()
    flux_cols = [n for n in df.columns if n.startswith('peak')]
    err_flux_cols = [n for n in df.columns if n.startswith('err_peak')]

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
    tab2.write(outfile, overwrite=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Calculate variability stats")
    group1.add_argument("--infile", dest='infile', type=str, default=None,
                        help="The input catalogue.")
    group1.add_argument("--outfile", dest='outfile', type=str, default=None,
                        help="Output catalogue")
    group1.add_argument("--ndof", dest='ndof', type=float, default=None,
                        help="Effective number of degrees of freedom. Defualt: N=epochs-1")

    results = parser.parse_args()

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    df = load_corrected_table(results.infile)
    add_stats(df, results.outfile, ndof=results.ndof)
