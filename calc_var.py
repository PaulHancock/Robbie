#! python
__author__ = 'Paul Hancock'
__date__ = ''

from astropy.table import Table
import pandas as pd
import numpy as np
from scipy import stats

def chisq(series, fluxes=[], errs=[]):
    f = series[fluxes].values
    e = series[errs].values
    m = e == e
    f = f[m]
    e = e[m]
    mean = np.mean(f)
    val = np.sum((f - mean) ** 2) / np.sum(e)
    return val


def pval(series, fluxes=[], errs=[]):
    chi = chisq(series, fluxes, errs)
    e = series[errs].values
    m = e == e
    npts = len(m)
    if npts < 2:
        return 0
    p = stats.chi2.sf(chi, npts - 1)
    return max(p, 1e-150)


def load_table(filename, clip=slice(0,None)):
    """

    Parameters
    ----------
    filename
    clip

    Returns
    -------

    """
    tab = Table.read(filename)
    df = tab.to_pandas()

    # ignore the last few epochs as they have messed up images.
    flux_cols = [n for n in tab.colnames if n.startswith('peak')][clip]
    err_flux_cols = [n for n in tab.colnames if n.startswith('err_peak')][clip]


    mean_flux = df[flux_cols].mean(axis=1)
    df['mean_peak_flux'] = pd.Series(mean_flux, index=df.index)

    std_flux = df[flux_cols].std(axis=1)
    df['std_peak_flux'] = pd.Series(std_flux, index=df.index)

    chi_flux = df.apply(chisq, axis=1, fluxes=flux_cols, errs=err_flux_cols)
    df['chisq_peak_flux'] = pd.Series(chi_flux, index=df.index)

    pval_flux = df.apply(pval, axis=1, fluxes=flux_cols, errs=err_flux_cols)
    df['pval_peak_flux'] = pd.Series(pval_flux, index=df.index)

    tab2 = Table.from_pandas(df)
    outfile = filename.split('.')[0] + '_var.fits'
    tab2.write(outfile, overwrite=True)


if __name__ == '__main__':
    load_table('154MHz_flux_table.fits', clip=slice(0, 10))
    load_table('185MHz_flux_table.fits', clip=slice(0, 10))
