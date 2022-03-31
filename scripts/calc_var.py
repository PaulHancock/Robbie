#! /usr/bin/env python3


__author__ = 'Paul Hancock'
__date__ = '2019/08/18'

import argparse
import numpy as np
from scipy import stats
import sys
import astropy.table
from astropy.table import Table, Column
from join_catalogues import write_table
import multiprocessing as mp


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
    src_stats = np.zeros(shape=(len(tab), 7))
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
            pval_ks = 0.
        else:
            # modulation index
            mean = np.mean(fluxes[mask])
            std = np.std(fluxes[mask])
            m = std/mean
            # chi squared
            chisq = np.sum((fluxes[mask] - mean)**2 / err[mask]**2)
            # pvalue from chi squared
            if ndof is None:
                ndof = max(1,npts - 1)
            else:
                ndof = max(1,ndof)
            pval = stats.chi2.sf(chisq, ndof)
            pval = max(pval, 1e-10)

            # Pvalue based on distribution of Z score
            Z = (fluxes[mask] - mean)/ err[mask]
            pval_ks = stats.kstest(Z,'norm').pvalue
            pval_ks = max(pval_ks, 1e-10)

            # debiased modulation index
            desc = np.sum((fluxes[mask] - mean)**2) - np.sum(err[mask]**2)
            #print(mean, desc, npts)
            md = 1./mean * np.sqrt(np.abs(desc)/npts)
            if desc < 0:
                md *= -1
        src_stats[i, :] = [mean, std, m,  md, chisq, pval, pval_ks]

    stats_tab = Table()
    stats_tab.add_column(tab['uuid'])
    stats_tab.add_column(Column(src_stats[:,0], name='mean_peak_flux'))
    stats_tab.add_column(Column(src_stats[:,1], name='std_peak_flux'))
    stats_tab.add_column(Column(src_stats[:,2], name='m'))
    stats_tab.add_column(Column(src_stats[:,3], name='md'))
    stats_tab.add_column(Column(src_stats[:,4], name='chisq_peak_flux'))
    stats_tab.add_column(Column(src_stats[:,5], name='pval_peak_flux_chisq'))
    stats_tab.add_column(Column(src_stats[:,6], name='pval_peak_flux_ks'))
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

    if results.table:
        if not results.out:
            print("ERROR: --table requires --out to be set")
            sys.exit(1)
        tab = calc_stats_table_parallel(results.table, results.ndof, nprocs=results.cores)
        write_table(tab, results.out)
    else:
        parser.print_help()
        sys.exit(1)

