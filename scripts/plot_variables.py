#! /usr/bin/env python
from __future__ import print_function

from astropy.table import Table
import dateutil
import dateutil.parser
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.cm import viridis_r
import argparse
import sys
import os
import multiprocessing as mp
import datetime

__author__ = ["Paul Hancock"]
__date__ = '2022/03/24'


def plot_summary_table(filename, plotfile):
    """
    Create a summary plot for all sources, identifying which are likely to be variable.

    parameters
    ----------
    filename : str
        Input table filename
    plotfile : str
        Filename for the output plot file
    """
    tab = Table.read(filename)
    pval_peak_flux = tab['pval_peak_flux_ks']
    md = tab['md']
    mean_peak_flux = tab['mean_peak_flux']

    kwargs = {'fontsize':14}
    fig = plt.figure(figsize=(5,8))


    ax = fig.add_subplot(1,1,1)
    cax = ax.scatter(md, np.log10(pval_peak_flux), c = np.log10(mean_peak_flux), cmap=viridis_r)
    cb = fig.colorbar(cax,ax=ax)
    cb.set_label("log10(Peak flux in epoch 1) (Jy)", **kwargs)

    ax.set_ylim((-11,1.001))
    ax.set_xlim((-0.3,0.3))
    ax.set_ylabel("log(p_val_ks)", **kwargs)
    ax.set_xlabel("Debiased modulation index ($m_d$)", **kwargs)
    ax.axhline(-3, c='k')
    ax.axvline(0.05, c='k')
    ax.text(0.1, -5, "variable", **kwargs)
    ax.fill_between([-0.3,0.05],-25, y2=2, color='k', alpha=0.2)
    ax.fill_betweenx([-3,2],0.05, x2=0.3, color='k', alpha=0.2)
    ax.text(-0.25, -5, "not variable", **kwargs)
    plt.savefig(plotfile)
    return


def plot_lc_table(flux_table, stats_table, start=0, stride=1, plot_dir="plots", dates=False):
    """
    Create individual light curve plots.
    Each plot is saved to plots/uuid.png

    parameters
    ----------
    flux_table : str
        Filename of the flux table

    stats_table : str
        Filename of the stats table

    start : int
        Starting row (default=0)

    stride : int
        Process every Nth row of the table. Default =1
    
    dates : bool = False
        If true then use dates for the x-axis value/format
    """
    flux_tab = Table.read(flux_table)
    stats_tab = Table.read(stats_table)
    epochs = [a for a in flux_tab.colnames if a.startswith('epoch')]
    fluxes = [a for a in flux_tab.colnames if a.startswith('peak_flux')]
    err_fluxes = [a for a in flux_tab.colnames if a.startswith('err_peak_flux')]
    for row in flux_tab[start::stride]:
        fname = '{0}/{1}.png'.format(plot_dir, row['uuid'])
        print(fname, end='')
        if os.path.exists(fname):
            print(" ... skip")
            continue
        srow = stats_tab[stats_tab['uuid'] == row['uuid']]

        # Sort date by date
        mask = np.where(['None' not in row[a] for a in epochs])[0]
        if len(mask) == 0 :
            print(" ... no data")
            continue
        epoch_mask = list(np.choose(mask, epochs))
        flux_mask = list(np.choose(mask, fluxes))
        err_flux_mask = list(np.choose(mask, err_fluxes))
        if dates:
            epoch_times = [datetime.datetime.strptime(a, "%Y-%m-%dT%H:%M:%S") for a in list(row[epochs][epoch_mask])]
        else:
            epoch_times = list(range(len(epoch_mask)))

        # Annotate with stats
        s = f"m={srow['m'][0]:5.3f}\nmd={srow['md'][0]:4.2f}\nchisq={srow['chisq_peak_flux'][0]:4.1f}"

        # convert epochs to datetime objects
        fig, ax = plt.subplots()
        ax.errorbar(epoch_times,
                    list(row[fluxes][flux_mask]), 
                    yerr=row[err_fluxes][err_flux_mask], 
                    label=s)
        ax.set_ylabel('Flux Density (Jy/Beam)')
        ax.set_xlabel('Epoch')
        ax.set_title('{0}'.format(row['uuid']))
        ax.legend()
        if dates:
            fig.autofmt_xdate()
            ax.fmt_xdata = mdates.DateFormatter("%Y-%m-%dT%H:%M:%S")
        plt.savefig(fname, bbox_inches='tight')
        plt.close(fig)
        print(" ... done")
    return


def plot_lc_table_parallel(flux_table, stats_table, light_curve_dir, dates, nprocs=1, debug=False):
    """
    parameters
    ----------
    flux_table : str
        Filename of the flux table

    stats_table : str
        Filename of the stats table
    
    light_curve_dir : str
        Location to store the plots

    dates : bool
        Whether to use dates for the plot horizontal axis (True) or epochs (False)

    nprocs : int
        Number of processes to use simultaneously
    """
    pool = mp.Pool(nprocs)
    results = []
    for i in range(nprocs):
        r=pool.apply_async(
                            plot_lc_table,
                            args=[
                                flux_table,
                                stats_table
                            ],
                            kwds={
                                'start':i,
                                'stride':nprocs,
                                'plot_dir':light_curve_dir,
                                'dates':dates
                            }
                       )
        if debug:
            r.get()
        else:
            results.append(r)
    pool.close()
    pool.join()
    if not debug:
        # This forces any raised exceptions within the apply_async to be re-raised here
        for r in results:
            r.get()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Create a variability plot")
    group1.add_argument("--ftable", dest='ftable', type=str, default=None,
                        help="flux table")
    group1.add_argument("--stable", dest='stable', type=str, default=None,
                        help="stats table")
    group1.add_argument("--plot", dest='plotfile', type=str, default=None,
                        help="output plot")
    group1.add_argument("--all", dest='all', action='store_true', default=False,
                        help="Also plot individual light curves. Default:False")
    group1.add_argument("--lc_dir", dest='light_curve_dir', type=str, default="plots",
                        help="The light curve plots output directory")
    group1.add_argument("--dates", dest='dates', action='store_true', default=False,
                        help="Individual plots have date on the horizontal axis.")
    group1.add_argument("--cores", dest='cores', type=int, default=None,
                        help="Number of cores to use: Default all")
    group1.add_argument("--debug", dest='debug', action='store_true', default=False,
                        help="Use debug mode")

    results = parser.parse_args()

    if results.cores is None:
        results.cores = mp.cpu_count()

    if results.ftable or results.stable:
        if not (results.ftable and results.stable):
            print("ERROR: --stable and --ftable are both required, only one supplied.")
        plot_summary_table(results.stable, results.plotfile)
        if results.all:
            plot_lc_table_parallel(results.ftable, 
                                   results.stable, 
                                   results.light_curve_dir, 
                                   results.dates,
                                   nprocs=results.cores,
                                   debug=results.debug)
    else:
        parser.print_help()
        sys.exit()

