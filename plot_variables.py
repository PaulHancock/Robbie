#! /usr/bin/env python
from __future__ import print_function

from astropy.table import Table
import numpy as np
import matplotlib
from matplotlib import pyplot
from matplotlib.patches import Ellipse
import sqlite3
import argparse
import sys
import os


def plot_summary(cur, plotfile):
    """
    Plot
    """
    cur.execute("""SELECT pval_peak_flux, md, mean_peak_flux FROM stats""")
    rows = cur.fetchall()

    pval_peak_flux, md, mean_peak_flux = map(np.array, zip(*rows))

    kwargs = {'fontsize':14}
    fig = pyplot.figure(figsize=(5,8))


    ax = fig.add_subplot(1,1,1)
    cax = ax.scatter(md, np.log10(pval_peak_flux), c = np.log10(mean_peak_flux), cmap=matplotlib.cm.viridis_r)
    cb = fig.colorbar(cax,ax=ax)
    cb.set_label("log10(Peak flux in epoch 1) (Jy)", **kwargs)

    ax.set_ylim((-11,1.001))
    ax.set_xlim((-0.3,0.3))
    ax.set_ylabel("log(p_val)", **kwargs)
    ax.set_xlabel("Debiased modulation index ($m_d$)", **kwargs)
    ax.axhline(-3, c='k')
    ax.axvline(0.05, c='k')
    ax.text(0.1, -5, "variable", **kwargs)
    ax.fill_between([-0.3,0.05],-25, y2=2, color='k', alpha=0.2)
    ax.fill_betweenx([-3,2],0.05, x2=0.3, color='k', alpha=0.2)
    ax.text(-0.25, -5, "not variable", **kwargs)
    pyplot.savefig(plotfile)
    return


def plot_lc(cur):
    cur.execute("""SELECT DISTINCT uuid FROM sources""")
    sources = cur.fetchall()

    for src in sources:
        uuid = src[0]
        fname = 'plots/{0}.png'.format(uuid)
        print(fname, end='')
        if os.path.exists(fname):
            print(".. skip")
            continue

        cur.execute("""SELECT peak_flux, err_peak_flux, epoch FROM sources WHERE uuid=? ORDER BY epoch """, (uuid,))
        peak_flux, err_peak_flux, _ = map(np.array, zip(*cur.fetchall()))
        cur.execute("""SELECT m, md, chisq_peak_flux FROM stats WHERE uuid=? """, (uuid,))
        m, md, chisq_peak_flux = cur.fetchone()

        pyplot.clf()
        pyplot.errorbar(range(len(peak_flux)), peak_flux, yerr=err_peak_flux)
        pyplot.xlabel('Epoch')
        pyplot.ylabel('Flux Density (Jy/Beam)')
        s = 'm={0:5.3f}\nmd={1:4.2f}\nchisq={2:4.1f}'.format(m, md, chisq_peak_flux)
        xlims = pyplot.xlim((-0.5, len(peak_flux)+5))
        ylims = pyplot.ylim()
        y = ylims[0] + (ylims[1]-ylims[0])*0.8
        pyplot.text(x=xlims[1]*0.8, y=y, s=s)
        pyplot.title('{0}'.format(uuid))
        pyplot.savefig(fname)
        print(".. done")
    return



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Create a variability plot")
    group1.add_argument("--name", dest='name', type=str, default=None,
                        help="The input catalogue.")
    group1.add_argument("--plot", dest='plotfile', type=str, default=None,
                        help="output plot")

    results = parser.parse_args()

    if None in (results.name, results.plotfile):
        parser.print_help()
        sys.exit()

    conn = sqlite3.connect(results.name)
    cur = conn.cursor()
    plot_summary(cur=cur, plotfile=results.plotfile)
    plot_lc(cur=cur)
    conn.close()
