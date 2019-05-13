#! /usr/bin/env python

from astropy.table import Table
import numpy as np
import matplotlib
from matplotlib import pyplot
from matplotlib.patches import Ellipse
import argparse
import sys


def plot(fitsfile, plotfile):
    """
    Plot
    """

    tab = Table.read(fitsfile)

    kwargs = {'fontsize':14}
    fig = pyplot.figure(figsize=(5,8))

    ax = fig.add_subplot(1,1,1)
    cax = ax.scatter(tab['md'], np.log10(tab['pval_peak_flux']), c = np.log10(tab['peak_flux_1']), cmap=matplotlib.cm.viridis_r)
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Create a variability plot")
    group1.add_argument("--in", dest='infile', type=str, default=None,
                        help="The input catalogue.")
    group1.add_argument("--plot", dest='plotfile', type=str, default=None,
                        help="output plot")

    results = parser.parse_args()

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()
    plot(fitsfile=results.infile, plotfile=results.plotfile)
