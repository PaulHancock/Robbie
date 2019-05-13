#! /usr/bin/env python

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
import astropy.units as u
import numpy as np
import matplotlib
from matplotlib import pyplot
from matplotlib.patches import Ellipse
import argparse
import sys

def grain_plot(fitsfile, plotfile):
    """

    """
    tab = Table.read(fitsfile)

    nepochs = np.max(tab['epoch'])*1.

    kwargs={'fontsize':14}

    cmap = pyplot.cm.plasma_r
    # define the bins and normalize
    bounds = np.linspace(3,15,13)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

    fig = pyplot.figure(figsize=(8,7))
    ax = fig.add_subplot(1,1,1)
    cax = ax.scatter(tab['ra'],tab['dec'],
                     c=tab['peak_flux']/tab['local_rms'],
                     cmap=cmap, norm=norm,
                     vmin=3, vmax=15,
                     zorder=100)
    for r in tab:
        ax.add_patch(Ellipse((r['ra'],r['dec']),
                             width=0.5, height=3, angle=r['epoch']/nepochs*360,
                             alpha=0.35,
                             color=cmap(norm(r['peak_flux']/r['local_rms'])),
                             edgecolor='none',
                             zorder=norm(r['peak_flux']/r['local_rms'])
                            ))

    cb = fig.colorbar(cax,ax=ax)

    cb.set_ticks(range(3,16,2))
    cb.set_label("SNR", **kwargs)
    ax.set_xlabel("RA J2000",**kwargs)
    ax.set_ylabel("Dec J2000", **kwargs)
    # flip the x axis so that RA increases to the left
    ax.set_xlim((ax.get_xlim()[1],ax.get_xlim()[0]))
    ax.grid()
    pyplot.savefig(plotfile)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Create a scatter grain plot")
    group1.add_argument("--in", dest='infile', type=str, default=None,
                        help="The input catalogue.")
    group1.add_argument("--plot", dest='plotfile', type=str, default=None,
                        help="output plot")

    results = parser.parse_args()

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()
    grain_plot(fitsfile=results.infile, plotfile=results.plotfile)
