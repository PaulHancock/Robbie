#! python
from __future__ import print_function
import numpy as np
import pandas as pd
import sys
import os

__author__ = 'Paul Hancock'
__date__ = ''


def add_stats_pd(infile, outfile):
    tab = pd.read_csv(infile)
    fluxcols = [c for c in tab.columns if c.startswith("peak_flux")]
    errcols = [c for c in tab.columns if c.startswith("err_peak_flux")]
    rmscols = [c for c in tab.columns if c.startswith("local_rms")]


    # replace the blank error entries with the local rms
    for i in range(len(rmscols)):
        tab[errcols[i]][np.isnan(tab[errcols[i]])] = tab[rmscols[i]]

    mean = tab[fluxcols].mean(axis=1)
    std = tab[fluxcols].std(axis=1)
    n = len(tab)

    tab['m'] = std/mean
    tab['chisq'] = np.sum((tab[fluxcols].values - mean.values[:, None])**2 / tab[errcols].values**2, axis=1)

    desc = 1./n * (np.sum((tab[fluxcols].values - mean.values[:, None])**2, axis=1) - np.sum(tab[errcols]**2, axis=1))
    md = 1/mean * np.sqrt(desc.abs())
    tab['md'] = md * ( (desc < 0) * -2 +1 )

    if os.path.exists(outfile):
        os.remove(outfile)
    tab.to_csv(outfile)
    return


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage {0} infile outfile".format(sys.argv[0]))
        sys.exit()
    infile, outfile = sys.argv[-2:]
    add_stats_pd(infile, outfile)
