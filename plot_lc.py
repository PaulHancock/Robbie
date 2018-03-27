#!/usr/bin/env python
__author__ = 'Paul Hancock'
__date__ = '27-03-2018'


import matplotlib.pyplot as plt
from astropy.table import Table
import sys
import os


def plot(n):
    row = cat.iloc[n]
    fname = 'plots/{0}.png'.format(row['uuid_0'])
    print fname,
    if os.path.exists(fname):
        print ".. skip"
        return
    fluxes = row[flux_cols]
    err_fluxes = row[err_flux_cols]
    plt.clf()
    plt.errorbar(range(len(fluxes)), fluxes, yerr=err_fluxes)
    plt.xlabel('Epoch')
    plt.ylabel('Flux Density (Jy/Beam)')
    s = 'm={0:5.3f}\nmd={1:4.2f}\nchisq={2:4.1f}'.format(row['m'], row['md'],row['chisq_peak_flux'])
    xlims = plt.xlim((0, len(fluxes)+5))
    ylims = plt.ylim()
    y = ylims[0] + (ylims[1]-ylims[0])*0.8
    plt.text(x=xlims[1]*0.8, y=y, s=s)
    plt.title('{0},{1}: {2}'.format(row['island_0'], row['source_0'], row['uuid_0']))
    plt.savefig(fname)
    print ".. done"
    return


fname = sys.argv[-1]
print "loading ",fname
cat = Table.read(fname).to_pandas()
flux_cols = [a for a in cat.columns if a.startswith('peak_flux')]
err_flux_cols = [a for a in cat.columns if a.startswith('err_peak_flux')]

for i in range(len(cat)):
    plot(i)