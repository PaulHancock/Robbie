#!/usr/bin/env python

"""
"""
__author__ = "Paul Hancock"
__date__ = '08/04/2016'

import argparse
from astropy.coordinates import SkyCoord, Angle, Latitude, Longitude
from astropy.table import Table, hstack
import astropy.units as u
import numpy as np
import os
from scipy import interpolate
import sys


def correct_tab(ref, tab, ref_cat=None, radius=2 / 60.):
    """
    Given a reference catalogue, and target catalogue:
    crossmatch the two catalogues
    calculate and remove the bulk offset
    crossmatch, model and remove smaller scale offsets (three times)
    return the corrected catalogue, and a map of crossmatches and separations
    """
    orig = tab.copy()
    # Allow for this to be precalculated
    if ref_cat is None:
        ref_cat = SkyCoord(ref['RAJ2000'], ref['DEJ2000'], unit=(u.degree, u.degree))
    
    # filter the catalog to use only high snr point sources
    # mask = np.where((tab['int_flux']/tab['peak_flux'] < 1.2) & (tab['peak_flux']/tab['local_rms']>10))
    mask = np.where(tab['peak_flux']/tab['local_rms'] > 10)
    # crossmatch the two catalogs
    cat = SkyCoord(tab['ra'][mask], tab['dec'][mask], unit=(u.degree, u.degree))
    idx, dist, _ = cat.match_to_catalog_sky(ref_cat)
    # accept only matches within radius
    distance_mask = np.where(dist.degree < radius) # this mask is into cat
    match_mask = idx[distance_mask] # this mask is into ref_cat
    # calculate the ra/dec offsets
    dra = ref_cat.ra.degree[match_mask] - cat.ra.degree[distance_mask]
    ddec = ref_cat.dec.degree[match_mask] - cat.dec.degree[distance_mask]
    
    # make a bulk correction as the first pass
    tab['ra'] += np.mean(dra)
    tab['dec'] += np.mean(ddec)
    
    # now do this again 3 more times but using the Rbf
    for i in range(3):
        # crossmatch the two catalogs
        cat = SkyCoord(tab['ra'][mask], tab['dec'][mask], unit=(u.degree, u.degree))
        idx, dist, _ = cat.match_to_catalog_sky(ref_cat)
        # accept only matches within radius
        distance_mask = np.where(dist.degree < radius) # this mask is into cat
        match_mask = idx[distance_mask] # this mask is into ref_cat
        if len(match_mask)<1:
            break
        # calculate the ra/dec offsets
        dra = ref_cat.ra.degree[match_mask] - cat.ra.degree[distance_mask]
        ddec = ref_cat.dec.degree[match_mask] - cat.dec.degree[distance_mask]
        
        # use the following to make some models of the offsets 
        dramodel  = interpolate.Rbf(cat.ra.degree[distance_mask], cat.dec.degree[distance_mask], dra, function='linear', smooth=3)
        ddecmodel = interpolate.Rbf(cat.ra.degree[distance_mask], cat.dec.degree[distance_mask], ddec, function='linear', smooth=3)

        # now adjust the origional source positions based on these models
        tab['ra'] += dramodel(tab['ra'], tab['dec'])
        tab['dec'] += ddecmodel(tab['ra'], tab['dec'])
        
    # final crossmatch to make the xmatch file
    cat = SkyCoord(tab['ra'][mask], tab['dec'][mask], unit=(u.degree, u.degree))
    idx, dist, _ = cat.match_to_catalog_sky(ref_cat)
    # accept only matches within radius
    distance_mask = np.where(dist.degree < radius) # this mask is into cat
    match_mask = idx[distance_mask] # this mask is into ref_cat
    # calculate the ra/dec offsets
    dra = ref_cat.ra.degree[match_mask] - cat.ra.degree[distance_mask]
    ddec = ref_cat.dec.degree[match_mask] - cat.dec.degree[distance_mask]
    ## cartesian!!
    separation = Table({'Separation':np.hypot(dra, ddec)}, names=('Separation',), dtype=(np.float32,))
    xmatch = hstack([orig[mask][distance_mask], ref[match_mask], separation])
    #xmatch['peak_flux'] = xmatch['peak_flux_1']
    #xmatch['local_rms'] = xmatch['local_rms_1']
    
    return tab, xmatch


def xmatch_one(infile, fout, refcat):
    print "reading", infile
    table = Table.read(infile)
    ref_table = Table.read(refcat)
    corr, xmatch = correct_tab(ref_table, table)
    print "writing", fout
    if os.path.exists(fout):
        os.remove(fout)
    xmatch.write(fout)
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("input/output files")
    group1.add_argument("--incat", dest='incat', type=str, default=None,
                        help="The catalogue to be corrected")
    group1.add_argument("--refcat", dest='refcat', type=str, default=None,
                        help="The reference catalogue")
    group1.add_argument("--xm", dest='xm', type=str, default=None,
                        help='Output file for the crossmatch between the reference and source catalogue.')
    group2 = parser.add_argument_group("catalog column names")
    group2.add_argument("--ra1", dest='ra1', type=str, default='ra',
                        help="The column name for ra  (degrees) for source catalogue.")
    group2.add_argument("--dec1", dest='dec1', type=str, default='dec',
                        help="The column name for dec (degrees) for source catalogue.")
    group2.add_argument("--ra2", dest='ra2', type=str, default='RAJ2000',
                        help="The column name for ra  (degrees) for reference catalogue.")
    group2.add_argument("--dec2", dest='dec2', type=str, default='DEJ2000',
                        help="The column name for dec (degrees) for reference catalogue.")
    group3 = parser.add_argument_group("Other")
    group3.add_argument('--smooth', dest='smooth', default=10, type=float,
                        help="Smoothness parameter to give to the radial basis function (default = 10 arcmin)")

    results = parser.parse_args()
    if len(sys.argv) <= 3:
        print("correct_astrometry.py refcat incat outcat")
        sys.exit(1)

    refcat, incat, outcat = sys.argv[1:]
    xmatch_one(incat, outcat, refcat)