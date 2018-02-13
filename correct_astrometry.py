#! python

"""
"""
__author__ = "Paul Hancock"
__date__ = '08/04/2016'

import astropy
from astropy.io import fits
from astropy.coordinates import SkyCoord, Angle, Latitude, Longitude
from astropy.table import Table, hstack
import astropy.units as u
from glob import glob
import matplotlib
from matplotlib import pyplot
import numpy as np
import os
import scipy
from scipy import interpolate
import sys


def correct_tab(ref, tab, ref_cat=None):
    """
    Create a model and apply ionospheric corrections to a table, given a reference catalogue.
    """

    # Allow for this to be precalculated
    if ref_cat is None:
        ref_cat = SkyCoord(ref['ra'], ref['dec'], unit=(u.degree, u.degree))
    
    # filter the catalog to use only high snr point sources
    mask = np.where( (tab['int_flux']/tab['peak_flux'] < 1.2) & (tab['peak_flux']/tab['local_rms']>10) )
    cat = SkyCoord(tab['ra'][mask], tab['dec'][mask], unit=(u.degree, u.degree))

    # crossmatch the two catalogs
    # calculate the ra/dec offsets
    idx, dist, _ = cat.match_to_catalog_sky(ref_cat)
    # accept only matches within 2 arcmin
    distance_mask = np.where(dist.arcmin < 2) # this mask is into cat
    match_mask = idx[distance_mask] # this mask is into ref_cat

    dra = ref_cat.ra.degree[match_mask] - cat.ra.degree[distance_mask]
    ddec = ref_cat.dec.degree[match_mask] - cat.dec.degree[distance_mask]

    # use the following to make some models of the offsets 
    dramodel  = interpolate.Rbf(cat.ra.degree[distance_mask], cat.dec.degree[distance_mask], dra, function='linear', smooth=3)
    ddecmodel = interpolate.Rbf(cat.ra.degree[distance_mask], cat.dec.degree[distance_mask], ddec, function='linear', smooth=3)

    # now adjust the origional source positions based on these models
    tab['ra'] += dramodel(tab['ra'], tab['dec'])
    tab['dec'] += ddecmodel(tab['ra'], tab['dec'])

    return tab


def correct_bulk_offset(ref, tab, ref_cat=None, radius=2/60.):
    """
    
    """
    orig = tab.copy()
    # Allow for this to be precalculated
    if ref_cat is None:
        ref_cat = SkyCoord(ref['ra'], ref['dec'], unit=(u.degree, u.degree))
    
    # filter the catalog to use only high snr point sources
    mask = np.where( (tab['int_flux']/tab['peak_flux'] < 1.2) & (tab['peak_flux']/tab['local_rms']>10) )
    
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
    #print "bulk", np.mean(dra),np.mean(ddec)
    
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
    separation = Table({'Separation': np.hypot(dra, ddec)}, names=('Separation',), dtype=(np.float32,))
    xmatch = hstack([orig[mask][distance_mask], ref[match_mask], separation])
    
    return tab, xmatch


def correct_tab2(ref, tab, ref_cat=None, radius=2/60.):
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
        ref_cat = SkyCoord(ref['ra'], ref['dec'], unit=(u.degree, u.degree))
    
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
    #print "bulk", np.mean(dra),np.mean(ddec)
    
    # now do this again 3 more times but using the Rbf
    for i in range(3):
        # crossmatch the two catalogs
        cat = SkyCoord(tab['ra'][mask], tab['dec'][mask], unit=(u.degree, u.degree))
        idx, dist, _ = cat.match_to_catalog_sky(ref_cat)
        # accept only matches within radius
        distance_mask = np.where(dist.degree < radius) # this mask is into cat
        match_mask = idx[distance_mask] # this mask is into ref_cat
        # calculate the ra/dec offsets
        dra = ref_cat.ra.degree[match_mask] - cat.ra.degree[distance_mask]
        ddec = ref_cat.dec.degree[match_mask] - cat.dec.degree[distance_mask]
        #print " iter",i, np.mean(dra),np.mean(ddec)
        
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
    separation = Table({'Separation':np.hypot(dra,ddec)}, names=('Separation',), dtype=(np.float32,))
    xmatch = hstack([orig[mask][distance_mask], ref[match_mask], separation])
    
    return tab, xmatch


def xmatch_one(infile, fout, ref_cat, ref_table):
    print "reading", infile
    table = Table.read(infile)
    corr, xmatch = correct_tab2(ref_table, table, ref_cat)
    print "writing", fout
    if os.path.exists(fout):
        os.remove(fout)
    xmatch.write(fout)
    return


def correct_all(refname):
    directories = glob('1*')
    ref_table = Table.read(refname)
    ref_cat = SkyCoord(ref_table['RAJ2000'],ref_table['DECJ2000'], unit=(u.degree,u.degree))
    for d in directories:
        print d
        f = "{0}/{0}_deeper-I-image_comp.fits".format(d)
        print "reading", f
        table = Table.read(f)
        cat = SkyCoord(table['ra'],table['dec'], unit=(u.degree,u.degree))
        corrected = correct_tab(None, table, ref_cat)
        # do a second round of corrections
        corrected = correct_tab(ref_table, corrected, ref_cat)
        fout = "{0}/{0}_deeper-I-image_comp_icorr.fits".format(d)
        print "writing", fout
        corrected.write(fout)
    return

        
def xmatch_all(refname, flist):
    directories = glob('1*')
    ref_table = Table.read(refname)
    ref_cat = SkyCoord(ref_table['RAJ2000'], ref_table['DEJ2000'], unit=(u.degree, u.degree))
    for d in directories:
        print d
        f = "{0}/{0}_deeper-I-image_comp.fits".format(d)
        fout = "{0}/{0}_xm.fits".format(d)
        xmatch_one(f, fout, ref_cat)
    return


def make_plots():
    directories = glob('1*')
    for d in directories:
        print d
        raw_data = fits.open(d+'/xmatched.fits')[1].data
        # filter the data to only include SNR>20 sources
        flux_mask = np.where(raw_data['peak_flux']/raw_data['local_rms']>20)
        data = raw_data[flux_mask]
        #calculate the offsets in the ra/dec directions
        raw_catalog = Longitude(raw_data['ra'],unit=u.degree), Latitude(raw_data['dec'], unit=u.degree)
        catalog = Longitude(data['ra'],unit=u.degree), Latitude(data['dec'], unit=u.degree)

        raw_reference = Longitude(raw_data['RAJ2000'],unit=u.degree), Latitude(raw_data['DECJ2000'], unit=u.degree)
        reference = Longitude(data['RAJ2000'],unit=u.degree), Latitude(data['DECJ2000'], unit=u.degree)
        dra = (reference[0]-catalog[0]).arcsec
        ddec = (reference[1]-catalog[1]).arcsec

        xmin,xmax = np.nanmin(catalog[0].degree), np.nanmax(catalog[0].degree)
        ymin,ymax = np.nanmin(catalog[1].degree), np.nanmax(catalog[1].degree)

        # make some plots of the offsets
        fig = pyplot.figure(figsize=(8,6))
        ax = fig.add_subplot(1,2,1)
        cax = ax.scatter(catalog[0].degree,catalog[1].degree,c=dra, edgecolor='',vmin=-60, vmax=60)
        ax.set_xlim((10,40))
        ax.set_ylim((-35,-15))
        ax.set_title("DRA")
        cbar = fig.colorbar(cax, orientation='horizontal')

        ax = fig.add_subplot(1,2,2)
        cax = ax.scatter(catalog[0].degree,catalog[1].degree,c=ddec, edgecolor='',vmin=-60, vmax=60)
        ax.set_xlim((10,40))
        ax.set_ylim((-35,-15))
        ax.set_title("DDec")
        cbar = fig.colorbar(cax, orientation='horizontal')
        pyplot.savefig('plots/offsets_'+d+'.png')
        pyplot.close()

        # Now make a model to capture the offsets.
        dramodel = interpolate.Rbf(catalog[0].degree,catalog[1].degree, dra, function='linear', smooth=3)
        ddecmodel = interpolate.Rbf(catalog[0].degree,catalog[1].degree, ddec, function='linear', smooth=3)

        # This grid approach is only needed for plotting and for some stupid reason I have to flip/transpose the data
        # in order to get the figure to plot properly
        grid = np.mgrid[xmin:xmax:0.8,ymin:ymax:0.8]
        drainterpolated = np.flipud(dramodel(grid[0],grid[1]).T)
        ddecinterpolated = np.flipud(ddecmodel(grid[0],grid[1]).T)

        extent = [xmin,xmax,ymin,ymax]
        kwargs = {'extent':extent,'interpolation':'nearest', 'aspect':2}
        fig = pyplot.figure(figsize=(9,6))
        ax = fig.add_subplot(1,4,1)
        cax = ax.imshow(drainterpolated, vmin=-30, vmax=30, **kwargs)
        ax.set_xlim((10,40))
        ax.set_ylim((-35,-15))
        ax.set_title('dra')
        cbar = fig.colorbar(cax, orientation='horizontal')

        ax = fig.add_subplot(1,4,2)
        cax = ax.imshow(ddecinterpolated, vmin=-30, vmax=30, **kwargs)
        cax = ax.scatter(catalog[0].degree, catalog[1].degree, 
           c=ddecmodel(catalog[0].degree, catalog[1].degree), 
           vmin=-30, vmax=30,
           edgecolors='', s=50)
        ax.set_xlim((10,40))
        ax.set_ylim((-35,-15))
        ax.set_title('ddec')
        cbar = fig.colorbar(cax, orientation='horizontal')
        ax.plot(catalog[0].degree, catalog[1].degree, 'k,')

        ax = fig.add_subplot(1,4,3)
        cax = ax.imshow(np.hypot(drainterpolated,ddecinterpolated), vmin=0, vmax=30, **kwargs)
        ax.set_xlim((10,40))
        ax.set_ylim((-35,-15))
        ax.set_title('offset')
        cbar = fig.colorbar(cax, orientation='horizontal')

        ax = fig.add_subplot(1,4,4)
        cax = ax.imshow(np.degrees(np.arctan2(drainterpolated,ddecinterpolated)), vmin=-180, vmax=180, **kwargs)
        ax.set_xlim((10,40))
        ax.set_ylim((-35,-15))
        ax.set_title('angle')
        cbar = fig.colorbar(cax, orientation='horizontal')

        pyplot.savefig('plots/astrometric_model'+d+'.png')
        pyplot.close()

        # compute the corrected ra/dec *offsets*
        # this should be ~0 if we have a good model
        raw_dra = (raw_reference[0]-raw_catalog[0]).arcsec
        raw_ddec = (raw_reference[1]-raw_catalog[1]).arcsec
        corrected_ra = raw_dra - dramodel(raw_catalog[0].degree, raw_catalog[1].degree)
        corrected_dec = raw_ddec - ddecmodel(raw_catalog[0].degree, raw_catalog[1].degree)

        # more plotting
        fig = pyplot.figure(figsize=(9,6))
        ax = fig.add_subplot(1,2,1)
        cax = ax.scatter(raw_catalog[0].degree,raw_catalog[1].degree,
           c=corrected_ra, edgecolor='',
           vmin=-60, vmax=60)
        ax.set_xlim((10,40))
        ax.set_ylim((-35,-15))
        ax.set_title("DRA")
        cbar = fig.colorbar(cax, orientation='horizontal')

        ax = fig.add_subplot(1,2,2)
        cax = ax.scatter(raw_catalog[0].degree,raw_catalog[1].degree,
           c=corrected_dec, edgecolor='',
           vmin=-60, vmax=60)
        ax.set_xlim((10,40))
        ax.set_ylim((-35,-15))
        ax.set_title("DDec")
        cbar = fig.colorbar(cax, orientation='horizontal')

        pyplot.savefig('plots/astrometric_corrected'+d+'.png')
        pyplot.close()


def xmatch_freq(freq):
    refname = 'median_{0}MHz_comp.fits'.format(freq)
    ref_table = Table.read(refname)
    ref_cat = SkyCoord(ref_table['ra'], ref_table['dec'], unit=(u.degree, u.degree))
    flist = 'K2_{0}MHz.dat'.format(freq)
    for line in open(flist).readlines():
        f_base = line.split()[0].split('.')[0]
        f_in = "{0}_comp.fits".format(f_base)
        f_out = "{0}_xm.fits".format(f_base)
        print "{0} -> {1}".format(f_in, f_out)
        xmatch_one(infile=f_in, fout=f_out, ref_cat=None, ref_table=ref_table)


if __name__ == '__main__':
    if 'xmatch' in sys.argv:
        freq = sys.argv[-1]
        if freq in ['154', '185']:
            xmatch_freq(freq)
            sys.exit()
    print "usage: python correct_astrometry.py xmatch freq"
    # elif 'correct' in sys.argv:
    #     correct_all('refcat.fits')
    # elif 'plot' in sys.argv:
    #     make_plots()
    # else:
    #     print "xmatching: {0}+{2}->{1}".format(*sys.argv[-3:])
    #     ref_table = Table.read(sys.argv[-1])
    #     ref_cat = SkyCoord(ref_table['RAJ2000'], ref_table['DEJ2000'], unit=(u.degree, u.degree))
    #     xmatch_one(sys.argv[-3], sys.argv[-2], ref_cat, ref_table)