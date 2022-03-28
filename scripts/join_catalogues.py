#! /usr/bin/env python3

import astropy
from astropy.table import Table, Column
from astropy.io.votable import from_table, writeto
import numpy as np
import glob
import argparse
import sys

__author__ = "Paul Hancock"
__date__ = "2019-11-21"

def get_epoch_catalogues(epochs_file):
    """
    Read a file which contains a list of the catalogues to be read

    parameters
    ----------
    epochs_files : str
        A file which has a list of catalogues, one per line.

    returns
    -------
    files : list
        A list of filenames
    """
    files = list(map(str.strip, open(epochs_file).readlines()))
    return files


def join_catalogues(reference, epochs):
    """
    Take a reference cataloge, strip all but the uuid/ra/dec columns and then
    join the flux/err_flux data from each of the epoch catalogues
    From each epoch we extract the peak_flux and err_peak_flux columns and rename
    them by appending _N where N is the epoch number

    parameters
    ----------
    reference : str
        Filename for the reference catalogue

    epochs : list
        A list of the individual epoch file names

    returns
    -------
    table : `astropy.table.Table`
        A joined table
    """
    # Read the reference catalogue and retain only the uuid and ra/dec columns
    # rename the ra/dec columns
    print("Using reference catalogue {0}".format(reference))
    ref = Table.read(reference)['uuid', 'ra','dec']
    ref.rename_column('ra', 'ref_ra')
    ref.rename_column('dec', 'ref_dec')
    ref.sort(keys='uuid')

    # make the empty columns
    new_cols =[]
    data = np.zeros(len(ref), dtype=np.float32)*np.nan
    str_data = np.full(len(ref),Table.read(files[0])['image'][0])
    str_data[:] = 'None'

    for i in range(len(files)):
        for colname in ['peak_flux_{0}','err_peak_flux_{0}', 'local_rms_{0}', 'background_{0}']:
            new_cols.append(Column(data=data.copy(), name=colname.format(i)))
        for colname in ['image_{0}', 'epoch_{0}']:
            new_cols.append(Column(data=str_data.copy(), name=colname.format(i)))
    print("ref table is {0} rows".format(len(ref)))
    
    # if we add all the columns at once the ordering is borked!
    for col in new_cols:
        ref.add_column(col)
    
    # now put the data into the new big table
    for i,f in enumerate(files):
        print("Joining epoch {0} catalogue {1}".format(i,f))
        new_cols = Table.read(f)['uuid', 'peak_flux', 'err_peak_flux', 'local_rms', 'background',
                                 'image', 'epoch']
        new_cols.sort(keys='uuid')
        # compute the order/presence
        ordering = np.argwhere(np.in1d(ref['uuid'], new_cols['uuid'], assume_unique=True))[:,0]
        ref['peak_flux_{0}'.format(i)][ordering] = new_cols['peak_flux']
        ref['err_peak_flux_{0}'.format(i)][ordering] = new_cols['err_peak_flux']
        ref['local_rms_{0}'.format(i)][ordering] = new_cols['local_rms']
        ref['background_{0}'.format(i)][ordering] = new_cols['background']
        ref['image_{0}'.format(i)][ordering] = new_cols['image']
        ref['epoch_{0}'.format(i)][ordering] = new_cols['epoch']

    return ref


def write_table(tab, filename):
    """
    Write a VOTable using a binary format.

    parameters
    ----------
    tab : `astropy.table.Table`
        The table to be written

    filename : str
        The output filename. Should end with .xml or .vot
    """
    vot = from_table(tab)
    vot.set_all_tables_format('binary')
    vot.to_xml(filename)
    print("Wrote {0}".format(filename))
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Combine catalogues into a flux table")
    group1.add_argument("--refcat", dest='ref', type=str, default=None,
                        help="The reference catalogue")
    group1.add_argument("--epochs", dest='epochs', type=str, default=None,
                        help="A file containing the list of epoch catalogues.")
    group1.add_argument("--out", dest='outfile', type=str, default='flux_table.vot',
                        help="The output table name. Default = flux_table.vot")
    results = parser.parse_args()

    if None in (results.ref, results.epochs, results.outfile):
        parser.print_help()
        sys.exit()

    files = get_epoch_catalogues(results.epochs)
    table = join_catalogues(results.ref, files)
    write_table(table, results.outfile)
