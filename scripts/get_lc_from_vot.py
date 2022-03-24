#! /usr/bin/env python3

import astropy
from astropy.table import Table, Column
from astropy.io.votable import from_table, writeto
from astropy.io import ascii
import numpy as np
import glob
import argparse
import sys


def get_lc(tab, uuid):
    """
    Open a table of fluxes and extract a light curve for a single row.
    Convert this row into a new table with columns = [Date, peak_flux, err_peak_flux, epoch, local_rms]

    parameters
    ----------
    tab : str
        Filename of the flux table
    uuid : str
        The name of the source to be extracted

    returns
    -------
    lc : astropy.table
       A table of the light curve
    """
    # load the table and extract the row of interest
    master = Table.read(tab)
    row_id = np.where(np.array(list(map(str.strip,master['uuid']))) == uuid)
    row = master[row_id]
    if len(row) == 0:
        print("uuid {0} not found".format(uuid))
        return None
    elif len(row)>1:
        print("uuid {0} found but not unique".format(uuid))
        return None
    # count the number of epochs and extract the data in epoch order
    nepochs = max([int(n.split('_')[-1]) for n in master.colnames if n.startswith('epoch')])
    date = list([row[name].tolist()[0] for name in ['epoch_{0}'.format(i) for i in range(0, nepochs+1)]])
    peak_flux = list([row[name].tolist()[0] for name in ['peak_flux_{0}'.format(i) for i in range(0, nepochs+1)]])
    err_peak_flux = list([row[name].tolist()[0] for name in ['err_peak_flux_{0}'.format(i) for i in range(0, nepochs+1)]])
    epoch = list(range(0,nepochs+1))
    local_rms = list([row[name].tolist()[0] for name in ['local_rms_{0}'.format(i) for i in range(0, nepochs+1)]])
    # convert the lists into a table
    lc = Table(data = [date, peak_flux, err_peak_flux, epoch, local_rms], names=['Date', 'peak_flux', 'err_peak_flux', 'epoch','local_rms'])
    return lc

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Extract a lightcurve from the flux table")
    group1.add_argument("--table", dest='tab', type=str, default=None,
                        help="The flux table")
    group1.add_argument("--name", dest='uuid', type=str, default=None,
                        help="The name (uuid) of the source to be extracted")
    group1.add_argument("--out", dest='outfile', type=str, default='out.csv',
                        help="The output file name (csv format). Default = out.csv")
    results = parser.parse_args()

    if None in (results.tab, results.uuid, results.outfile):
        parser.print_help()
        sys.exit()

    lc = get_lc(tab=results.tab, uuid=results.uuid)
    if lc is None:
        sys.exit(1)
    ascii.write(lc, results.outfile, format='csv', overwrite=True)
    print("Extracted light curve: {0} -> {1}".format(results.tab, results.outfile))
