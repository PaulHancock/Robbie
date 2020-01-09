#! /usr/bin/env python
from __future__ import print_function

from astropy.table import Table
from astropy.io import fits
import numpy as np
import sqlite3
import argparse
import sys
import os


def get_date(filename):
    """
    Get date from the file

    parameters
    ----------
    filename : str
        Filename of fits image

    returns
    -------
    date : str
       String repr of the date (from DATE-OBS or similar)
    """
    header = fits.getheader(filename)
    date = header.get('DATE-OBS', 'NULL')
    return date


def add_catalogue(cur, date, cat):
    """
    Add the given catalogue to the database.

    parameters
    ----------
    cur : sqlite3.connection.cursor
        Cursor for accessing the database

    date: str
        Date for epoch.

    cat : str
        Catalogue file.
    """
    # open the catalogue
    # add the date to the epoch table and get the epoch number
    cur.execute("INSERT INTO epochs (date, file) VALUES (?,?)", (date,cat))
    cur.execute("SELECT epoch FROM epochs WHERE date=? and file=?", (date,cat))
    epoch = cur.fetchone()[0]
    print("File {0} has date {1} = epoch {2}".format(cat, date, epoch))
    c = Table.read(cat)
    # map catalogue columns into db columns, append epoch number
    cols = c.colnames
    # construct a query without having to hard code everything
    qry = "INSERT INTO sources ({0}, epoch) VALUES ({1},{2})".format(','.join(cols), ','.join(['?']*len(cols)), epoch)
    for row in c.as_array():
        # convert numpy types into python types
        data = list(map(lambda x: x.item(), row))
        #convert np.nan (in err fields) into a -1
        data = [ -1.0 if x is np.nan else x for x in data]
        cur.execute(qry, data)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("")
    group1.add_argument("--name", dest='name', type=str, default=None,
                        help="Database name")
    group1.add_argument("--cat",dest="cat", type=str, default=None,
                        help='Catalogue filename')
    group1.add_argument("--image", dest="image", type=str, default=None,
                        help='Image from which the catalogue was made.')
    results = parser.parse_args()

    if None in (results.name, results.cat, results.image):
        parser.print_help()
        sys.exit(1)

    conn = sqlite3.connect(results.name)
    c = conn.cursor()
    date = get_date(results.image)
    add_catalogue(c, date, results.cat)
    conn.commit()
    conn.close()

