#! /usr/bin/env python
from __future__ import print_function

from astropy.table import Table
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
    header = fits.gethead(filename)
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
    c = table.read(cat)
    # add the date to the epoch table and get the epoch number
    cur.execute("INSERT INTO epochs (date) VALUES (?)", date)
    cur.execute("SELECT epoch FROM epochs WHERE date=?", date)
    epoch = cur.fetchone()
    
    # map catalogue columns into db columns, append epoch number
    # input all rows
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

    conn = sqlite3.connect(results.name)
    c = conn.cursor()

    conn.commit()
    conn.close()

