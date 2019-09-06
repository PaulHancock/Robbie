#! /usr/bin/env python

from __future__ import print_function
__author__ = 'Paul Hancock'
__date__ = '2019/09/06'

import argparse
from astropy import table
import numpy as np
import os
import sys


def add_epoch(table, epoch):
    """
    Add an epoch column to the table.

    :param table:
    :param epoch:
    :return:
    """
    col = table.Column(data=np.ones(len(table), dtype=np.int32)*epoch, name='Epoch')
    table.add_column(col)
    return table


def join_all(flist):
    """
    Concatenate all the tables and add an extra column which identifies the epoch for each row

    parameters
    ==========
    flist : [str, ...]
        A list of the catalogue file names to be read

    return
    ======
    tab : astropy.table.Table
        A joined table
    """
    print("Reading {0}".format(flist[0]))
    tab = table.Table.read(flist[0])
    tab = add_epoch(tab, 0)

    if len(flist) <= 1:
        return tab

    for i in range(1, len(flist)):
        print("Appending {0}".format(flist[i]))
        tab2 = add_epoch(table.Table.read(flist[i]), i)
        tab = table.vstack([tab, tab2])

    return tab


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Collect catalogues")
    group1.add_argument("--infile", dest='infile', type=str, default=None,
                        help="A list of catalogues in a file. [optional]")
    group1.add_argument("--in", dest='files', type=str, default=None, nargs='+',
                        help="Explicit list of catalogues to include. [optional]")
    group1.add_argument("--out", dest='outfile', type=str, default=None,
                        help="Output filename.")
    results = parser.parse_args()

    if results.outfile is None:
        parser.print_help()
        sys.exit(1)

    flist = []
    if results.infile:
        flist.extend([a.strip() for a in open(results.infile).readlines()])
    if results.files:
        flist.extend(results.files)

    tab = join_all(flist)
    if os.path.exists(results.outfile):
        os.remove(results.outfile)
    tab.write(results.outfile)
