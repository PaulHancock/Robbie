#! /usr/bin/env python
from __future__ import print_function

import sqlite3
import argparse
import sys
import os


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("")
    group1.add_argument("--name", dest='name', type=str, default=None,
                        help="Database name")
    results = parser.parse_args()

    if os.path.exists(results.name):
        os.remove(results.name)

    conn = sqlite3.connect(results.name)
    c = conn.cursor()
    c.execute("""
    CREATE TABLE epochs
    (date TEXT,
    epoch INTEGER PRIMARY KEY ASC AUTOINCREMENT)
    """)
    c.execute("""
    CREATE TABLE sources
    (island INTEGER, source INTEGER,
    background NUMERIC, local_rms NUMERIC,
    ra_str TEXT, dec_str TEXT,
    ra NUMERIC, err_ra NUMERIC, dec NUMERIC, err_dec NUMERIC,
    peak_flux NUMERIC, err_peak_flux NUMERIC,
    INTEGER_flux NUMERIC, err_INTEGER_flux NUMERIC,
    a NUMERIC, err_a NUMERIC, b NUMERIC, err_b NUMERIC,
    pa NUMERIC, err_pa NUMERIC,
    flags INTEGER,
    residual_mean NUMERIC, residual_std NUMERIC,
    uuid TEXT,
    psf_a NUMERIC, psf_b NUMERIC, psf_pa NUMERIC,
    epoch INTEGER REFERENCES epochs(epoch),
    UNIQUE (uuid, epoch))
    """)
    conn.commit()
    conn.close()

