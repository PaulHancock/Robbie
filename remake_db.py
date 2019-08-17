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
    (date text, epoch int)
    """)
    c.execute("""
    CREATE TABLE sources
    (island int, source int,
    background float, local_rms float,
    ra_str text, dec_str text,
    ra float, err_ra float, dec float, err_dec float,
    peak_flux float, err_peak_flux float,
    int_flux float, err_int_flux float,
    a float, err_a float, b float, err_b float,
    pa float, err_pa float,
    flags int,
    residual_mean float, residual_std float,
    uuid text,
    psf_a float, psf_b float, psf_pa float,
    epoch int REFERENCES epochs(epoch),
    UNIQUE (uuid, epoch))
    """)
    conn.commit()
    conn.close()

