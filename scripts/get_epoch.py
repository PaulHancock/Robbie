#! /usr/bin/env python3

from astropy.io import fits
import argparse


def get_obs_date(fname):
    """
    Get the observation date from the fits header of a given file
    """
    return fits.getheader(fname)['DATE-OBS']


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("file", type=str,
                        help="Filename (fits image)")
    results = parser.parse_args()


    print(get_obs_date(results.file))
