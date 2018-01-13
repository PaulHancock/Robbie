#! /usr/bin/env bash

#make clean
make input
./trim_all.sh
# ./bkg_all.sh
# make K2_154MHz.dat K2_185MHz.dat
# make cube_154MHz.fits cube_185MHz.fits
# make median_154MHz.fits median_185MHz.fits
make median_154MHz_comp.fits median_185MHz_comp.fits
./priorize_all.sh
