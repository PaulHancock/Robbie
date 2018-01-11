#! /usr/bin/env bash

#
# Make all the priorized catalogues for the 154 and 185 MHz images.
#

for freq in 154 185
do
    files=$(cat K2_${freq}MHz.dat | awk '{print $1}')
    if [[ -e median_${freq}MHz_comp.fits ]]
    then
	for f in ${files}
	do
	    make ${f%%.fits}_bkg.fits
	    aegean ${f} --autoload --table ${f} --priorized 2 --input median_${freq}MHz_comp.fits
	done
    fi
done