#! /usr/bin/env bash

#
# Make all the priorized catalogues for the 154 and 185 MHz images.
#
for freq in 154 185
do
    files=($( cat K2_${freq}MHz.dat | awk '{print $1}' ))
	for f in ${files[@]}
	do
	    make ${f%%.fits}_warped_blanked_comp_filtered.fits
	done
done