#!/usr/bin/env bash

for freq in 154 185
do
    files=($( cat K2_${freq}MHz.dat | awk '{print $1}' ))
	for f in ${files[@]}
	do
	    make ${f%%.fits}_warped_blanked.fits
	done
done