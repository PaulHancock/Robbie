#! /usr/bin/env bash

#
# Make all the priorized catalogues for the 154 and 185 MHz images.
#

for freq in 154 185
do
    files=($(cat K2_${freq}MHz.dat | awk '{print $1}' ))
    for f in ${files[@]}
    do
	base=${f%%.fits}
	image=${base}_warped.fits
	table=${base}_warped_comp.fits
	bkg=${base}_bkg.fits
	rms=${base}_rms.fits

	if [[ ! -e ${table} ]]
	then
	    aegean ${image} --background ${bkg} --noise ${rms}\
                   --table ${image} --priorized 1 \
		   --input median_${freq}MHz_comp.fits --noregroup
	fi
    done
done