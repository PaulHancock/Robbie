#! /usr/bin/env bash

lowfreq=$(cat K2_154MHz.dat | awk '{print $1}')

if [[ -e median_154MHz_comp.fits ]]
then 
    for f in ${lowfreq}
    do
	echo aegean ${f} --autoload --table ${f} --priorized 2 --input median_154MHz_comp.fits
    done
fi

highfreq=$(cat K2_185MHz.dat | awk '{print $1}')

if [[ -e median_185MHz_comp.fits ]]
then
    for f in ${highfreq}
    do
	echo aegean ${f} --autoload --table ${f} --priorized 2 --input median_185MHz_comp.fits
    done
fi