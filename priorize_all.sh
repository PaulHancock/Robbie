#! /usr/bin/env bash

lowfreq=$(cat K2_154MHz.dat | awk '{print $1}')

for f in lowfreq
do
    echo aegean ${f} --autoload --table ${f} --priorized 2 --input median_154MHz.fits
done

highfreq=$(cat K2_186MHz.dat | awk '{print $1}')

for f in highfreq
do
    echo aegean ${f} --autoload --table ${f} --priorized 2 --input median_154MHz.fits
done