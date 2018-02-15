#! /usr/bin/env bash

files=($( cat K2_1??MHz.dat | awk '{print $1}' ))
xmfiles=($( cat K2_1??MHz.dat | awk '{print $1}' | sed 's/.fits/_xm.fits/g' ))

rsync --ignore-existing --progress ${files[@]} galaxy:/astro/mwasci/phancock/warping/.
rsync --ignore-existing --progress ${xmfiles[@]} galaxy:/astro/mwasci/phancock/warping/.
