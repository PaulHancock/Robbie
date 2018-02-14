#! /usr/bin/env bash

files=($( cat K2_1??MHz.dat | awk '{print $1}' ))
xmfiles=($( cat K2_1??MHz.dat | awk '{print $1}' | sed 's/.fits/_xm.fits/g' ))

rsync --ignore-existing --progress ${files[@]1::4} galaxy:/astro/mwasci/phancock/warping/.
rsync --ignore-existing --progress ${xmfiles[@]1::4} galaxy:/astro/mwasci/phancock/warping/.
