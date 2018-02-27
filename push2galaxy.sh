#! /usr/bin/env bash

files=($( cat K2_1??MHz.dat | awk '{print $1}' ))
bkg=($( cat K2_1??MHz.dat | awk '{print $1}' | sed 's/.fits/_bkg.fits/g'  ))
rms=($( cat K2_1??MHz.dat | awk '{print $1}' | sed 's/.fits/_rms.fits/g' ))
xmfiles=($( cat K2_1??MHz.dat | awk '{print $1}' | sed 's/.fits/_xm.fits/g' ))

rsync --ignore-existing --progress ${files[@]} ${xmfiles[@]} ${bkg[@]} ${rms[@]} galaxy:/astro/mwasci/phancock/warping/.
rsync --ignore-existing --progress K2_1??MHz.dat galaxy:/astro/mwasci/phancock/warping/.
rsync --ignore-existing --progress median_{154,185}MHz{,_bkg,_rms,_comp}.fits galaxy:/astro/mwasci/phancock/warping/.
rsync --progress galaxy_codes/* galaxy:/astro/mwasci/phancock/warping/.