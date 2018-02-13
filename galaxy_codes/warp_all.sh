#! /bin/bash

field=$1

nimg=(`ls HGL_${field}/1*image.fits | wc`)
nimg=${nimg[0]}

echo ${field} ${nimg}

sed s/NNN/${field}/g warp_all.template | sed s/NIMG/${nimg}/g > warp_${field}.sh
echo "sbatch warp_${field}.sh"
