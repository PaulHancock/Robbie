#! /bin/bash

field=$1

nimg=(`ls HGL_${field}/1*warped.fits | wc`)
nimg=${nimg[0]}

echo ${field} ${nimg}

sed s/NNN/${field}/g priorize_all.template | sed s/NIMG/${nimg}/g > priorize_${field}.sh
echo "sbatch priorize_${field}.sh"
