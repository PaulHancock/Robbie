#! /bin/bash

field="K2"

nimg=$( ls K2_trim_[0-9]*[0-9]_warped.fits | wc -l )

echo ${field} ${nimg}

sed s/NNN/${field}/g priorize_all.template | sed s/NIMG/${nimg}/g > priorize_${field}.sh
echo $( sbatch priorize_${field}.sh )
