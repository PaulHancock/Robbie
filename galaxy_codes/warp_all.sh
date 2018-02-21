#! /bin/bash

field="K2"

nimg=($( ls K2_trim_[0-9]*[0-9].fits | wc -l ))

echo ${field} ${nimg}

sed s/NNN/${field}/g warp_all.template | sed s/NIMG/${nimg}/g > warp_${field}.sh
echo $( sbatch warp_${field}.sh )
