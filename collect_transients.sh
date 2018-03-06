#!/usr/bin/env bash


stilts="java -jar /home/hancock/Software/stilts.jar"
for freq in 154 185
do
    files=($(cat K2_${freq}MHz.dat | awk '{print $1}' | sed 's/.fits/_warped_blanked_comp_filtered.fits/g' | xargs ls ))
    cmd="${stilts} tcatn nin=${#files[@]}"
    for i in $( seq 1 1 ${#files[@]} )
    do
	    j=$( echo "${i} -1" | bc )
	    time=$( echo ${files[${j}]} | sed 's/.*_\([0-9]*\)_.*/\1/g' )
	    cmd="${cmd} in${i}=${files[${j}]} icmd${i}='addcol epoch ${i}'"
    done
    cmd="${cmd} out=${freq}MHz_transients.fits ofmt=fits"
    echo ${cmd} | bash
done