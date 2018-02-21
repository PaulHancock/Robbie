#!/bin/bash -l

stilts="java -jar /home/hancock/Software/stilts.jar"

for freq in 154 185
do
    files=($(cat K2_${freq}MHz.dat | awk '{print $1}' | sed 's/.fits/_warped_comp.fits/g' | xargs ls ))
    cmd="${stilts} tmatchn nin=${#files[@]} matcher=exact"
    for i in $( seq 0 1 ${#files[@]} )
    do
	if [[ -e ${files[${i}]} ]]
	then
	    j=$( echo "${i} +1" | bc )
	    cmd="${cmd} in${j}=${files[${i}]} suffix${j}=_${j} values${j}='uuid'"
	fi
    done
    cmd="${cmd} out=joined_${freq}MHz.csv ofmt=csv"
    $( ${cmd} )
done