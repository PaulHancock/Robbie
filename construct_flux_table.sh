#! /usr/bin/env bash

for freq in 154 185
do
    files=( $( cat K2_${freq}MHz.dat | awk '{print $1}' | sed 's/.fits/_comp.fits/g' | xargs ls ) )
    cmd="java -jar /home/hancock/Software/stilts.jar tmatchn nin=${#files[@]} matcher=exact out=${freq}MHz_flux_table.fits"
    for n in ${!files[@]}
    do
	m=$( echo "${n}+1" | bc )
	time=$( echo ${files[${n}]} | sed 's/.*_\([0-9]*\)_.*/\1/g' )
	cmd="${cmd} in${m}=${files[${n}]} values${m}='uuid' suffix${m}=_${time}"
    done
    `${cmd}`
done