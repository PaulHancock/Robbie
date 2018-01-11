#! /usr/bin/env bash

[[ -e freqs.dat ]] && rm freqs.dat

for f in $(ls K2_trim_[0-9]*[0-9].fits);
do
    echo ${f} `gethead ${f} CRVAL3` >> freqs.dat
done

grep ' 15' freqs.dat > K2_154MHz.dat
grep ' 18' freqs.dat > K2_185MHz.dat
