#! /usr/bin/env bash

files=$(ls K2_final_*.fits)
for f in ${files}
do
  trim=$(echo ${f} | sed 's/final/trim/g')
  make ${trim}
done
