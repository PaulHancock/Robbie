#! /usr/bin/env bash

files=$(ls K2_trim_*.fits)
for f in ${files}
do
  bane=${f%%.fits}_bkg.fits
  make ${bane}
done