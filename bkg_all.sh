#! /usr/bin/env bash

files=$(ls K2_trim_[0-9]*[0-9].fits)
for f in ${files}
do
  bane=${f%%.fits}_bkg.fits
  make ${bane}
done