#! /usr/bin/env bash

rsync --ignore-existing --progress galaxy:'/astro/mwasci/phancock/warping/*_warped.fits' .
rsync --ignore-existing --progress galaxy:'/astro/mwasci/phancock/warping/*_warped_comp.{fits,reg}' .