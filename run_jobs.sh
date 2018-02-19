#! /usr/bin/env bash
ssh galaxy 'bash -s' << 'ENDSSH'
cd /astro/mwasci/phancock/warping
./warp_all.sh
ENDSSH
