#! /usr/bin/env bash
ssh galaxy 'bash -s' << 'ENDSSH'
cd /astro/mwasci/phancock/warping
./priorize_all.sh
ENDSSH
