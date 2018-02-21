#!/usr/bin/env bash

for f in $( cat bad_files.dat )
do
    rm ${f}
done
