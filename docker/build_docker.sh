#! /usr/bin/env bash

version=$(grep "^version" ../main.nf | awk '{print $3}' | sed s:\'::g )

# build the container using
docker build . -t "robbie-next:${version}"
# tag this build as latest
docker tag robbie-next:${version} robbie-next:latest