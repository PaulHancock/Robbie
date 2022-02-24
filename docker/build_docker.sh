#! /usr/bin/env bash

version=$(grep "^version" ../main.nf | awk '{print $3}' | sed s:\"::g )

# build the container using
docker build . -t "paulhancock/robbie-next:${version}"
# tag this build as latest
docker tag paulhancock/robbie-next:${version} paulhancock/robbie-next:latest