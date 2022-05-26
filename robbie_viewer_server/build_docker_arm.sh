#! /usr/bin/env bash

# if directory doesn't exist, copy from parent Robbie results
if [ ! -d "./results" ] 
then
    mkdir results
    cp -r ../results/reprojected_images/* ./results/
    cp ../results/*.vot* ./results/
else
    rm -r results
    mkdir results
    cp -r ../results/reprojected_images/* ./results/
    cp ../results/*.vot* ./results/
fi

# build the container and add a target
docker buildx build . --platform linux/amd64 -f ./Dockerfile -t "robbie/robbie-viewer"