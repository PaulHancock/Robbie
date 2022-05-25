#! /usr/bin/env bash

# if directory doesn't exist, copy from parent Robbie results
if [ ! -d "$results" ]; then
    cp -r ../results/ .
    cp ./results/*.vot* ./results/reprojected_images/
fi

# build the container and add a tage
docker build -f ./Dockerfile -t "robbie/robbie-viewer" .