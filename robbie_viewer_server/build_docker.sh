#! /usr/bin/env bash

# build the container and add a target
docker build . -f ./Dockerfile -t "robbie/robbie-viewer"