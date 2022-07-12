#! /usr/bin/env bash

# build the container and add a target
docker buildx build . --platform linux/amd64 -f ./Dockerfile -t "paulhancock/robbie-viewer:latest"
