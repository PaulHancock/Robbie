#! /usr/bin/env bash

# build the container and add a target
docker buildx build . --platform linux/amd64 -f ./Dockerfile -t "cjproud/robbie_viewer:latest"