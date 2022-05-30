#! /usr/bin/env bash

# run robbie viewer
docker run -v ${PWD}/results:/data -p 5006:5006 robbie/robbie-viewer:latest