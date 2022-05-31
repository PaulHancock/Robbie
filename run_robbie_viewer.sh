#! /usr/bin/env bash

# Run robbie viewer
# If argument exists, map the volume from this path to the docker data path
if [ $1 ]
then
    # Convert path to absolute path
    abs_path=$(readlink -f $1)
    docker run -v $abs_path:/data -p 5006:5006 robbie/robbie-viewer:latest
else
    docker run -v ${PWD}/results:/data -p 5006:5006 robbie/robbie-viewer:latest
fi
