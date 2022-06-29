#! /usr/bin/env bash

# Run robbie viewer
# If argument exists, map the volume from this path to the docker data path
# Regexp for checking number
re='^[+-]?[0-9]+([.][0-9]+)?$'

# Check if only single argument exists and its not a number
if [ $1 ] && [ $# == 1 ] && ! [[ $1 =~ $re ]]
    then
        echo "Only specifying results path"
        # Convert path to absolute path
        abs_path=$(readlink -f $1)
        docker run -it -v $abs_path:/data -p 5006:5006 -it robbie/robbie-viewer:latest

# Check if first argument is string and also provided number args
elif [ $1 ] && ! [[ $1 =~ $re ]] && [ $# -gt 1 ]
    then
        echo "Specifying results path and RA, DEC and padding for cutout."
        # Convert path to absolute path
        abs_path=$(readlink -f $1)
        docker run --env RA=$2 --env DEC=$3 --env CUT=$4 -it -v $abs_path:/data -p 5006:5006 -it robbie/robbie-viewer:latest


# Check if first argument exists and it's a number
elif [ $1 ] && [[ $1 =~ $re ]]
    then
        echo "Specifying only RA, DEC and padding for cutout."
        # Convert path to absolute path
        abs_path=$(readlink -f $1)
        docker run --env RA=$1 --env DEC=$2 --env CUT=$3 -it -v ${PWD}/results:/data -p 5006:5006 -it robbie/robbie-viewer:latest

# If no arguments are given
elif [ $# -lt 1 ]
    then
        echo "Specifying default results path and no cutout."
        docker run -it -v ${PWD}/results:/data -p 5006:5006 -it robbie/robbie-viewer:latest
fi
