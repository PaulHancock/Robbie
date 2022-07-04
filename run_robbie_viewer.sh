#! /usr/bin/env bash

# Echo usage if something isn't right.
usage() { 
    echo "
    Usage: $0 [-c <ra,dec,pad>] [-p <path_to_results> ]

    where:
        -c centre position of image in RA and DEC with a padding value in degrees around the centre. For example: -c 330,-20,5
        -p path to the results directory
    " 1>&2; exit 1; 
}

while getopts ":c:p:h" opt; do
case ${opt} in 
    c )
        IFS=',' array=($OPTARG)
        ;;
    p ) 
        FILEPATH=$OPTARG
        ;;
    h)
        usage
    ;;
    \?)
    echo "ERROR: Invalid option -$OPTARG"
    usage
    ;;

    esac
done

# If coordinates and filepath
if ! [ -z "${array}" ] && ! [ -z "${FILEPATH}" ]; then
    echo "Specifying results path and RA, DEC and padding for cutout."
    # Convert path to absolute path
    abs_path=$(readlink -f ${FILEPATH})
    docker run -it -v $abs_path:/data -p 5006:5006 robbie/robbie-viewer:latest bokeh serve . --args /data ${array[0]} ${array[1]} ${array[2]}

# If only coordinates
elif  ! [ -z "${array}" ] && [ -z "${FILEPATH}" ]; then
    echo "Specifying only RA, DEC and padding for cutout."  
    docker run -it -v ${PWD}/results:/data -p 5006:5006 robbie/robbie-viewer:latest bokeh serve . --args /data ${array[0]} ${array[1]} ${array[2]}

# If only filepath
elif   [ -z "${array}" ] && ! [ -z "${FILEPATH}" ]; then
    # Convert path to absolute path
    abs_path=$(readlink -f ${FILEPATH})
    docker run -it -v $abs_path:/data -p 5006:5006 robbie/robbie-viewer:latest bokeh serve . --args /data

# If neither
elif  [ -z "${array}" ] && [ -z "${FILEPATH}" ]; then
    echo "Using defaults"
    docker run -it -v ${PWD}/results:/data -p 5006:5006 robbie/robbie-viewer:latest bokeh serve . --args /data 
fi
