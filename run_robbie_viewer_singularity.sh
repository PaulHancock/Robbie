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
        #array=$OPTARG
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

if [ -z "${FILEPATH}" ]; then
    abs_path=${PWD}/results/
else
    abs_path=$(readlink -f ${FILEPATH})
fi

if [ -z "${array}" ]; then
    pos_command=""
else
    pos_command=(${array[@]})
fi

# docker run -it -v $abs_path:/data -p 5006:5006 robbie/robbie-viewer:latest bokeh serve . --args /data ${pos_command[@]}
singularity exec --bind $abs_path:/data robbie_viewer.sif bokeh serve . --args /data ${pos_command[@]}
