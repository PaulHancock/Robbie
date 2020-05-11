#! /bin/bash -l
#SBATCH --export=NONE
#SBATCH --time=2:00:00
#SBATCH --nodes=1

set -ex
#ensure data is owned by mwasci
newgrp mwasci

module load singularity 
# see https://stackoverflow.com/a/39454426 for the following hack
latest=$( wget -q https://registry.hub.docker.com/v1/repositories/paulhancock/robbie-next/tags -O -  | sed -e 's/[][]//g' -e 's/"//g' -e 's/ //g' | tr '}' '\n'  | awk -F: '{print $3}' | sort -V | tail -n1 )

cd /group/mwasci/phancock/.singularity
singularity pull docker://paulhancock/robbie-next:${latest}
rm paulhancock-robbie-next.img
ln -s robbie-next_${latest}.sif paulhancock-robbie-next.img
cd ../
# update the permissions
chmod ugo+rx .singularity
