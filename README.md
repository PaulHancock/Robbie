# Robbie: A batch processing workflow for the detection of radio transients and variables

[![Documentation Status](https://readthedocs.org/projects/robbie/badge/?version=latest)](https://robbie.readthedocs.io/en/latest/?badge=latest)
## Description

Robbie automates the process of cataloguing sources, finding variables, and identifying transients.

The workflow is described in [Hancock et al. 2018](https://ui.adsabs.harvard.edu/abs/2019A%26C....27...23H/abstract) and carries out the following steps:
- Preprocessing:
  - Find sources in images
  - Compare these catalogues to a reference catalogue
  - Use the offsets to model image based distortions
  - Make warped/corrected images
- Persistent source catalogue creation:
  - Stack the warped images into a cube and form a mean image
  - Source find on the mean image to make a master catalogue
  - Priorized fit this catalogue into each of the individual images
  - Join the catalogues into a single table and calculate variability stats
- Transient candidate identification:
  - Use the persistent source to mask known sources from the individual images
  - Source find on the masked images to look for transients
  - Combine transients tables into a single catalogue, identifying the epoch of each detection

## Dependencies
Robbie relies on the following software:
- [AegeanTools](https://github.com/PaulHancock/Aegean)
- [fits_warp](https://github.com/nhurleywalker/fits_warp)
- [Stils/TOPCAT](http://www.star.bris.ac.uk/~mbt/topcat/)
- [Nextflow](https://www.nextflow.io/)
- [SWarp](https://www.astromatic.net/software/swarp/)

All dependencies except for Nextflow will be installed in the docker image.

## Installation
The best way to use Robbie is via a docker container that has all the software dependencies installed. Ensure docker is running, then build the container using:
```
docker build -t paulhancock/robbie-next -f docker/Dockerfile .
```

or by pulling the latest build from [DockerHub](https://hub.docker.com/r/paulhancock/robbie-next) via
```
docker pull paulhancock/robbie-next
```

Then, install Nextflow with a package management system such as Conda:

```
conda install -c bioconda nextflow
```

Once Nextflow is installed, add robbie.nf to your path with
```
python setup.py install
```

## Quickstart
Robbie now uses Nextflow to manage the workflow and can be run on a local system or a supercomputing cluster. You can use a container via singularity, docker, or the host's software. The current development cycle tests Robbie using singularity on an HPC with the Slurm executor - other setups *should* work but haven't been extensively tested.

### `images.txt`
Before running Robbie, you will need to create a text file that contains the paths to each image to be processed. For example, if within the "Robbie" parent directory there exists a folder named "images" containing the `.fits` files:

```
ls images/* > images.txt
```

will populate `images.txt` with the image paths relative to the parent directory.

### `robbie.nf`
This file describes the workflow and can be inspected but shouldn't be edited directly. To describe the command line arguments, use
```
robbie.nf --help
```

### `nextflow.config`
This file is the configuration setup and contains all the command line arguments' default values. You can change these defaults by copying the `nextflow.config` and editing the relevant params.\<argument\>. You can then use your custom config via:
```
nextflow -C my.config run robbie.nf
```
The `-C my.config` directs Nextflow to use *only* the configuration described in `my.config`. If you use `-c`, then it will also read the `nextflow.config` file.

### `-profile`

If you're running Robbie on your local machine, you should use the `-profile local` option to use the Robbie docker image. For example:

```
nextflow -C my.config run robbie.nf -profile local
```

If you're running Robbie on a supercomputing cluster (HPC), you should use the relevant cluster profile (`-profile zeus` or `-profile magnus`) to assure you're using the cluster's job queue (such as Slurm). If there isn't a profile for your cluster (check in `nextflow.config`), you may have to make your own.

Additional configuration files are stored in the `./config` directory and may be useful templates for your work.

## Visualisation

Firstly, build the Docker image located in the robbie_viewer_server directory:

```
./build_docker.sh
```

once this is complete, run the viewer in the main Nextflow directory via:

```
./run_robbie_viewer.sh
```

This will run the viewer using the images within the default ``results`` directory. If your directory is different to the default, you can add either the relative or absolute path as an optional argument:

```
./run_robbie_viewer.sh -p path_to_dir
```

## Credit
If you use Robbie as part of your work, please cite [Hancock et al. 2018](http://adsabs.harvard.edu/abs/2019A%26C....27...23H), and link to this repository.

## Links
You can obtain a docker image with the Robbie dependencies installed at [DockerHub](https://hub.docker.com/r/paulhancock/robbie-next/)
