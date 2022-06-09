## Installation

### Robbie workflow

The best way to use Robbie is via a docker container that has all the software dependencies installed. Ensure docker is running, then build the container using:

``` bash
docker build -t paulhancock/robbie-next -f docker/Dockerfile .
```

or by pulling the latest build from [DockerHub](https://hub.docker.com/r/paulhancock/robbie-next) via

``` bash
docker pull paulhancock/robbie-next
```

Then, install Nextflow with a package management system such as Conda:

``` bash
conda install -c bioconda nextflow
```

Once Nextflow is installed, add robbie.nf to your path with

``` bash
python setup.py install
```

Please see [Quickstart](quickstart) on how to run Robbie once the setup is complete.

### Robbie visualisation

To construct the Docker image for visualisation of the Robbie results, run the following:

``` bash
./build_docker.sh
``` 

within the ``robbie_viewer_server`` directory. Please see [Visualisation](visualisation) on how to visualise the results of Robbie using the aforementioned Docker image.