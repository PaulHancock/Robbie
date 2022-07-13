## Installation

Robbie relies on 3 core technologies to run:
- nextflow
- python
- docker or singularity containers (optional)

### Nextflow

Nextflow can be installed in a few ways:
- By following instructions at [nextflow.io](https://www.nextflow.io/docs/latest/getstarted.html)
  - `wget -qO- https://get.nextflow.io | bash`
- By using a packaged manager such as [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html):
  - `conda install -c bioconda nextflow`
- If you are using Pawsey you can just run:
  - `module load nextflow`

### Robbie scripts
The Robbie scripts are all bundled as a python package which can be downloaded and installed via github.

Using pip: `pip install https://github.com/ADACS-Australia/Robbie.git`

Using Conda: `conda install git pip` (and then run the above)

### Docker containers

The best way to use Robbie is via a docker container that has all the software dependencies installed. Ensure docker is running, then build the container using:

``` bash
docker build -t paulhancock/robbie-next -f docker/Dockerfile .
```

or by pulling the latest build from [DockerHub](https://hub.docker.com/r/paulhancock/robbie-next) via

``` bash
docker pull paulhancock/robbie-next
```

Please see [Quickstart](quickstart) on how to run Robbie once the setup is complete.

### Robbie visualisation

To construct the Docker image for visualisation of the Robbie results, run the following:

``` bash
./build_docker.sh
``` 

within the ``robbie_viewer_server`` directory. Please see [Visualisation](visualisation) on how to visualise the results of Robbie using the aforementioned Docker image.