(quickstart)=
## Quickstart
Robbie now uses Nextflow to manage the workflow and can be run on a local system or a supercomputing cluster. You can use a container via singularity, docker, or the host's software. The current development cycle tests Robbie using singularity on an HPC with the Slurm executor - other setups *should* work but haven't been extensively tested.

### `images.txt`
Before running Robbie, you will need to create a text file that contains the paths to each image to be processed. For example, if within the "Robbie" parent directory there exists a folder named "images" containing the `.fits` files:

``` bash
ls images/* > images.txt
```

will populate `images.txt` with the image paths relative to the parent directory.

### `robbie.nf`
This file describes the workflow and can be inspected but shouldn't be edited directly. To describe the command line arguments, use
``` bash
robbie.nf --help
```

### `nextflow.config`
This file is the configuration setup and contains all the command line arguments' default values. You can change these defaults by copying the `nextflow.config` and editing the relevant params.\<argument\>. You can then use your custom config via:
``` bash
nextflow -C my.config run robbie.nf
```
The `-C my.config` directs Nextflow to use *only* the configuration described in `my.config`. If you use `-c`, then it will also read the `nextflow.config` file.

### `-profile`

If you're running Robbie on your local machine, you should use the `-profile local` option to use the Robbie docker image. For example:

``` bash
nextflow -C my.config run robbie.nf -profile local
```

If you're running Robbie on a supercomputing cluster (HPC), you should use the relevant cluster profile (`-profile zeus` or `-profile magnus`) to assure you're using the cluster's job queue (such as Slurm). If there isn't a profile for your cluster (check in `nextflow.config`), you may have to make your own.

Additional configuration files are stored in the `./config` directory and may be useful templates for your work.
