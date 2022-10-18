(visualisation)=
## Visualisation

### Running locally with Docker
To start the Docker container containing the Bokeh server, run the following script in the main Nextflow directory:

``` bash
./run_robbie_viewer.sh
```

This will run the viewer using the images output from Robbie within the default ``results`` directory. If your output directory is different to the default, you can add either the relative or absolute path as an optional argument:

``` bash
./run_robbie_viewer.sh -p path_to_dir
```

When plotting large images, it is recommended to also specify an RA and DEC position, as well as a size in coordinate units, to cutout a portion of the image for plotting. For example, if we want to plot an image with centre position of RA 335°, DEC -15° and size of 5°:

``` bash
./run_robbie_viewer.sh -p path_to_dir -c 335,-15,5
```

### Running on a cluster with Singularity

The Robbie Viewer is available on Pawsey as a part of the SHPC (Singularity Recipe HPC). To install it, we will load the SHPC module, install the viewer as a module and then load it:

``` bash
module load shpc/0.0.53
shpc install cjproud/robbie_viewer
module load cjproud/robbie_viewer/latest/module
```

Now, the viewer will be available on our path and we can run it as normal:

``` bash
./run_robbie_viewer.sh -p path_to_dir -c RA,DEC,PAD
```


### Visualising transients 

Visualising different transient candidates can be done in multiple ways. For example, the transient candidate can be selected using the table, sky plot or variables plot as shown below:

![Cycling through different transients](images/Transient_Selector_Example.gif)

### Visualising transient Epoch's

Cycling through Epoch's for each transient candidate is just as easy, for example you can use either the Epoch slider or select each Epoch in the Peak Flux vs. Epoch graph:

![Cycling through different Epochs](images/Epoch_Example.gif)

### Transient candidate selection

Bokeh has multiple ways to interact with the data shown in the plots and table. To select multiple transient candidates, one option is to hold ``shift`` and click on the table entries. Once we zoom out, we can see all the selected transients on each plot:

![Transient selection 1](images/Selection_Example1.gif)

The Box select tool can also be used to select transient candidates. After drawing the bounding box for selection, the transient candidates are highlighted in the other plots as well as the table below:

![Transient selection 2](images/Selection_Example2.gif)