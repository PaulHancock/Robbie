# Robbie: A batch processing work-flow for the detection of radio transients and variables

## Description

Robbie automates the process of cataloguing sources, finding variables, and identifying transients.

The workflow is described in Hancock et al. 2018 (submitted), and carries out the following steps:
- Find sources in images
- Compare these catalogues to a reference catalogue
- Use the offsets to model image based distortions
- Make warped/corrected images
- Stack the warped images into a cube and form a mean image
- Source find on the mean image to make a master catalogue
- Priorized fit this catalogue into each of the individual images
- Join the catalogues into a single table and calculate variability stats
- Use the master catalogue to mask known sources from the individual images
- Source find on the masked images to look for transients
- Combine transients tables into a single catalogue, identifying the epoch of each detection

## Configuration
You need to have the following software installed in order to use Robbie:
- [AegeanTools](https://github.com/PaulHancock/Aegean)
- [fits_warp](https://github.com/nhurleywalker/fits_warp)
- [Stils/TOPCAT](http://www.star.bris.ac.uk/~mbt/topcat/)

The included `Makefile` should be edited to set up some custom parameters.
In particular you need to set:
- `STILTS` = `<however you would run stilts from the command line>`
- `IMAGEFILE` = a file that contains a list of all the images in epoch order (default=all_images.txt)
- `REFCAT` = /path/to/your/external/reference/catalogue.fits
- `REGION` = a [MIMAS](https://github.com/PaulHancock/Aegean/wiki/MIMAS) region file describing the region of interest.

## Usage
```
Usage is: make [command | file]
 files:
  refcat.fits - a masked version of the external reference catalogue
  cube.fits - a stack of astrometry corrected images
  mean.fits - a mean image from the above stack
  flux_table_var.fits - light curves and variability stats for all persistent sources
  transients.fits - a catalogue of all candidate transient events
  transients.png - a visualisation of transients.fits

 commands:
  transients = transients.png
  variables = flux_table_var.fits
  sceince = variables + transients

I recommend that you `make science`
```

## Credit
If you make use of Robbie as part of your work please cite [Hancock et al. 2018](http://adsabs.harvard.edu/abs/2019A%26C....27...23H), and link to this repository.

## Links
You can obtain a docker image with the Robbie dependencies installed at [DockerHub](https://hub.docker.com/repository/docker/paulhancock/robbie-next/)