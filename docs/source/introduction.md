# Description

Robbie automates the process of cataloguing sources, finding variables, and identifying transients.

The workflow was described initially in [Hancock et al. 2018](https://ui.adsabs.harvard.edu/abs/2019A%26C....27...23H/abstract) however the current workflow is shown below:
- Preprocessing:
  - Convolve all images to a common psf (optional)
  - Create background and noise maps (if they are not found)
  - Correct astrometry using fitswarp (optional)
- Variabile/persistent source detection:
  - Stack the warped images to form a mean image
  - Source find on the mean image to make a reference catalogue
  - Priorized fit this catalogue into each of the individual images
  - Join the epoch catalogues to make a persistent source catalogue
  - Calculate variability stats and generate a light curve for each source
- Transient candidate identification:
  - Use the persistent source to mask known sources from the individual epochs
  - Source find on the masked images to find transients
  - Concatenate transients into a single catalogue, identifying the epoch of each detection

See [workflow](workflow) for a diagram of how Robbie works.

## Dependencies
Robbie relies on the following software:
- [AegeanTools](https://github.com/PaulHancock/Aegean)
- [fits_warp](https://github.com/nhurleywalker/fits_warp)
- [Stils/TOPCAT](http://www.star.bris.ac.uk/~mbt/topcat/)
- [Nextflow](https://www.nextflow.io/)
- [SWarp](https://www.astromatic.net/software/swarp/)
- [Astropy](https://www.astropy.org/)
- [reproject](https://reproject.readthedocs.io/en/stable/)

All dependencies except for Nextflow will be installed in the docker image.

## Credit
If you use Robbie as part of your work, please cite [Hancock et al. 2018](http://adsabs.harvard.edu/abs/2019A%26C....27...23H), and link to this repository.

## Links
You can obtain a docker image with the Robbie dependencies installed at [DockerHub](https://hub.docker.com/r/paulhancock/robbie-next/)
