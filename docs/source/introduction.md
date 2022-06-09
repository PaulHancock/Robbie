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
- [Astropy](https://www.astropy.org/)
- [reproject](https://reproject.readthedocs.io/en/stable/)

All dependencies except for Nextflow will be installed in the docker image.

## Credit
If you use Robbie as part of your work, please cite [Hancock et al. 2018](http://adsabs.harvard.edu/abs/2019A%26C....27...23H), and link to this repository.

## Links
You can obtain a docker image with the Robbie dependencies installed at [DockerHub](https://hub.docker.com/r/paulhancock/robbie-next/)
