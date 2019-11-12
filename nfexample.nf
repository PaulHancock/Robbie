#! /usr/bin/env nextflow

params.input_dir = 'data/'
params.ref_catalogue = "reference.fits"
params.region_file = 'square.mim'
params.output_dir = 'results/'

ref_catalogue_ch = Channel.value("${params.ref_catalogue}")
region_ch = Channel.value("${params.region_file}")

raw_image_ch = Channel.fromPath("${params.input_dir}/Epoch0?.fits").map{ it -> [it.baseName, it]}


process bane_raw {
  publishDir params.output_dir, mode:'copy', overwrite:true

  input:
  set val(basename), file(image) from raw_image_ch

  output:
  set val(basename), file(image), file("${basename}_bkg.fits"), \
     file("${basename}_rms.fits") into raw_image_with_bkg_ch

  script:
  """
  echo input ${image}
  touch ${basename}_bkg.fits
  touch ${basename}_rms.fits
  """
}

// raw_image_with_bkg_ch.view()

process initial_sfind {

  publishDir params.output_dir, mode:'copy', overwrite:true

  input:
  set val(basename), file(image), file(bkg), file(rms) from raw_image_with_bkg_ch

  output:
  set val(basename), file(image), file("${basename}_comp.fits") into initial_catalogue_ch

  script:
  """
  echo aegean --bkg ${bkg} --rms ${rms} --table ${image} ${image}
  touch ${basename}_comp.fits
  """
}

// initial_catalogue_ch.view()

process fits_warp {
  publishDir params.output_dir, mode:'copy', overwrite:true
  input:
  set val(basename), file(image), file(catalogue) from initial_catalogue_ch
  file(ref) from ref_catalogue_ch

  output:
  set val("${basename}_warped}"), file("${basename}_warped.fits") into warped_images_ch

  script:
  """
  echo fits_warp on ${image} with ${catalogue} and ${ref}
  touch ${basename}_warped.fits
  """
}

