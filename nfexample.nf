#! /usr/bin/env nextflow

params.input_dir = 'data/'
params.ref_catalogue = "reference.fits"
params.region_file = "$baseDir/square.mim"
params.output_dir = 'results/'

ref_catalogue_ch = Channel.value("${params.ref_catalogue}")

raw_image_ch = Channel.fromPath("${params.input_dir}/Epoch0?.fits").map{ it -> [it.baseName, it]}


process bane_raw {
  publishDir params.output_dir, mode:'copy', overwrite:false

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

  publishDir params.output_dir, mode:'copy', overwrite:false

  input:
  set val(basename), file(image), file(bkg), file(rms) from raw_image_with_bkg_ch

  output:
  set val(basename), file(image), file("${basename}_comp.fits") into initial_catalogue_ch

  script:
  """
  echo aegean --background=${bkg} --noise=${rms} --table=${image} ${image}
  touch ${basename}_comp.fits
  """
}

// initial_catalogue_ch.view()

process fits_warp {
  publishDir params.output_dir, mode:'copy', overwrite:false

  input:
  set val(basename), file(image), file(catalogue) from initial_catalogue_ch
  path rfile from params.region_file

  output:
  file("${basename}_warped.fits") into warped_images_ch
  set val("${basename}_warped"), file("${basename}_warped.fits") into warped_images_ch2

  script:
  """
  echo fits_warp on ${image} with ${catalogue} and ${rfile}
  touch ${basename}_warped.fits
  """
}


process make_mean_image {
  publishDir params.output_dir, mode:'copy', overwrite:false

  input:
  file(image) from warped_images_ch.collect()

  output:
  set val('mean'), file('mean.fits') into mean_image_ch

  script:
  """
  echo "do stuff with ${image}"
  touch mean.fits
  """

}

process bane_mean_image {
  publishDir params.output_dir, mode:'copy', overwrite:false

  input:
  set val(basename), file(mean) from mean_image_ch

  output:
  set val(basename), file(mean), file("${basename}_bkg.fits"), file("${basename}_rms.fits") into bane_mean_image_ch

  script:
  """
  echo bane on ${mean}
  touch ${basename}_bkg.fits
  touch ${basename}_rms.fits
  """
}

process sfind_mean_image {
  publishDir params.output_dir, mode:'copy', overwrite:false

  input:
  set val(basename), file(mean), file(bkg), file(rms) from bane_mean_image_ch

  output:
  file("${basename}_comp.fits") into mean_catalogue_sh

  script:
  """
  echo aegean --background=${bkg} --noise=${rms} --table=${mean} ${mean}
  touch ${basename}_comp.fits
  """
}

process source_monitor {
  publishDir params.output_dir, mode:'copy', overwrite:false

  input:
  file(mean_cat) from mean_catalogue_sh
  set val(basename), file(image) from warped_images_ch2

  output:
  file("${basename}_comp.fits") into priorized_catalogue_ch

  script:
  """
  echo aegean --background=${basename}_bkg.fits --noise=${basename}_rms.fits \
              --table=${image} ${image} --priorized 1 --input=${mean_cat}
  touch ${basename}_comp.fits
  """
}


process create_db {
  input:
  file(catalogue) from priorized_catalogue_ch

  output:
  val('done') into db_finished_ch

  script:
  """
  echo ingest ${catalogue} into db
  """
}


process other {
  input:
  val(whatever) from db_finished_ch.collect()

  output:
  val('done') into (stats_finished_ch, stats_finished_ch2)

  script:
  """
  echo analyse db
  """
}

process plot_lc {
  publishDir params.output_dir, mode:'copy', overwrite:false

  input:
  val(whatever) from stats_finished_ch

  output:
  path('lc_plots') into plots_ch

  script:
  """
  echo make lots of plots !
  mkdir lc_plots
  cd lc_plots
  for i in \$(seq 1 30); do touch plot\${i}.png;done
  """
}

process variable_summary_plot {
  publishDir params.output_dir, mode:'copy', overwrite:true

  input:
  val(whatever) from stats_finished_ch2

  output:
  file('variables.png') into summary_ch

  script:
  """
  echo do summary plot
  touch variables.png
  """
}