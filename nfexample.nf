#! /usr/bin/env nextflow

params.ref_catalogue = "master.fits"
params.region_file = "$baseDir/square.mim"
params.output_dir = 'results/'
params.image_file = 'images.txt'

params.refcat_ra = 'ra'
params.refcat_dec = 'dec'
params.warp = true
params.monitor = 'monitor.fits'

// Read the image names from a text file
image_ch = Channel
  .fromPath(params.image_file)
  .splitText()
  .map{ it -> tuple(file(it).baseName, file(it.trim()))}


process bane_raw {
  label 'bane'

  input:
  set val(basename), file(image) from image_ch

  output:
  set val(basename), file(image), file("*_{bkg,rms}.fits") into raw_image_with_bkg_ch

  script:
  """
  echo BANE --cores ${task.cpus} ${image}
  touch ${basename}_{bkg,rms}.fits
  """
}


process initial_sfind {
  label 'aegean'

  input:
  set val(basename), file(image), file('*') from raw_image_with_bkg_ch

  output:
  set val(basename), file(image), file("*_{bkg,rms,comp}.fits") into initial_catalogue_ch

  script:
  """
  echo ${task.process}
  echo aegean --cores ${task.cpus} --background=*_bkg.fits --noise=*_rms.fits --table=${image} ${image}
  touch ${basename}_comp.fits
  """
}


process fits_warp {
  label 'warp'

  echo true
  input:
  set val(basename), file(image), file('*') from initial_catalogue_ch
  path rfile from params.region_file

  output:
  file("*_warped.fits") into warped_images_ch // -> to mean image
  set val("${basename}_warped"), file("*_warped.fits") into warped_images_ch2 // to monitoring
  set val("${basename}_warped"), file("*_{warped,bkg,rms}.fits") into warped_images_ch3 // to masking

  script:
  if (params.warp == true)
  """
  fits_warp.py --cores ${task.ncpus} --refcat ${params.ref_catalogue} --incat ${catalogue} \
               --ra1 ra --dec1 dec --ra2 ${params.refcat_ra} --dec2 ${params.refcat_dec} \
               --xm ${basename}_xm.fits
  fits_warp.py --infits ${image} --xm ${basename}_xm.fits --suffix warped \
               --ra2 ${params.refcat_ra} --dec2 ${params.refcat_dec} \
               --plot

  ${image} with ${catalogue} and ${rfile}
  """
  else
  """
  ln -s ${image} ${basename}_warped.fits
  echo ${task.process}
  ls *.fits
  """
}


process make_mean_image {
  input:
  file(image) from warped_images_ch.collect()

  output:
  set val('mean'), file('mean.fits') into mean_image_ch

// TODO: How do we avoid a command with an argument list of 3k files?
// TODO: Can we write the list into a text file as we did with Make?
  script:
  """
  echo "do stuff with ${image}"
  touch mean.fits
  """
}

process bane_mean_image {
  label 'bane'

  input:
  set val(basename), file(mean) from mean_image_ch

  output:
  set val(basename), file(mean), file("${basename}_{bkg,rms}.fits") into bane_mean_image_ch

  script:
  """
  echo bane on ${mean}
  touch ${basename}_bkg.fits
  touch ${basename}_rms.fits
  """
}

// Create persistent catalogue and optionally join on a monitoring list

process sfind_mean_image {
  label 'aegean'
  echo true

  input:
  set val(basename), file(mean), file('*') from bane_mean_image_ch

  output:
  file("${basename}_comp.fits") into (mean_catalogue_ch,  // to source finding 
                                      mean_catalogue_ch2) // to masking

  script:
  def mon="""
  echo ${params.stilts} tcatn nin=2 in1=${mean} in2=${params.monitor} out=persistent_sources.fits
  touch persistent_sources.fits
  """
  """
  echo aegean --background *_bkg.fits --noise *_rms.fits --table=${mean} ${mean}
  ${ (params.monitor) ? "${mon}"  : "touch persistent_sources.fits" } 
  """
}


process source_monitor {
  label 'aegean'

  input:
  file(mean_cat) from mean_catalogue_ch
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


// TODO: Future problem is that sqlite db is not good for 3M rows
// TODO: What other options are there.

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


process compute_stats {
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

process mask_images {

  echo true
  
  input:
  file(mean_cat) from mean_catalogue_ch2
  set val(basename), file('*') from warped_images_ch3

  output:
  set val("${basename}_masked"), file("*_{masked,bkg,rms}.fits") into masked_images_ch

  script:
  """
  echo ${task.process}
  ls *.fits
  echo aeres -c ${mean_cat} -f *_warped.fits -r ${basename}_masked.fits --add
  touch ${basename}_masked.fits
  """
}

process sfind_masked {
  label 'aegean'

  echo true
  
  input:
  path file from params.region_file
  set val(basename), file('*') from masked_images_ch

  output:
  file("${basename}_comp.fits") into masked_catalogue_ch
  
  script:
  """
  echo \$(ls *.fits)
  echo ${task.process}: aegean --background *_bkg.fits --noise *_rms.fits --table *_masked.fits *_masked.fits
  touch ${basename}_comp.fits
  """
}

process compile_transients_candidates {

  input:
  file(catalogue) from masked_catalogue_ch.collect()

  output:
  val('done') into transients_imported_ch

  script:
  """
  for f in ${catalogue}
  do
    echo filter on \${f}
    echo import \${f} into db
  done
  """
}

process transients_plot {

  input:
  val(whatever) from transients_imported_ch

  output:
  file('transients.png') into transients_plot_ch

  script:
  """
  echo plot transients.png
  touch transients.png
  """
}