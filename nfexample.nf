#! /usr/bin/env nextflow

/* CONFIGURATION STAGE */

// input images are listed in this file, one image per line
params.image_file = "$baseDir/images.txt"

// database to use
params.db = 'sqlite3'
params.db_file = "$baseDir/flux_table.db"

// Warping stage
params.warp = true
params.ref_catalogue = "$baseDir/master.fits"
params.refcat_ra = 'ra'
params.refcat_dec = 'dec'

// name of monitoring file - set to null if not required
params.monitor='monitor.fits'

// calling stilts
params.stilts = "java -jar /data/nextflow/stilts/stilts.jar"

// Source finding params
params.region_file = "$baseDir/square.mim"

// output directory
params.output_dir = 'results/'

log.info """\
         ROBBIE the Space Detective 
         ==========================
         images from  : ${params.image_file}
         using db     : ${params.db_file} (type=${params.db})
         do warping?  : ${params.warp}
         ref cat      : ${params.ref_catalogue}
         minotor cat  : ${params.monitor}
         region file  : ${params.region_file}
         output to    : ${params.output_dir}
         """
         .stripIndent()

// Read the image names from a text file
image_ch = Channel
  .fromPath(params.image_file)
  .splitText()
  .map{ it -> tuple(file(it).baseName, file(it.trim()))}


process bane_raw {
  label 'bane'

  // echo true
  input:
  tuple val(basename), path(image) from image_ch

  output:
  tuple val(basename), path(image), path("*_{bkg,rms}.fits") into raw_image_with_bkg_ch

  script:
  """
  echo BANE --cores ${task.cpus} ${image}
  touch ${basename}_{bkg,rms}.fits
  ls *.fits
  echo \${HOSTNAME}
  """
}


process initial_sfind {
  label 'aegean'
  // echo true
  input:
  tuple val(basename), path(image), path('*') from raw_image_with_bkg_ch

  output:
  //tuple val(basename), path(image), path("*_{bkg,rms,comp}.fits") into initial_catalogue_ch
  tuple val(basename), path('*.fits', includeInputs:true) into initial_catalogue_ch

  script:
  """
  echo ${task.process}
  echo aegean --cores ${task.cpus} --background=*_bkg.fits --noise=*_rms.fits --table=${image} ${image}
  touch ${basename}_comp.fits
  ls *.fits
  echo \${HOSTNAME}
  """
}


process fits_warp {
  label 'warp'

  //echo true
  input:
  tuple val(basename), path('*') from initial_catalogue_ch
  path rfile from params.region_file

  output:
  path("*_warped.fits") into warped_images_ch // -> to mean image
  tuple val("${basename}_warped"), path("*.fits", includeInputs:true) into warped_images_ch2 // to monitoring
  tuple val("${basename}_warped"), path('*.fits', includeInputs:true) into warped_images_ch3 // to mask_images

  script:
  if (params.warp == true)
  """
  fits_warp.py --cores ${task.ncpus} --refcat ${params.ref_catalogue} --incat ${catalogue} \
               --ra1 ra --dec1 dec --ra2 ${params.refcat_ra} --dec2 ${params.refcat_dec} \
               --xm ${basename}_xm.fits
  fits_warp.py --infits ${basename}.fits --xm ${basename}_xm.fits --suffix warped \
               --ra2 ${params.refcat_ra} --dec2 ${params.refcat_dec} \
               --plot

  echo ${basename}.fits with ${catalogue} and ${rfile}
  """
  else
  """
  ln -s ${basename}.fits ${basename}_warped.fits
  echo ${task.process}
  ls *.fits
  echo \${HOSTNAME}
  echo \${HOSTNAME}
  """
}


process make_mean_image {

  // echo true
  input:
  path(image) from warped_images_ch.collect()

  output:
  tuple val('mean'), path('mean.fits') into mean_image_ch

  script:
  """
  ls *_warped.fits > images.txt
  echo make_mean.py --out mean.fits --infile images.txt
  touch mean.fits
  echo \${HOSTNAME}
  """
}

process bane_mean_image {
  label 'bane'

  input:
  tuple val(basename), path(mean) from mean_image_ch

  output:
  tuple val(basename), path(mean), path("${basename}_{bkg,rms}.fits") into bane_mean_image_ch

  script:
  """
  echo bane on ${mean}
  touch ${basename}_bkg.fits
  touch ${basename}_rms.fits
  echo \${HOSTNAME}
  """
}

process sfind_mean_image {
  label 'aegean'

  input:
  tuple val(basename), path(mean), path('*') from bane_mean_image_ch

  output:
  path("persistent_sources.fits") into (mean_catalogue_ch,  // to source finding 
                                       mean_catalogue_ch2)  // to masking

  script:
  def mon="""
  echo ${params.stilts} tcatn nin=2 in1=${mean} in2=${params.monitor} out=persistent_sources.fits
  touch persistent_sources.fits
  """
  """
  echo aegean --background *_bkg.fits --noise *_rms.fits --table=${mean} ${mean}
  ${ (params.monitor) ? "${mon}"  : "touch persistent_sources.fits" } 
  echo \${HOSTNAME}
  """
}


process source_monitor {
  label 'aegean'

  input:
  path(mean_cat) from mean_catalogue_ch
  tuple val(basename), path(image) from warped_images_ch2

  output:
  path("${basename}_comp.fits") into priorized_catalogue_ch

  script:
  """
  echo aegean --background=${basename}_bkg.fits --noise=${basename}_rms.fits \
              --table=${image} ${image} --priorized 1 --input=${mean_cat}
  touch ${basename}_comp.fits
  echo \${HOSTNAME}
  """
}


// TODO: Future problem is that sqlite db is not good for 3M rows
// TODO: What other options are there.

process create_db {
  input:
  path(catalogue) from priorized_catalogue_ch

  output:
  val('done') into db_finished_ch

  script:
  """
  echo ingest ${catalogue} into db
  echo \${HOSTNAME}
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
  echo \${HOSTNAME}
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
  echo \${HOSTNAME}
  """
}

process variable_summary_plot {

  input:
  val(whatever) from stats_finished_ch2

  output:
  path('variables.png') into summary_ch

  script:
  """
  echo do summary plot
  touch variables.png
  echo \${HOSTNAME}
  """
}

process mask_images {

  // echo true
  
  input:
  path(mean_cat) from mean_catalogue_ch2
  tuple val(basename), path('*') from warped_images_ch3

  output:
  tuple val("${basename}_masked"), path("*.fits", includeInputs:true) into masked_images_ch

  script:
  """
  echo ${task.process}
  ls *.fits
  echo aeres -c ${mean_cat} -f *_warped.fits -r ${basename}_masked.fits --add
  touch ${basename}_masked.fits
  echo \${HOSTNAME}
  """
}

process sfind_masked {
  label 'aegean'

  // echo true
  
  input:
  path file from params.region_file
  tuple val(basename), path('*') from masked_images_ch

  output:
  path("${basename}_comp.fits") into masked_catalogue_ch
  
  script:
  """
  echo ${task.process}
  echo  aegean --background *_bkg.fits --noise *_rms.fits --table *_masked.fits *_masked.fits
  touch ${basename}_comp.fits
  ls *.fits
  echo \${HOSTNAME}
  """
}

process compile_transients_candidates {

  input:
  path(catalogue) from masked_catalogue_ch.collect()

  output:
  val('done') into transients_imported_ch

  script:
  """
  for f in ${catalogue}
  do
    echo filter on \${f}
    echo import \${f} into db
  done
  echo \${HOSTNAME}
  """
}

process transients_plot {

  input:
  val(whatever) from transients_imported_ch

  output:
  path('transients.png') into transients_plot_ch

  script:
  """
  echo plot transients.png
  touch transients.png
  echo \${HOSTNAME}
  """
}