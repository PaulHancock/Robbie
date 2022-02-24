#! /usr/bin/env nextflow
author = "Paul Hancock"
version = "2.2.0"
prev_git_hash = "7d6c6"
date = "2022-02-24"
/* CONFIGURATION STAGE */

/* for executing without a container */
params.codeDir = "${baseDir}/"
params.stilts = "stilts"

// output directory
params.output_dir = "${baseDir}/results/"

// input images are listed in this file, one image per line
params.image_file = "$baseDir/images.txt"

// Warping stage
params.warp = true
params.ref_catalogue = "$baseDir/master.fits"
params.refcat_ra = 'ra'
params.refcat_dec = 'dec'

// Plotting params
params.by_epoch = true

// monitoring of a pre-determined source
params.use_monitoring_src_file = false
params.monitoring_src_file = ""

// Source finding region file
params.use_region_file = false
params.region_file = ""


log.info """\
         ROBBIE the Space Detective 
         ==========================
         version      : ${version} (${prev_git_hash}) - ${date}
         images from  : ${params.image_file}
         do warping   : ${params.warp}
         warp ref cat : ${params.ref_catalogue}
         minotor src  : ${params.use_monitoring_src_file} / ${params.monitoring_src_file}
         region file  : ${params.use_region_file} / ${params.region_file}
         output to    : ${params.output_dir}
         --
         run as       : ${workflow.commandLine}
         config files : ${workflow.configFiles}
         container    : ${workflow.containerEngine}:${workflow.container}
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
  echo ${task.process} on \${HOSTNAME}
  BANE --cores ${task.cpus} ${image}
  """
}


process initial_sfind {
  label 'aegean'
  // echo true
  input:
  tuple val(basename), path(image), path('*') from raw_image_with_bkg_ch
  file('region.mim') from file(params.region_file)

  output:
  //tuple val(basename), path(image), path("*_{bkg,rms,comp}.fits") into initial_catalogue_ch
  tuple val(basename), path('*.fits', includeInputs:true) into initial_catalogue_ch

  script:
  def region = params.region_file != 'NO_FILE' ? "--region region.mim" : ''
  """
  echo ${task.process} on \${HOSTNAME}
  aegean --cores ${task.cpus} --background *_bkg.fits --noise *_rms.fits --table ${image} ${region} ${image}
  ls *.fits
  """
}


process fits_warp {
  label 'warp'

  //echo true
  input:
  tuple val(basename), path('*') from initial_catalogue_ch

  output:
  path("*_warped.fits") into warped_images_ch // -> to mean image
  tuple val("${basename}_warped"), path("*.fits", includeInputs:true) into warped_images_ch2 // to monitoring
  tuple val("${basename}_warped"), path('*.fits', includeInputs:true) into warped_images_ch3 // to mask_images

  script:
  suff1=(params.refcat_ra=='ra' ? '_1':'')
  suff2=(params.refcat_ra=='ra' ? '_2':'')
  
  if (params.warp == true)
  """
  echo ${task.process} on \${HOSTNAME}
  fits_warp.py --cores ${task.cpus} --refcat ${params.ref_catalogue} --incat ${basename}_comp.fits \
               --ra1 ra --dec1 dec --ra2 ${params.refcat_ra} --dec2 ${params.refcat_dec} \
               --xm ${basename}_xm.fits
  fits_warp.py --infits ${basename}.fits --xm ${basename}_xm.fits --suffix warped \
               --ra1 ra${suff1} --dec1 dec${suff1} --ra2 ${params.refcat_ra}${suff2} --dec2 ${params.refcat_dec}${suff2} \
               --plot
  ls *.fits
  """
  else
  """
  echo ${task.process} on \${HOSTNAME}
  ln -s ${basename}.fits ${basename}_warped.fits
  ls *.fits
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
  echo ${task.process} on \${HOSTNAME}
  ls *_warped.fits > images.txt
  ${params.codeDir}make_mean.py --out mean.fits --infile images.txt
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
  echo ${task.process} on \${HOSTNAME}
  BANE --cores ${task.cpus} ${mean}
  """
}

process sfind_mean_image {
  label 'aegean'

  input:
  tuple val(basename), path(mean), path('*') from bane_mean_image_ch
  path('region.mim') from params.region_file
  path('monitor.fits') from params.monitoring_src_file
  
  output:
  path("persistent_sources.fits") into (mean_catalogue_ch,  // to source_monitor
                                       mean_catalogue_ch2,  // to mask_images
                                       mean_catalogue_ch3)  // to join_fluxes

  script:
  def region = params.use_region_file ? "--region region.mim":''
  def mon = params.use_monitoring_src_file ? """
  ${params.stilts} tcatn nin=2 in1=mean_comp.fits in2=monitor.fits out=persistent_sources.fits ofmt=fits
  """ : "mv *_comp.fits persistent_sources.fits" 

  """
  echo ${task.process} on \${HOSTNAME}
  aegean --cores ${task.cpus} --background *_bkg.fits --noise *_rms.fits --table ${mean} ${region} ${mean}
  ${mon} 
  """
}


process source_monitor {
  label 'aegean'

  input:
  path(mean_cat) from mean_catalogue_ch.collect()
  tuple val(basename), path(image) from warped_images_ch2

  output:
  tuple path("${basename}_comp.fits"), path(image, includeInputs:true) into priorized_catalogue_ch
  path("${basename}_comp.fits") into priorized_catalogue_ch2

  script:
  """
  echo ${task.process} on \${HOSTNAME}
  aegean --cores ${task.cpus} --background *_bkg.fits --noise *_rms.fits --noregroup\
         --table ${basename}.fits --priorized 1 --input ${mean_cat} ${basename}.fits

  # super hack to get stilts to play nice and add two columns of strings
  epoch=\$(${params.codeDir}get_epoch.py ${basename}.fits)
  epoch="\\\\\\\"\${epoch}\\\\\\\""
  filename="\\\\\\\"${basename}.fits\\\\\\\""
  ${params.stilts} tpipe in=${basename}_comp.fits cmd="addcol image \${filename}" \
                                                  cmd="addcol epoch \${epoch}" \
                                                  ofmt=fits out=temp.fits
  mv temp.fits ${basename}_comp.fits
  """
}


process join_fluxes {
  label 'python'

  input:
  path('*') from priorized_catalogue_ch2.collect()
  path("reference.fits") from mean_catalogue_ch3

  output:
  path("flux_table.vot") into flux_table_ch

  script:
  """
  echo ${task.process} on \${HOSTNAME}
  ls
  ls *_comp.fits > epochs.txt
  ${params.codeDir}join_catalogues.py --refcat reference.fits --epochs epochs.txt --out flux_table.vot
  """
}


process compute_stats {
  label 'python'

  input:
  path("flux_table.vot") from flux_table_ch

  output:
  tuple path("flux_table.vot", includeInputs:true), path("stats_table.vot") into flux_stats_ch

  script:
  """
  echo ${task.process} on \${HOSTNAME}
  ls
  NDOF=(\$(${params.codeDir}auto_corr.py --table flux_table.vot))
  echo \${NDOF[@]}
  echo \${NDOF[@]} > NDOF.txt
  ${params.codeDir}calc_var.py --table flux_table.vot --ndof \${NDOF[-1]} --out stats_table.vot --cores ${task.cpus}
  """
}

process plot_lc {
  label 'python'

  input:
  tuple path("flux_table.vot"), path("stats_table.vot") from flux_stats_ch

  output:
  path('plots') into plots_ch

  script:
  """
  echo ${task.process} on \${HOSTNAME}
  mkdir plots
  ${params.codeDir}plot_variables.py --ftable flux_table.vot --stable stats_table.vot --plot plots/variables.png --all --cores ${task.cpus}
  """
}


process mask_images {
  label 'python'
  // echo true
  
  input:
  path(mean_cat) from mean_catalogue_ch2.collect()
  tuple val(basename), path('*') from warped_images_ch3

  output:
  tuple val("${basename}_masked"), path("*.fits", includeInputs:true) into masked_images_ch

  script:
  """
  echo ${task.process} on \${HOSTNAME}
  ls *.fits
  AeRes -c ${mean_cat} -f *_warped.fits -r ${basename}_masked.fits --mask --sigma 0.1
  """
}

process sfind_masked {
  label 'aegean'

  // echo true
  
  input:
  tuple val(basename), path('*') from masked_images_ch
  path('region.mim') from file(params.region_file)

  output:
  path("${basename}_comp.fits") optional true into masked_catalogue_ch
  
  script:
  def region = (params.region_file != 'NO_FILE' ? "--region region.mim":'')
  """
  echo ${task.process} on  \${HOSTNAME}
  ls *
  aegean --cores ${task.cpus} --background *_bkg.fits --noise *_rms.fits --table ${basename}.fits ${basename}.fits ${region}
  # Don't filter if there is no output
  if [[ -e ${basename}_comp.fits ]]
  then
    ${params.codeDir}filter_transients.py --incat ${basename}_comp.fits --image ${basename}.fits --outcat out.fits
    if [[ -e out.fits ]]
    then
      mv out.fits ${basename}_comp.fits
    fi
  fi
  ls *.fits
  """
}

process compile_transients_candidates {
  label 'python'

  input:
  path(catalogue) from masked_catalogue_ch.collect()

  output:
  path('transients.fits') into transients_imported_ch

  script:
  """
  echo ${task.process} on \${HOSTNAME}
  ls *_comp.fits > temp.dat
  ${params.codeDir}collect_transients.py --infile temp.dat --out transients.fits --ignoremissing
  """
}

process transients_plot {
  label 'python'
  
  input:
  path(transients) from transients_imported_ch

  output:
  path('transients.png') into transients_plot_ch

  script:
  """
  echo ${task.process} on \${HOSTNAME}
  ${params.codeDir}plot_transients.py --in ${transients} --plot transients.png
  """
}
