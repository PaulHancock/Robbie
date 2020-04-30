#! /usr/bin/env nextflow

/* CONFIGURATION STAGE */

/* for executing without a container */
params.codeDir = "${baseDir}/"
params.stilts = ""

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
params.monitoring_src_file = ''

// Source finding region file
params.use_region_file = false
params.region_file = ""

// set NO_FILE as the filename if region/warping is turned off
if (params.use_monitoring_src_file == false) {
   params.monitoring_src_file='NO_FILE'
}
if (params.use_region_file == false) {
   params.region_file='NO_FILE'
}

log.info """\
         ROBBIE the Space Detective 
         ==========================
         minotor src  : ${params.use_monitoring_src_file} / ${params.monitoring_src_file}
         region file  : ${params.use_region_file} / ${params.region_file}
         output to    : ${params.output_dir}
         --
         run as       : ${workflow.commandLine}
         config files : ${workflow.configFiles}
         container    : ${workflow.containerEngine}:${workflow.container}
         """
         .stripIndent()

/*
mean_catalogue_ch = Channel
  .fromPath("/astro/mwasci/phancock/nxf/Robbie/results/${params.timescale}/persistent_sources.fits")

mean_catalogue_ch3 = Channel
  .fromPath("/astro/mwasci/phancock/nxf/Robbie/results/${params.timescale}/persistent_sources.fits")

warped_images_ch = Channel
  .fromPath("/astro/mwasci/phancock/nxf/Robbie/results/${params.timescale}/*warped.fits")
  .map{ it -> tuple(file(it).baseName, file(it))}
*/

all_images_ch = Channel
  .fromFilePairs("/astro/mwasci/phancock/nxf/Robbie/results/${params.timescale}/*_warped{,_comp}.fits", size:-1)


process update_catalogues {
  label 'python'

  input:
  tuple val(ignore), path('*') from all_images_ch

  output:
  path("*_comp_fix.fits") into priorized_catalogue_ch2

  script:
  """
  echo ${task.process} on \${HOSTNAME}
  # super hack to get stilts to play nice and add two columns of strings
  basename=\$(ls *_warped.fits)
  epoch=\$(get_epoch.py \${basename})
  epoch="\\\\\\\"\${epoch}\\\\\\\""
  filename="\\\\\\\"\${basename}.fits\\\\\\\""
  ${params.stilts} tpipe in=\${basename%%.fits}_comp.fits cmd="addcol image \${filename}" \
                                         cmd="addcol epoch \${epoch}" \
                                         ofmt=fits out=temp.fits
  mv temp.fits \${basename%%.fits}_comp_fix.fits
  """
}


process join_fluxes {
  label 'python'

  input:
  path('*') from priorized_catalogue_ch2.collect()
  path("reference.fits") from file("/astro/mwasci/phancock/nxf/Robbie/results/${params.timescale}/persistent_sources.fits")

  output:
  path("flux_table.vot") into flux_table_ch

  script:
  """
  echo ${task.process} on \${HOSTNAME}
  ls
  ls *_comp_fix.fits > epochs.txt
  ${params.codeDir}join_catalogues.py --refcat reference.fits --epochs epochs.txt --out flux_table.vot
  """
}

process make_GRB_lc {
  label 'python'

  input:
  path(flux_table) from flux_table_ch
 
  output:
  path('*.csv') into lc_ch

  script:
  """
  echo ${task.process} on \${HOSTNAME}
  ls
  ${params.codeDir}get_lc_from_vot.py --table ${flux_table} --name ${params.grb_name} --out ${params.grb_name}.csv
  """
}