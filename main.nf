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

// Plotting params
params.by_epoch = true

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
  echo ${task.process}
  echo \${HOSTNAME}
  ls -lh *
  pwd
  BANE --cores ${task.cpus} ${image}
  # touch ${basename}_{bkg,rms}.fits
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
  aegean --cores ${task.cpus} --background *_bkg.fits --noise *_rms.fits --table ${image} ${image}
  # touch ${basename}_comp.fits
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
  echo ${task.process}
  fits_warp.py --cores ${task.ncpus} --refcat ${params.ref_catalogue} --incat ${catalogue} \
               --ra1 ra --dec1 dec --ra2 ${params.refcat_ra} --dec2 ${params.refcat_dec} \
               --xm ${basename}_xm.fits
  fits_warp.py --infits ${basename}.fits --xm ${basename}_xm.fits --suffix warped \
               --ra2 ${params.refcat_ra} --dec2 ${params.refcat_dec} \
               --plot
  ls *.fits
  echo \${HOSTNAME}
  """
  else
  """
  echo ${task.process}
  ln -s ${basename}.fits ${basename}_warped.fits
  ls *.fits
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
  echo ${task.process}
  ls *_warped.fits > images.txt
  make_mean.py --out mean.fits --infile images.txt
  # touch mean.fits
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
  echo ${task.process}
  BANE --cores ${task.cpus} ${mean}
  # touch ${basename}_{bkg,rms}.fits
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
  ${params.stilts} tcatn nin=2 in1=mean_comp.fits in2=${params.monitor} out=persistent_sources.fits ofmt
  # touch persistent_sources.fits
  """

  """
  echo ${task.process}
  aegean --background *_bkg.fits --noise *_rms.fits --table ${mean} ${mean}
  ${ (params.monitor) ? "${mon}"  : "mv *_comp.fits persistent_sources.fits" } 
  echo \${HOSTNAME}
  """
}


process source_monitor {
  label 'aegean'

  input:
  path(mean_cat) from mean_catalogue_ch
  tuple val(basename), path(image) from warped_images_ch2

  output:
  tuple path("${basename}_comp.fits"), path(image, includeInputs:true) into priorized_catalogue_ch

  script:
  """
  echo ${task.process}
  aegean --background *_bkg.fits --noise *_rms.fits \
         --table ${basename}.fits --priorized 1 --input ${mean_cat} ${basename}.fits
  # touch ${basename}_comp.fits
  echo \${HOSTNAME}
  """
}


// TODO: Future problem is that sqlite db is not good for 3M rows
// TODO: What other options are there.

process make_db {
  label 'python'

  output:
  path '*.db' into (db_file_ch1, db_file_ch2, db_file_ch3)

  script:
  """
  remake_db.py --name ${params.db_file}
  """
}

process populate_db {
  label 'python'

  input:
  tuple path(cat), path(image) from priorized_catalogue_ch
  path 'database.db' from db_file_ch1

  output:
  val('done') into db_finished_ch

  script:
  """
  echo ${task.process}
  add_cat_to_db.py --name database.db --cat ${cat} --image *_warped.fits
  echo \${HOSTNAME}
  """
}

// TODO: Figure out how to make this not break when we do -resume

process compute_stats {
  label 'python'

  input:
  val(whatever) from db_finished_ch.collect()
  path 'database.db' from db_file_ch2

  output:
  val('done') into stats_finished_ch
  path('NDOF.txt') into ndof_publish_ch

  script:
  """
  echo ${task.process}
  NDOF=(\$(auto_corr.py --dbname database.db))
  echo \${NDOF[@]} > NDOF.txt
  calc_var.py --name database.db --ndof \${NDOF[-1]} 
  echo \${HOSTNAME}
  """
}

process plot_lc {
  label 'python'

  input:
  val(whatever) from stats_finished_ch
  path 'database.db' from db_file_ch3

  output:
  path('plots') into plots_ch

  script:
  """
  echo ${task.process}
  mkdir plots
  plot_variables.py --name database.db --plot plots/variables.png --all ${ (params.by_epoch)? '': '--dates'}
  echo \${HOSTNAME}
  """
}


process mask_images {
  label 'python'
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
  AeRes -c ${mean_cat} -f *_warped.fits -r ${basename}_masked.fits --add
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
  path("${basename}_comp.fits") optional true into masked_catalogue_ch
  
  script:
  """
  echo ${task.process}
  echo \${HOSTNAME}
  aegean --background *_bkg.fits --noise *_rms.fits --table ${basename}.fits ${basename}.fits
  filter_transients.py --incat ${basename}_comp.fits --image ${basename}.fits --outcat ${basename}_comp.fits
  # Remove  the output file if it's empty
  nsrc=(\$( stilts tcat omode=count in=${basename}_comp.fits ))
  if [[ \${nsrc[-1]} -lt 1 ]]
  then
    rm *_comp.fits
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
  echo ${task.process}
  echo \${HOSTNAME}
  ls *_comp.fits > temp.dat
  collect_transients.py --infile temp.dat --out transients.fits --ignoremissing
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
  echo ${task.process}
  plot_transients.py --in ${transients} --plot transients.png
  echo \${HOSTNAME}
  """
}