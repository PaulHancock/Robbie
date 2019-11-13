#! /usr/bin/env nextflow

params.input_dir = 'data/'
params.ref_catalogue = "master.fits"
params.region_file = "$baseDir/square.mim"
params.output_dir = 'results/'
params.image_file = 'images.txt'

params.refcat_ra = 'ra'
params.refcat_dec = 'dec'
params.warp = true

image_ch = Channel
  .fromPath(params.image_file)
  .splitText()
  .map{ it -> tuple(file(it).baseName, file(it.trim()))}
//image_ch.subscribe { print it }

process bane_raw {
  label 'bane'

// TODO: See if we can set this directive as the default behaviour
  publishDir params.output_dir, mode:'symlink', overwrite:true

  input:
  set val(basename), file(image) from image_ch
// TODO: Can we pas a dict instead of a set/list?
  output:
  set val(basename), file(image), file("*_{bkg,rms}.fits") into raw_image_with_bkg_ch

  script:
  """
  echo BANE --cores ${task.cpus} ${image}
  touch ${basename}_{bkg,rms}.fits
  """
  //echo "$image $basename"
}

// raw_image_with_bkg_ch.view()

process initial_sfind {
  label 'aegean'

  publishDir params.output_dir, mode:'symlink', overwrite:true

  input:
  set val(basename), file(image), file('*') from raw_image_with_bkg_ch

  output:
  set val(basename), file(image), file("${basename}_comp.fits") into initial_catalogue_ch

  script:
  """
  echo aegean --cores ${task.cpus} --background=*_bkg.fits --noise=*_rms.fits --table=${image} ${image}
  touch ${basename}_comp.fits
  """
}

// initial_catalogue_ch.view()


// TODO: Figure out how to skip this step
process fits_warp {
  label 'warp'

  publishDir params.output_dir, mode:'symlink', overwrite:false

  input:
  set val(basename), file(image), file(catalogue) from initial_catalogue_ch
  path rfile from params.region_file

  output:
  file("${basename}_warped.fits") into warped_images_ch
  set val("${basename}_warped"), file("${basename}_warped.fits") into (warped_images_ch2, warped_images_ch3)

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
  """
}


process make_mean_image {
  publishDir params.output_dir, mode:'symlink', overwrite:false

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

  publishDir params.output_dir, mode:'symlink', overwrite:false

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
  label 'aegean'

  publishDir params.output_dir, mode:'symlink', overwrite:false

  input:
  set val(basename), file(mean), file(bkg), file(rms) from bane_mean_image_ch

  output:
  file("${basename}_comp.fits") into (mean_catalogue_ch, mean_catalogue_ch2)

  script:
  """
  echo aegean --background=${bkg} --noise=${rms} --table=${mean} ${mean}
  touch ${basename}_comp.fits
  """
}

// TODO: Join the monitoring list to the persistent source catalogue

process source_monitor {
  label 'aegean'

  publishDir params.output_dir, mode:'symlink', overwrite:false

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
  publishDir params.output_dir, mode:'symlink', overwrite:false

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
  publishDir params.output_dir, mode:'symlink', overwrite:false

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
  publishDir params.output_dir, mode:'symlink', overwrite:false

  input:
  file(mean_cat) from mean_catalogue_ch2
  set val(basename), file(warped) from warped_images_ch3

  output:
  set val("${basename}_masked"), file("${basename}_masked.fits") into masked_images_ch

  script:
  """
  echo aeres -c ${mean_cat} -f ${warped} -r ${basename}_masked.fits --add
  touch ${basename}_masked.fits
  """
}

process sfind_masked {
  label 'aegean'

  publishDir params.output_dir, mode:'symlink', overwrite:false

  input:
  path rfile from params.region_file
  set val(basename), file(masked) from masked_images_ch

  output:
  file("${basename}_comp.fits") into masked_catalogue_ch

  script:
  // Do a hack to figure out the name of the bkg/rms files since recomputing them isn't an option
  """
  base="${basename}"
  base="\${base%%_warped_masked}"
  echo aegean --background=\${base}_bkg.fits --noise=\${base}_rms.fits --table=${masked} ${masked}
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
  publishDir params.output_dir, mode:'symlink', overwrite:false

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