#! /usr/bin/env nextflow

nextflow.enable.dsl=2

author = "Paul Hancock"

if ( params.help ) {
    help = """Robbie.nf: A batch processing work-flow for the detection of radio transients and
             |        variables
             |Standard argurments:
             |  --image_file  A text file where each line is the location of an image fits file.
             |                [default: ${params.image_file}]
             |  --stilts      The command required to run stilts.
             |                Eg. "java -jar ~/Downloads/stilts.jar"
             |                [default: ${params.stilts}]
             |  --convolve
             |                Determine the smallest psf common to all input images and then convolve all images
             |                to this psf prior to any other processing [default: ${params.convolve}]
             |
             |Warping arguments:
             |  --fits_warp   Use astrometric correction via fits_warp.py.
             |                [default: ${params.fits_warp}]
             |  --flux_warp   Use flux density correction via flux_warp.
             |                [default: ${params.flux_warp}]
             |  --ref_catalogue
             |                The reference catalogue to warp your images to match.
             |                [default: will download and use GLEAM catalogue]
             |  --refcat_ra   The label the reference catalogue uses for Right Acension.
             |                [default: ${params.refcat_ra}]
             |  --refcat_dec  The label the reference catalogue uses for Declination.
             |                [default: ${params.refcat_dec}]
             |
             |Monitoring arguments:
             |  --use_monitoring_src_file
             |                Use the monitoring source file. [default: ${params.use_monitoring_src_file}]
             |  --monitoring_src_file
             |                The location of the monitoring source file. [default: ${params.monitoring_src_file}]
             |
             |Source finding arguments:
             | --use_region_file
             |                Use a source finding file. [default: ${params.use_region_file}]
             | --region_file  The location of the source finding file. [default: ${params.region_file}]
             |
             |Directory arguments:
             |  --output_dir  The directory to output the results to.
             |                [default: ${params.output_dir}]
             |  --keep_epoch_images
             |                Keep the epoch images after warping.
             |                [default: ${params.keep_epoch_images}]
             |  -w            The Nextflow work directory. Delete the directory once the processs
             |                is finished [default: ${workDir}]""".stripMargin()
    println(help)
    exit(0)
}


log.info """\
         ROBBIE the Space Detective
         ==========================
         images from  : ${params.image_file}
         convolve img : ${params.convolve}
         fits warp    : ${params.fits_warp}
         flux warp    : ${params.flux_warp}
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
  .fromPath( params.image_file )
  .splitCsv()
  .map{ it -> tuple(file(it[0]).baseName, file(it[0])) }

// Set up optional commands
if ( params.use_monitoring_src_file ) {
  monitoring_src_file = Channel.fromPath( params.monitoring_src_file )
  monitoring_command = "${params.stilts} tcatn nin=2 in1=mean_comp.fits in2=${monitoring_src_file} out=persistent_sources.fits ofmt=fits"
}
else {
  monitoring_command =  "mv *_comp.fits persistent_sources.fits"
}
if ( params.use_region_file ) {
  region_file = Channel.fromPath( params.region_file )
  region_command = "--region ${region_file}"
}
else {
  region_command = ""
}


process get_version {
  publishDir params.output_dir, mode: 'copy'

  output:
  path "version.txt"

  """
  robbie_version.sh > version.txt
  """
}


process convolve_beams {
  input:
  path(image)

  output:
  path("*_convolved.fits")

  script:
  """
  echo ${task.process} on \${HOSTNAME}

  convol_common_resolution.py --in ${image}
  """
}


process bane_raw {
  label 'bane'

  input:
  tuple val(basename), path(image)

  output:
  tuple val(basename), path(image), path("*_{bkg,rms}.fits")

  script:
  """
  echo ${task.process} on \${HOSTNAME}
  BANE --cores ${task.cpus} ${image}
  """
}


process initial_sfind {
  label 'aegean'

  input:
  tuple val(basename), path(image), path(bkg_rms_fits)

  output:
  tuple val(basename), path('*.fits', includeInputs:true)

  script:
  """
  echo ${task.process} on \${HOSTNAME}
  aegean --background *_bkg.fits --noise *_rms.fits --table ${image} ${region_command} ${image}
  ls *.fits
  """
}

process download_gleam_catalogue {
  output:
  path("*fits")

  script:
  """
  #! /usr/bin/env python

  import os
  from robbie import data_load
  from astroquery.vizier import Vizier

  if os.path.exists(data_load.REF_CAT):
    # Already exists so just sym link
    os.symlink(data_load.REF_CAT, "GLEAM_ref_cat.fits")
  else:
    # Download it from Vizier
    cat = Vizier(catalog="VIII/100/gleamegc", columns=['GLEAM', 'RAJ2000', 'DEJ2000', 'Fpwide', 'Fintwide'], row_limit=-1).query_constraints()[0]
    try:
      cat.write(data_load.REF_CAT, format='fits')
      os.symlink(data_load.REF_CAT, "GLEAM_ref_cat.fits")
    except OSError:
      # No permission so dump it here
      cat.write("GLEAM_ref_cat.fits", format='fits')
  """
}


process fits_warp {
  label 'warp'
  publishDir params.output_dir, mode: 'copy', pattern: "*_warped.fits", enabled: params.keep_epoch_images

  input:
  tuple val(basename), path(initial_catalogue)
  each ref_catalogue

  output:
  tuple val(basename), path("*_warped.fits")

  script:
  suff1=(params.refcat_ra=='ra' ? '_1':'')
  suff2=(params.refcat_ra=='ra' ? '_2':'')

  """
  echo ${task.process} on \${HOSTNAME}
  fits_warp.py --cores ${task.cpus} --refcat ${ref_catalogue} --incat ${basename}_comp.fits \
               --ra1 ra --dec1 dec --ra2 ${params.refcat_ra} --dec2 ${params.refcat_dec} \
               --xm ${basename}_xm.fits
  fits_warp.py --infits ${basename}.fits --xm ${basename}_xm.fits --suffix warped \
               --ra1 ra${suff1} --dec1 dec${suff1} --ra2 ${params.refcat_ra}${suff2} --dec2 ${params.refcat_dec}${suff2} \
               --plot
  ls *.fits
  """
}


process flux_warp {
  label 'warp'
  publishDir params.output_dir, mode: 'copy', pattern: "*_warped.fits", enabled: params.keep_epoch_images

  input:
  tuple val(basename), path(initial_catalogue)
  each ref_catalogue

  output:
  tuple val(basename), path("*_warped.fits")

  script:
  suff1=(params.refcat_ra=='ra' ? '_1':'')
  suff2=(params.refcat_ra=='ra' ? '_2':'')

  """
  echo ${task.process} on \${HOSTNAME}
  match_catalogues ${initial_catalogue} ${ref_catalogue} -o matched.fits --ra2 RAJ2000 --dec2 DEJ2000
  flux_warp matched.fits ${initial_catalogue} -o ${basename}_flux.fits
  """
}

process make_mean_image {
  publishDir params.output_dir, mode: 'copy'

  input:
  path image

  output:
  tuple val('mean_image'), path('mean_image.fits')
  // This is just to publish the (warped) image
  path(image)

  script:
  """
  echo ${task.process} on \${HOSTNAME}

  ls *.fits > images.txt
  ${params.swarp} -d > swarp.config
  ${params.swarp} @images.txt -c swarp.config \
                  -SUBTRACT_BACK N \
                  -PROJECTION_TYPE SIN \
                  -COMBINE_TYPE MEDIAN \
                  -IMAGEOUT_NAME mean_image.fits \
                  -COPY_KEYWORDS BPA,BMAJ,BMIN,FREQ
  """
}

process make_sky_coverage {
  publishDir params.output_dir, mode: 'copy'

  input:
  path image

  output:
  path('sky_coverage.fits')

  script:
  """
  echo ${task.process} on \${HOSTNAME}

  ls *.fits > images.txt

  for f in \$(ls *.fits); do make_weights.py \${f}; done

  ${params.swarp} -d > swarp.config
  ${params.swarp} @images.txt -c swarp.config \
                  -SUBTRACT_BACK N \
                  -PROJECTION_TYPE SIN \
                  -COMBINE_TYPE WEIGHTED \
                  -WEIGHTOUT_NAME sky_coverage.fits \
                  -WEIGHT_TYPE MAP_WEIGHT \
                  -RESCALE_WEIGHTS N 
  rm *.weight.fits

  python - <<EOF
  from astropy.io import fits 
  import numpy as np 
  hdu = fits.open('sky_coverage.fits') 
  hdu[0].data = np.array(np.round(hdu[0].data), dtype=np.int8) 
  hdu.writeto('sky_coverage.fits', overwrite=True) 
  EOF
  """
}

process bane_mean_image {
  label 'bane'
  publishDir params.output_dir, mode: 'copy'

  input:
  tuple val(basename), path(mean)

  output:
  tuple val(basename), path(mean), path("${basename}_{bkg,rms}.fits")

  script:
  """
  echo ${task.process} on \${HOSTNAME}
  BANE --cores ${task.cpus} ${mean}
  """
}


process sfind_mean_image {
  label 'aegean'

  input:
  tuple val(basename), path(mean), path(bkg_rms_fits)

  output:
  path "persistent_sources.fits"

  script:
  """
  echo ${task.process} on \${HOSTNAME}
  aegean --cores ${task.cpus} --background *_bkg.fits --noise *_rms.fits --table ${mean} ${region_command} ${mean}
  ${monitoring_command}
  """
}


process source_monitor {
  label 'aegean'
  maxForks 1

  input:
  path mean_cat
  tuple val(basename), path(image), path(bkg_rms)

  output:
  tuple path("${image.baseName}_comp.fits"), path(image, includeInputs:true)
  path "${image.baseName}_comp.fits"

  script:
  """
  echo ${task.process} on \${HOSTNAME}
  aegean --cores ${task.cpus} --background *_bkg.fits --noise *_rms.fits --noregroup\
         --table ${image} --priorized 1 --input ${mean_cat} ${image}

  # super hack to get stilts to play nice and add two columns of strings
  epoch=\$(get_epoch.py ${image})
  epoch="\\\\\\\"\${epoch}\\\\\\\""
  filename="\\\\\\\"${image}\\\\\\\""
  ${params.stilts} tpipe in=${image.baseName}_comp.fits cmd="addcol image \${filename}" \
                                                  cmd="addcol epoch \${epoch}" \
                                                  ofmt=fits out=temp.fits
  mv temp.fits ${image.baseName}_comp.fits
  """
}


process join_fluxes {
  label 'python'

  input:
  path source_monitor_cat
  path reference_fits

  output:
  path "flux_table.vot"

  script:
  """
  echo ${task.process} on \${HOSTNAME}
  ls
  ls *_comp.fits > epochs.txt
  join_catalogues.py --refcat ${reference_fits} --epochs epochs.txt --out flux_table.vot
  """
}


process compute_stats {
  label 'python'
  publishDir params.output_dir, mode: 'copy'

  input:
  path flux_table

  output:
  tuple path("flux_table.vot", includeInputs:true), path("stats_table.vot")

  script:
  """
  echo ${task.process} on \${HOSTNAME}
  ls
  NDOF=(\$(auto_corr.py --table ${flux_table}))
  echo \${NDOF[@]}
  echo \${NDOF[@]} > NDOF.txt
  calc_var.py --table flux_table.vot --ndof \${NDOF[-1]} --out stats_table.vot --cores ${task.cpus}
  """
}


process plot_lc {
  label 'python'
  publishDir params.output_dir, mode: 'copy'

  input:
  tuple path(flux_table), path(stats_table)

  output:
  path 'variables.png'
  path 'light_curve_plots'

  script:
  dates=(params.plotdates?"--dates":"")
  """
  echo ${task.process} on \${HOSTNAME}
  mkdir light_curve_plots
  plot_variables.py --ftable ${flux_table} \
                    --stable ${stats_table} \
                    --plot variables.png \
                    --lc_dir light_curve_plots \
                    --all \
                    --cores ${task.cpus} ${dates}
  """
}


process mask_images {
  label 'python'

  input:
  path mean_cat
  tuple val(basename), path(image), path(bkg_rms)

  output:
  tuple val(basename), path("${basename}_masked.fits"), path(bkg_rms)

  script:
  """
  echo ${task.process} on \${HOSTNAME}
  ls *.fits
  AeRes -c ${mean_cat} -f ${image} -r ${basename}_masked.fits --mask --sigma 0.1
  """
}


process sfind_masked {
  label 'aegean'

  input:
  tuple val(basename), path(masked_images), path(bkg_rms)

  output:
  path "${masked_images.baseName}_comp.fits" optional true

  script:
  """
  echo ${task.process} on  \${HOSTNAME}
  ls *
  aegean --cores ${task.cpus} --background *_bkg.fits --noise *_rms.fits --table ${masked_images} ${masked_images} ${region_command}
  # Don't filter if there is no output
  if [[ -e ${masked_images.baseName}_comp.fits ]]
  then
    filter_transients.py --incat ${masked_images.baseName}_comp.fits --image ${masked_images} --outcat out.fits
    if [[ -e out.fits ]]
    then
      mv out.fits ${masked_images.baseName}_comp.fits
    fi
  fi
  ls *.fits
  """
}


process compile_transients_candidates {
  label 'python'
  publishDir params.output_dir, mode: 'copy'

  input:
  path catalogue

  output:
  path 'transients.fits'

  script:
  """
  echo ${task.process} on \${HOSTNAME}
  ls *_comp.fits > temp.dat
  collect_transients.py --infile temp.dat --out transients.fits --ignoremissing
  """
}


process transients_plot {
  label 'python'
  publishDir params.output_dir, mode: 'copy'

  input:
  path transients

  output:
  path 'transients.png'

  script:
  """
  #! /usr/bin/env python

  from astropy.table import Table
  import numpy as np
  import matplotlib
  from matplotlib import pyplot
  from matplotlib.patches import Ellipse
  import socket

  print(f"${task.process} on {socket.gethostname()}")

  tab = Table.read("${transients}")
  nepochs = np.max(tab['epoch'])*1.
  kwargs={'fontsize':14}

  cmap = pyplot.cm.plasma_r
  # define the bins and normalize
  bounds = np.linspace(3,15,13)
  norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

  fig = pyplot.figure(figsize=(8,7))
  ax = fig.add_subplot(1,1,1)
  cax = ax.scatter(tab['ra'],tab['dec'],
                    c=tab['peak_flux']/tab['local_rms'],
                    cmap=cmap, norm=norm, zorder=100)
  for r in tab:
      ax.add_patch(Ellipse((r['ra'],r['dec']),
                            width=0.5, height=3, angle=r['epoch']/nepochs*360,
                            alpha=0.35, edgecolor='none',
                            color=cmap(norm(r['peak_flux']/r['local_rms'])),
                            zorder=norm(r['peak_flux']/r['local_rms'])
                          ))

  cb = fig.colorbar(cax,ax=ax)
  cb.set_ticks(range(3,16,2))
  cb.set_label("SNR", **kwargs)
  ax.set_xlabel("RA J2000",**kwargs)
  ax.set_ylabel("Dec J2000", **kwargs)
  # flip the x axis so that RA increases to the left
  ax.set_xlim((ax.get_xlim()[1],ax.get_xlim()[0]))
  ax.grid()
  pyplot.savefig("transients.png")
  """
}

process reproject_images {
  label 'python'
  publishDir params.output_dir, mode: 'copy'

  input:
  tuple val(basename), path(mean_img)
  path fits

  output:
  path 'reprojected_images'

  script:
  """
  echo ${task.process} on \${HOSTNAME}
  mkdir reprojected_images
  ls *.fits | grep -v "mean_image" > temp_epochs.txt
  ls *mean_image.* > temp_mean.txt
  reprojection.py --epochs temp_epochs.txt --mean temp_mean.txt --reproj_dir reprojected_images
  """

}

workflow {
  get_version( )
  // image_ch = epoch_label, image_fits
  if (params.convolve) {
    convolve_beams(image_ch.map{it->it[1]}.collect())
    image_ch = convolve_beams.out.flatten().map{it -> tuple(it.baseName.split('_')[0], it)}
    //image_ch.view()
  }
  bane_raw( image_ch )
  // image_bkg_rms = epoch_label, image_fits, [bkg_fits, rms_fits]
  image_bkg_rms = bane_raw.out
  if ( params.fits_warp ) {
    if ( params.ref_catalogue == null ) {
      // No ref catalogue supplied so download default one
      download_gleam_catalogue()
      ref_cat = download_gleam_catalogue.out
    }
    else {
      ref_cat = Channel.fromPath( params.ref_catalogue )
    }
    initial_sfind( image_bkg_rms )
    fits_warp(
      initial_sfind.out,
      ref_cat,
    )
    image_ch = fits_warp.out
    image_bkg_rms = fits_warp.out.concat(bane_raw.out.map{ it -> [it[0], [it[1..2]]]}).groupTuple().map{ it -> [ it[0], it[1][0], it[1][1][0][1]]}
  }

  make_mean_image( image_ch.map{ it -> it[1] }.collect() )
  make_sky_coverage( image_ch.map{ it -> it[1] }.collect() )
  bane_mean_image( make_mean_image.out[0] )
  sfind_mean_image( bane_mean_image.out )
  source_monitor(
    sfind_mean_image.out,
    image_bkg_rms,
  )
  join_fluxes(
    source_monitor.out[1].collect(),
    sfind_mean_image.out,
  )
  compute_stats( join_fluxes.out )
  // plot_lc( compute_stats.out )
  mask_images(
    sfind_mean_image.out,
    image_bkg_rms,
  )
  sfind_masked( mask_images.out )
  compile_transients_candidates( sfind_masked.out.collect() )
  transients_plot( compile_transients_candidates.out )
  reproject_images(make_mean_image.out)
}