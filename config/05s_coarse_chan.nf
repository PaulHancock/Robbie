#! /usr/bin/env nextflow

log.info """\
         RUBY
         ====
         output to    : ${params.output_dir}
         --
         run as       : ${workflow.commandLine}
         config files : ${workflow.configFiles}
         container    : ${workflow.containerEngine}:${workflow.container}
         """
         .stripIndent()

channel_ch = Channel.from(0..23)
             .map{it -> String.format("%04d",it)}
             .view()

obsid_ch1 = Channel.create()
obsid_ch2 = Channel.create()
trange_ch1 = Channel.create()
trange_ch2 = Channel.create()

Channel.from(1217495184,
             1217495304,
             1217495424,
             1217495544,
             1217495664,
             1217495784,
             1217495904,
             1217496024,
             1217496144,
             1217496264,
             1217496384,
             1217496504,
             1217496624,
             1217496744,
             1217496864)
       .into(obsid_ch1, obsid_ch2)

Channel.from(9..228)
       .map{t -> tuple(t, t+1)}
       .into(trange_ch1, trange_ch2)

// create a set which is one xm_ch per group of image_ch
group_ch = obsid_ch1.combine(trange_ch1).combine(channel_ch)
           .map{it -> tuple(
                   file("${params.xm_base}regrid_${it[0]}-0.5s-t${it[1]}-${it[2]}-pbcorr-I_xm.fits"),
                   file("${params.image_base}${it[0]}/${it[0]}-0.5s-24c-t${it[1]}-${it[2]}-${it[3]}-pbcorr-I.fits")
                )}        /* creates every xm/chan pair */
           .groupTuple()  /* groups all the chan together in xm groups */


/* WORKFLOW DESCRIPTION */
process fits_warp {
  label 'warp'

  input:
  tuple(path(xm), path('*')) from group_ch

  output:
  path('*_warped.fits') into done_ch

  script:
  """
  echo ${task.process} on \${HOSTNAME}
  ls *.fits
  set -x
  for f in \$( ls *I.fits )
  do 
    fits_warp.py --infits \${f} --xm ${xm} --suffix warped --cores ${task.cpus}
  done
  ls *.fits
  """
}