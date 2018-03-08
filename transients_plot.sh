#!/usr/bin/env bash


for freq in 154 185
do

java -jar ~/Software/topcat/topcat-full.jar -stilts plot2plane \
   xpix=645 ypix=563 \
   xflip=true xlabel=RAJ2000 ylabel=DEJ2000 grid=true texttype=antialias \
    fontsize=14 fontstyle=serif fontweight=bold \
   xmin=315 xmax=360 ymin=-35 ymax=5 \
   auxmap=sron auxquant=12 auxmin=3 auxmax=15 \
   auxvisible=true auxlabel=peak_flux/local_rms \
   legend=false \
   layer=Size \
      in=${freq}MHz_transients.fits \
      x=ra y=dec size=epoch+2 aux=peak_flux/local_rms \
      shading=aux shape=open_circle scale=1.5 autoscale=false \
    out=${freq}MHz_transients.png

done