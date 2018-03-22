.SECONDARY:
.ONESHELL:
SHELL:=/bin/bash

IMFILE:=all_images.txt
IMAGES:=$(shell cat $(IMFILE))

help:
	echo "help!"

clean:
	rm *.fits *.dat k2.mim k2.reg

test:
	f=(${IMAGES}) ;\
	for a in $${f[@]};\
	do \
	echo "$${a}" ;\
	done

# dummy rules to indicate that these files are pre-existing
$(IMAGES) region.mim:

GLEAM_SUB.fits: region.mim
	MIMAS --maskcat $< ~/alpha/DATA/GLEAM_EGC.fits $@ --colnames RAJ2000 DEJ2000 --negate

# Background and noise maps for the sub images
$(IMAGES:.fits=_bkg.fits): %_bkg.fits : %.fits
	 BANE $<
# This duplication is required since I can't mash it into the above rule whislt also using target:pattern:prereq
# However it causes BANE to run twice so we first check to see if the file is already built. arg!
$(IMAGES:.fits=_rms.fits): %_rms.fits : %.fits
	test -f $@ || echo BANE again $<

# Blind source finding
$(IMAGES:.fits=_comp.fits): %_comp.fits : %.fits %_bkg.fits %_rms.fits region.mim
	aegean $< --autoload --island --table $<,$*.reg --region=region.mim


# cross matching
$(IMAGES:.fits=_xm.fits): %_xm.fits : %_comp.fits GLEAM_SUB.fits
	./correct_astrometry.py GLEAM_SUB.fits $< $@

# warping
$(IMAGES:.fits=_warped.fits): %_warped.fits : %.fits %_xm.fits
	./fits_warp.py --infits $< --xm $*_xm.fits --suffix warped --ra1 ra --dec1 dec --ra2 RAJ2000 --dec2 DEJ2000

# Create cube
cube.fits: $(IMAGES:.fits=_warped.fits)
	./make_cube.py $@ $^

# create mean image
mean.fits: cube.fits
	./make_mean.py $^ $@

mean_bkg.fits mean_rms.fits: mean.fits
	BANE $<

# create master (mean) catalogue
mean_comp.fits: mean.fits mean_bkg.fits mean_rms.fits
	aegean $< --autoload --island --table $<,$*.reg

# priorize to make light curves from warped images
$(IMAGES:.fits=_warped_prior_comp.fits): %_warped_prior_comp.fits : %_warped.fits %_bkg.fits %_rms.fits mean_comp.fits
	aegean ${image} --background $*_bkg.fits --noise $*_rms.fits \
                    --table $*_warped_prior.fits,$*_warped_prior.reg --priorized 2 \
		            --input mean_comp.fits --noregroup

# joine all priorized sources into a single table
flux_table.fits: $(IMAGES:.fits=_warped_prior_comp.fits)
	files=($^) ;\
	cmd="java -jar /home/hancock/Software/stilts.jar tmatchn nin=$${#files[@]} matcher=exact out=$@ " ;\
	for n in $${!files[@]} ;\
	do ;\
	m=$$( echo "$${n}+1" | bc ) ;\
	cmd="$${cmd} in$${m}=$${files[${n}]} values$${m}='uuid' suffix$${m}=_$${n}" ;\
	done ;\
	echo $${cmd} | bash

# add variability stats to the flux table
flux_table_var.fits: flux_table.fits
	./calc_var.py $< $@


# blank the warped images
$(IMAGES:.fits=_warped_blanked.fits): %_blanked.fits : %.fits %_comp.fits
	AeRes -c $*_comp.fits -f $< -r $@ --mask --sigma=0.1

# blind source find warped/blanked images
$(IMAGES:.fits=_warped_blanked_comp.fits): %_warped_blanked_comp.fits : %_warped_blanked.fits $_rms.fits $_bkg.fits region.mim
	aegean $< --background $*_bkg.fits --noise $*_rms.fits \
	--table $*_warped_blanked.fits,$*_warped_blanked.reg --island --region region.mim


# remove the bad transients
$(IMAGE:.fits=_warped_blanked_comp_filtered.fits): %_warped_blanked_comp_filtered.fits: %_warped_blanked_comp.fits
	./filter_transients.py $^ $@

# join all transients into one catalogue
transients.fits: $(IMAGE:.fits=_warped_blanked_comp_filtered.fits)
	files=($^) ;\
	cmd="java -jar /home/hancock/Software/stilts.jar tcatn nin=$${#files[@]}" ;\
	for i in $$( seq 1 1 $${#files[@]} ) ;\
	do ;\
	j=$$( echo "$${i} -1" | bc ) ;\
	cmd="$${cmd} in$${i}=$${files[$${j}]} icmd$${i}='addcol epoch $${i}'" ;\
	done ;\
	cmd="$${cmd} out=$@ ofmt=fits" ;\
	$$($${cmd})

# plot the transients into a single image
transients.png: transients.fits
 	java -jar ~/Software/topcat/topcat-full.jar -stilts plot2plane \
   xpix=645 ypix=563 \
   xflip=true xlabel=RAJ2000 ylabel=DEJ2000 grid=true texttype=antialias \
    fontsize=14 fontstyle=serif fontweight=bold \
   auxmap=sron auxquant=12 auxmin=3 auxmax=15 \
   auxvisible=true auxlabel=peak_flux/local_rms \
   legend=false \
   layer=Size \
      in=$< \
      x=ra y=dec size=epoch+2 aux=peak_flux/local_rms \
      shading=aux shape=open_circle scale=1.5 autoscale=false \
    out=$@

makefile2dot.py:
	wget https://github.com/vak/makefile2dot/raw/master/makefile2dot.py

vis.png: Makefile makefile2dot.py
	python makefile2dot.py < $< | dot -Tpng > $@
