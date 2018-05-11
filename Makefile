.SECONDARY:
.ONESHELL:
SHELL:=/bin/bash
# input images
IMFILE:=all_images.txt
IMAGES:=$(shell cat $(IMFILE))
# reference catalogue
REFCAT:=/home/hancock/alpha/DATA/GLEAM_EGC.fits
# invocation of stilts
STILTS:=java -jar /home/hancock/Software/topcat/topcat-full.jar -stilts
# prefix for outputfiles
PREFIX:=
CUBE:=$(PREFIX)cube.fits
MEAN:=$(PREFIX)mean.fits
# region file to use
REGION:=region.mim


help:
	echo "help!"

variables: $(PREFIX)flux_table_var.fits
transients: $(PREFIX)transients.png
science: variables transients

# dummy rules to indicate that these files are pre-existing
$(IMAGES) $(REGION):

$(PREFIX)GLEAM_SUB.fits: $(REGION)
	MIMAS --maskcat $< $(REFCAT) $@ --colnames RAJ2000 DEJ2000 --negate

# Background and noise maps for the sub images
$(IMAGES:.fits=_bkg.fits): %_bkg.fits : %.fits
	 BANE $<
# This duplication is required since I can't mash it into the above rule whislt also using target:pattern:prereq
# However it causes BANE to run twice so we first check to see if the file is already built. arg!
$(IMAGES:.fits=_rms.fits): %_rms.fits : %.fits
	test -f $@ || BANE $<

# Blind source finding
$(IMAGES:.fits=_comp.fits): %_comp.fits : %.fits %_bkg.fits %_rms.fits $(REGION)
	aegean $< --autoload --island --table $<,$*.reg --region=$(REGION)


# cross matching
$(IMAGES:.fits=_xm.fits): %_xm.fits : %_comp.fits $(PREFIX)GLEAM_SUB.fits
	./fits_warp.py --refcat $(PREFIX)GLEAM_SUB.fits --incat $< --xm $@

# warping
$(IMAGES:.fits=_warped.fits): %_warped.fits : %.fits %_xm.fits
	./fits_warp.py --infits $< --xm $*_xm.fits --suffix warped --ra1 ra --dec1 dec --ra2 RAJ2000 --dec2 DEJ2000 --plot

# Create cube
$(CUBE): $(IMAGES:.fits=_warped.fits)
	./make_cube.py $@ $^

# create mean image
$(MEAN): $(CUBE)
	./make_mean.py $^ $@

$(MEAN:.fits=_bkg.fits) $(MEAN:.fits=_rms.fits): $(MEAN)
	BANE $<

# create master (mean) catalogue
$(MEAN:.fits=_comp.fits): %_comp.fits : %.fits %_bkg.fits %_rms.fits
	aegean $< --autoload --island --table $*.fits,$*.reg --region=$(REGION)

# priorize to make light curves from warped images
$(IMAGES:.fits=_warped_prior_comp.fits): %_warped_prior_comp.fits : %_warped.fits %_bkg.fits %_rms.fits $(PREFIX)mean_comp.fits
	aegean $< --background $*_bkg.fits --noise $*_rms.fits \
                    --table $*_warped_prior.fits,$*_warped_prior.reg --priorized 2 \
		            --input $(PREFIX)mean_comp.fits --noregroup

# join all priorized sources into a single table
$(PREFIX)flux_table.fits: $(IMAGES:.fits=_warped_prior_comp.fits)
	files=($^) ;\
	cmd="java -jar /home/hancock/Software/stilts.jar tmatchn nin=$${#files[@]} matcher=exact out=$@ " ;\
	for n in $${!files[@]} ;\
	do \
	m=$$( echo "$${n}+1" | bc ) ;\
	cmd="$${cmd} in$${m}=$${files[$${n}]} values$${m}='uuid' suffix$${m}=_$${n}" ;\
	done ;\
	echo $${cmd} | bash

# add variability stats to the flux table
$(PREFIX)flux_table_var.fits: $(PREFIX)flux_table.fits
	./calc_var.py $< $@
	./plot_lc.py $@


# blank the warped images
$(IMAGES:.fits=_warped_blanked.fits): %_warped_blanked.fits : %_warped.fits $(PREFIX)mean_comp.fits
	AeRes -c $(PREFIX)mean_comp.fits -f $< -r $@ --mask --sigma=0.1

# blind source find warped/blanked images
$(IMAGES:.fits=_warped_blanked_comp.fits): %_warped_blanked_comp.fits : %_warped_blanked.fits %_rms.fits %_bkg.fits $(REGION)
	aegean $< --background $*_bkg.fits --noise $*_rms.fits \
	--table $*_warped_blanked.fits,$*_warped_blanked.reg --island --region $(REGION)


# remove the bad transients
$(IMAGES:.fits=_warped_blanked_comp_filtered.fits): %_warped_blanked_comp_filtered.fits : %_warped_blanked_comp.fits
	./filter_transients.py $^ $*_warped_blanked.fits $@

# join all transients into one catalogue
$(PREFIX)transients.fits: $(IMAGES:.fits=_warped_blanked_comp_filtered.fits)
	files=(${^}) ;\
	cmd="${STILTS} tcatn nin=$${#files[@]}" ;\
	for i in $$( seq 1 1 $${#files[@]} ) ;\
	do \
	j=$$( echo "$${i} -1" | bc ) ;\
	cmd="$${cmd} in$${i}=$${files[$${j}]} icmd$${i}='addcol epoch $${i}'" ;\
	done ;\
	cmd="$${cmd} out=$@ ofmt=fits" ;\
	echo $${cmd} | bash

# plot the transients into a single image
$(PREFIX)transients.png: $(PREFIX)transients.fits
	$(STILTS) plot2plane \
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