.SECONDARY:
.ONESHELL:
SHELL:=/bin/bash

###
# Setup and configuration parameters
###

# input images should be listed in epoch order in this file
IMFILE:=all_images.txt
IMAGES:=$(shell cat $(IMFILE))
NEPOCH:=$(shell cat $(IMFILE) | wc -l)
# set warp to be empty to avoid running astrometry corrections
WARP:=
# (external) reference catalogue used for astrometry correction via fits_warp
REFCAT:=/home/hancock/alpha/DATA/GLEAM_EGC.fits
REFCAT_RA:=RAJ2000
REFCAT_DEC:=DEJ2000
# This variable is used to invoke stilts
STILTS:=java -jar /home/hancock/Software/topcat/topcat-full.jar -stilts
# prefix for output files
PREFIX:=
# The name of the mean image and image cube
MEAN:=$(PREFIX)mean.fits
CUBE:=$(PREFIX)cube.fits
# region file to use for masking the source finding and reference catalogue.
REGION:=region.mim

# HELP!
help:
	@echo "Hello this is Robbie."
	@echo "Usage is: make [command | file]"
	@echo " files:"
	@echo "  refcat.fits - a masked version of the external reference catalogue"
	@echo "  cube.fits - a stack of astrometry corrected images"
	@echo "  mean.fits - a mean image from the above stack"
	@echo "  flux_table_var.fits - light curves and variability stats for all persistent sources"
	@echo "  transients.fits - a catalogue of all candidate transient events"
	@echo "  transients.png - a visualisation of transients.fits"
	@echo ""
	@echo " commands:"
	@echo "  transients = transients.png"
	@echo "  variables = flux_table_var.fits"
	@echo "  sceince = variables + transients"
	@echo ""
	@echo 'I recommend that you `make science`'

# Shorcuts for easy processing
variables: $(PREFIX)flux_table_var.fits
transients: $(PREFIX)transients.png
science: variables transients

# dummy rules to indicate that these files are pre-existing
$(IMAGES) $(REGION):

# Create a masked version of the reference catalogue
$(PREFIX)refcat.fits: $(REGION)
	if [[ -n "$(WARP)" ]] ;\
	then touch $@ ;\
	else \
	MIMAS --maskcat $< $(REFCAT) $@ --colnames $(REFCAT_RA) $(REFCAT_DEC) --negate ;\
	fi

###
# Apply astrometric correction to all the input images
###

# Background and noise maps for the sub images
$(IMAGES:.fits=_bkg.fits): %_bkg.fits : %.fits
	 BANE $<
# This duplication is required since I can't mash it into the above rule whislt also using target:pattern:prereq
# However it causes BANE to run twice so we first check to see if the file is already built. arg!
$(IMAGES:.fits=_rms.fits): %_rms.fits : %.fits
	test -f $@ || BANE $<

# Blind source finding on the input images
$(IMAGES:.fits=_comp.fits): %_comp.fits : %.fits %_bkg.fits %_rms.fits $(REGION)
	if [[ -n "$(WARP)" ]] ;\
	then touch $@ ;\
	else \
	aegean $< --autoload --island --table $<,$*.reg --region=$(REGION);\
	fi


# cross matching of blind catalogues with the reference catalogue, in order to create astrometry solutions
$(IMAGES:.fits=_xm.fits): %_xm.fits : %_comp.fits $(PREFIX)refcat.fits
	if [[ -n "$(WARP)" ]] ;\
	then touch $@ ;\
	else ./fits_warp.py --refcat $(PREFIX)refcat.fits --incat $< --xm $@ ;\
	fi

# warp the input images using the astrometry solutions, and create warped versions of the files
$(IMAGES:.fits=_warped.fits): %_warped.fits : %.fits %_xm.fits
	if [[ -n "$(WARP)" ]] ;\
	then rm $@; ln -s $< $@ ;\
	else \
	./fits_warp.py --infits $< --xm $*_xm.fits --suffix warped --ra1 ra --dec1 dec --ra2 $(REFCAT_RA) --dec2 $(REFCAT_DEC) --plot ;\
	fi

###
# Create a master catalogue of persistent sources from a mean of the warped images
###

# create an image cube from all the warped images
$(CUBE): $(IMAGES:.fits=_warped.fits)
	./make_cube.py $@ $^

# create mean image from the image cube
$(MEAN): $(CUBE)
	./make_mean.py $^ $@

# create background and noise mapse for the mean image
$(MEAN:.fits=_bkg.fits) $(MEAN:.fits=_rms.fits): $(MEAN)
	BANE $<

# create a catalogue from the mean image which will be used as a master catalogue
$(MEAN:.fits=_comp.fits): %_comp.fits : %.fits %_bkg.fits %_rms.fits
	aegean $< --autoload --island --table $*.fits,$*.reg --region=$(REGION)

###
# Persistent sources
###

# priorize source finding to make light curves from warped images based on the master catalogue
$(IMAGES:.fits=_warped_prior_comp.fits): %_warped_prior_comp.fits : %_warped.fits %_bkg.fits %_rms.fits $(PREFIX)mean_comp.fits
	aegean $< --background $*_bkg.fits --noise $*_rms.fits \
                    --table $*_warped_prior.fits,$*_warped_prior.reg --priorized 2 \
		            --input $(PREFIX)mean_comp.fits --noregroup

# join all priorized sources into a single table based on the UUID column
$(PREFIX)flux_table.fits: $(IMAGES:.fits=_warped_prior_comp.fits)
	files=($^) ;\
	cmd="$(STILTS) tmatchn nin=$${#files[@]} matcher=exact out=$@ " ;\
	for n in $${!files[@]} ;\
	do \
	m=$$( echo "$${n}+1" | bc ) ;\
	cmd="$${cmd} in$${m}=$${files[$${n}]} values$${m}='uuid' suffix$${m}=_$${n}" ;\
	done ;\
	echo $${cmd} | bash

# add variability stats to the flux table
$(PREFIX)flux_table_var.fits: $(PREFIX)flux_table.fits
	ndof=($$(./auto_corr.py $(PREFIX)cube.fits)) ;\
	./calc_var.py --infile $< --outfile $@ --ndof $(ndof[-1])
	./plot_lc.py $@

$(PREFIX)variables.png: $(PREFIX)flux_table_var.fits
	./plot_variables.py --in $< --plot $@

###
# Transient candidates
###

# blank the warped images
$(IMAGES:.fits=_warped_blanked.fits): %_warped_blanked.fits : %_warped.fits $(PREFIX)mean_comp.fits
	AeRes -c $(PREFIX)mean_comp.fits -f $< -r $@ --mask --sigma=0.1

# blind source find warped/blanked images
$(IMAGES:.fits=_warped_blanked_comp.fits): %_warped_blanked_comp.fits : %_warped_blanked.fits %_rms.fits %_bkg.fits $(REGION)
	aegean $< --background $*_bkg.fits --noise $*_rms.fits \
	--table $*_warped_blanked.fits,$*_warped_blanked.reg --island --region $(REGION)


# remove the bad transients
$(IMAGES:.fits=_warped_blanked_comp_filtered.fits): %_warped_blanked_comp_filtered.fits : %_warped_blanked_comp.fits
	files=$$( ls $^ ) ;\
	if [[ ! -z $${files} ]];\
	then ./filter_transients.py $${files} $*_warped_blanked.fits $@ ;\
	fi

# join all transients into one catalogue
$(PREFIX)transients.fits: $(IMAGES:.fits=_warped_blanked_comp_filtered.fits)
	files=($$( ls $^ )) ;\
	if [[ -z $${files} ]];\
	then touch $@;\
	else \
	cmd="$(STILTS) tcatn nin=$${#files[@]}" ;\
	for i in $$( seq 1 1 $${#files[@]} ) ;\
	do \
	j=$$( echo "$${i} -1" | bc ) ;\
	cmd="$${cmd} in$${i}=$${files[$${j}]} icmd$${i}='addcol epoch $${i}'" ;\
	done ;\
	cmd="$${cmd} out=$@ ofmt=fits" ;\
	echo $${cmd} | bash ;\
	fi

# plot the transients into a single image
$(PREFIX)transients.png: $(PREFIX)transients.fits
	./plot_transients.py --in $< --plot $@
