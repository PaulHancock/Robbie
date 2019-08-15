.SECONDARY:
.ONESHELL:
SHELL:=/bin/bash

###
# Setup and configuration parameters
###

# input images should be listed in epoch order in this file
IMFILE:=sim_images.txt
IMAGES:=$(shell cat $(IMFILE))
NEPOCH:=$(shell cat $(IMFILE) | wc -l)
# set warp to be empty to run astrometry corrections
WARP:=NO
# (external) reference catalogue used for astrometry correction via fits_warp
REFCAT:=Reference_sim.fits
REFCAT_RA:=ra
REFCAT_DEC:=dec
# This variable is used to invoke stilts
STILTS:=java -jar /home/hancock/Software/topcat/topcat-full.jar -stilts
# prefix for output files
PREFIX:=sim_
# The name of the mean image and image cube
MEAN:=$(PREFIX)mean.fits
CUBE:=$(PREFIX)cube.fits
# region file to use for masking the source finding and reference catalogue.
REGION:=square.mim

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
	@echo "  variables = variables.png"
	@echo "  sceince = variables + transients"
	@echo ""
	@echo 'I recommend that you `make science`'

cite:
	@echo '% If you make use of Robbie please cite the following work:'
	@echo '% Hancock et al. 2019:'
	@echo '@ARTICLE{2019A&C....27...23H,'
	@echo '   author = {{Hancock}, P.~J. and {Hurley-Walker}, N. and {White}, T.~E.},'
	@echo '    title = "{ROBBIE: A batch processing work-flow for the detection of radio transients and variables}",'
	@echo '  journal = {Astronomy and Computing},'
	@echo 'archivePrefix = "arXiv",'
	@echo '   eprint = {1902.06956},'
	@echo ' primaryClass = "astro-ph.IM",'
	@echo ' keywords = {Methods, Data analysis, Techniques, Radio astronomy, Variability, Transients},'
	@echo '     year = 2019,'
	@echo '    month = apr,'
	@echo '   volume = 27,'
	@echo '      eid = {23},'
	@echo '    pages = {23},'
	@echo '      doi = {10.1016/j.ascom.2019.02.004},'
	@echo '   adsurl = {http://adsabs.harvard.edu/abs/2019A%26C....27...23H},'
	@echo '  adsnote = {Provided by the SAO/NASA Astrophysics Data System}'
	@echo '}'

# Shorcuts for easy processing
variables: $(PREFIX)variables.png
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
	then rm -f $@; ln -s $$(basename $<) $@ ;\
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

$(PREFIX)flux_table.vot: $(IMAGES:.fits=_warped_prior_comp.fits)
	files=($^) ;\
	cmd="$(STILTS) tmatchn nin=$${#files[@]} matcher=exact out=$@ " ;\
	for n in $${!files[@]} ;\
	do \
	m=$$( echo "$${n}+1" | bc ) ;\
	cmd="$${cmd} in$${m}=$${files[$${n}]} values$${m}='uuid' suffix$${m}=_$${n}" ;\
	done ;\
	echo $${cmd} | bash

$(PREFIX)flux_table_var.vot: $(PREFIX)flux_table.vot
	ndof=($$(./auto_corr.py $(PREFIX)cube.fits)) ;\
	./calc_var.py --infile $< --outfile $@ --ndof $${ndof[-1]}
	./plot_lc.py $@

#$(PREFIX)variables.png: $(PREFIX)flux_table_var.fits
#	./plot_variables.py --in $< --plot $@

$(PREFIX)variables.png: $(PREFIX)flux_table_var.vot
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
	if [[ -e $^ ]];\
	  then ./filter_transients.py --incat $^ --image $*_warped_blanked.fits --outcat $@ ;\
	  nsrc=($$( $(STILTS) tcat omode=count in=$@ ));\
	  nsrc=$${nsrc[-1]};\
	  if [[ $${nsrc} -lt 1 ]];\
	    then rm $@;\
	  fi;\
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
