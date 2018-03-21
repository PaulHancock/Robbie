.PHONY: help all clean input trimmed split bkg cubes medians priors stats xmatch push_warp pull_warp submit_jobs
.SECONDARY:

IMFILE:=all_images.txt
IMAGES:=$( shell cat $(IMFILE))

help:
	echo "help!"

clean:
	rm *.fits *.dat k2.mim k2.reg

# dummy rules to indicate that these files are pre-existing
$(IMAGES):

# Background and noise maps for the sub images
$(IMAGES:.fits=_bkg.fits): %_bkg.fits : %.fits
	 BANE $<
# This duplication is required since I can't mash it into the above rule whislt also using target:pattern:prereq
# However it causes BANE to run twice so we first check to see if the file is already built. arg!
$(IMAGES:.fits=_rms.fits): %_rms.fits : %.fits
	test -f $@ || echo BANE again $<

# Blind source finding
$(IMAGES:.fits=_comp.fits): %_comp.fits : %.fits %_bkg.fits %_rms.fits
	aegean $< --autoload --island --table $<,$*.reg


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

# blank the warped images
$(IMAGES:.fits=_warped_blanked.fits): %_blanked.fits : %.fits %_comp.fits
	AeRes -c $*_comp.fits -f $< -r $@ --mask --sigma=0.1

# source find warped/blanked images
$(IMAGES:.fits=_warped_blanked_comp.fits): %_warped_blanked_comp.fits : %_warped_blanked.fits $_rms.fits $_bkg.fits
	aegean $< --background $*_bkg.fits --noise $*_rms.fits \
	--table $*_warped_blanked.fits,$*_warped_blanked.reg --island



stats: 154MHz_flux_table_var.fits 185MHz_flux_table_var.fits

transients:
	./filter_transents.sh


154MHz_flux_table.fits 185MHz_flux_table.fits:
	./construct_flux_table.sh

154MHz_flux_table_var.fits 185MHz_flux_table_var.fits: 154MHz_flux_table.fits 185MHz_flux_table.fits
	python calc_var.py


makefile2dot.py:
	wget https://github.com/vak/makefile2dot/raw/master/makefile2dot.py

vis.png: Makefile makefile2dot.py
	python makefile2dot.py < $< | dot -Tpng > $@

gleam: k2.mim
	MIMAS --maskcat k2.mim ~/alpha/DATA/GLEAM_EGC.fits GLEAM_SUB.fits --colnames RAJ2000 DEJ2000 --negate

%MHz_transients.fits:
	./collect_transients.sh

%MHz_transients.png: %MHz_transients.fits
	./transients_plot.sh