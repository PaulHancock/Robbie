.PHONY: help all clean input trimmed split bkg cubes medians priors stats xmatch push_warp pull_warp submit_jobs
.SECONDARY:

help:
	echo "Use clean, input,.... "

all: input rmbad trimmed split bkg cubes medians sfind xmatch push_warp

clean:
	rm *.fits *.dat k2.mim k2.reg

input:
	rsync --ignore-existing --progress hancock@mwa-process02.ivec.org:'~stingay/MWA-SkyMapper/*.fits' .

rmbad: bad_files.dat
	./kill_bad.sh

trimmed:
	./trim_all.sh

split: K2_154MHz.dat K2_185MHz.dat

bkg:
	./bkg_all.sh

cubes: cube_154MHz.fits cube_185MHz.fits

medians: median_154MHz.fits median_185MHz.fits

sfind: median_154MHz_comp.fits median_185MHz_comp.fits
	./sfind_all.sh

priors: median_154MHz_comp.fits median_185MHz_comp.fits
	./priorize_all.sh

stats: 154MHz_flux_table.fits 185MHz_flux_table.fits

joined_154MHz.csv joined_185MHz.csv: K2_154MHz.dat K2_185MHz.dat
	./join_catalogues.sh

K2_trim_%.fits: K2_final_%.fits
	getfits -o $@ -x 4800 4800 $< 3000 3400

K2_trim_%_bkg.fits: K2_trim_%.fits
	BANE $<

K2_trim_%_warped_comp.fits: K2_trim_%_warped.fits K2_trim_%_bkg.fits K2_trim_%_rms.fits k2.mim
	aegean $< --background K2_trim_$*_bkg.fits --noise K2_trim_$*_rms.fits --table $< --island --region k2.mim

K2_trim_%_comp.fits: K2_trim_%.fits K2_trim_%_bkg.fits k2.mim
	aegean $< --autoload --table $< --island --region k2.mim

k2.mim k2.reg: 
	MIMAS +c 337.5 -14.5 18 -o k2.mim
	MIMAS --mim2reg k2.mim k2.reg

K2_154MHz.dat K2_185MHz.dat:
	./split_freqs.sh

cube_154MHz.fits cube_185MHz.fits: K2_154MHz.dat K2_185MHz.dat
	python make_cube.py

median_154MHz.fits median_185MHz.fits: cube_154MHz.fits cube_185MHz.fits
	python make_median.py

median_%MHz_bkg.fits: median_%MHz.fits
	BANE $<

median_%MHz_comp.fits: median_%MHz.fits median_%MHz_bkg.fits k2.mim
	aegean $< --autoload --table $< --island --region k2.mim

154MHz_flux_table.fits 185MHz_flux_table.fits:
	./construct_flux_table.sh

154MHz_flux_table_var.fits: 154MHz_flux_table.fits
	python calc_var.py

xmatch:
	python correct_astrometry.py xmatch 154
	python correct_astrometry.py xmatch 185

push_warp:
	./push2galaxy.sh

pull_warp:
	./pullFromGalaxy.sh

submit_jobs:
	./run_jobs.sh

makefile2dot.py:
	wget https://github.com/vak/makefile2dot/raw/master/makefile2dot.py

vis.png: Makefile makefile2dot.py
	python makefile2dot.py < $< | dot -Tpng > $@
