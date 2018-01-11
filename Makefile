.PHONY: input
.SECONDARY:

input:
	rsync --progress hancock@mwa-process02.ivec.org:'~stingay/MWA-SkyMapper/*.fits' .

K2_trim_%.fits: K2_final_%.fits
	getfits -o $@ -x 4800 4800 $< 3000 3400

K2_trim_%_bkg.fits: K2_trim_%.fits
	BANE $<

K2_trim_%_comp.fits: K2_trim_%.fits K2_trim_%_bkg.fits k2.mim
	aegean $< --autoload --table $< --island --region k2.mim

k2.mim k2.reg: 
	MIMAS +c 337.5 -14.5 18 -o k2.mim
	MIMAS --mim2reg k2.mim k2.reg

K2_154MHz.dat K2_185MHz.dat:
	./split_freqs.sh

cube_%MHz.fits: K2_%MHz.dat
	python make_cube.py

median_%MHz.fits: cube_%MHz.fits
	python make_median.py

median_%MHz_bkg.fits: median_%MHz.fits
	BANE $<

median_%MHz_comp.fits: median_%MHz.fits median_%MHz_bkg.fits k2.mim
	aegean $< --autoload --table $< --island --region k2.mim