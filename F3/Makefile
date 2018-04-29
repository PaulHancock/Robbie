.SECONDARY:
.ONESHELL:
SHELL:=/bin/bash

help:
	echo "help "
input: trimmed K2_154MHz.dat K2_185MHz.dat

trimmed:
	for f in $$( ls K2_final_*.fits | sed 's:final:trim:g' ); \
	do \
	make $${f} ;\
	done


K2_trim_%.fits: K2_final_%.fits
	getfits -o $@ -x 4800 4800 $< 3000 3400

k2.mim k2.reg: 
	MIMAS +c 337.5 -14.5 18 -o k2.mim
	MIMAS --mim2reg k2.mim k2.reg

freqs.dat: 
	for f in $$(ls K2_trim_[0-9]*[0-9].fits); \
	do \
		echo $${f} `gethead $${f} CRVAL3` >> freqs.dat ;\
	done

K2_154MHz.dat K2_185MHz.dat: freqs.dat
	grep ' 15' freqs.dat | awk '{print $$1}' > K2_154MHz.dat
	grep ' 18' freqs.dat | awk '{print $$1}' > K2_185MHz.dat
