test.csv:
	python system2.py -T 223,303,3 -L BR2+10C_,BR2+10E_,BR2+10SE,BR2+10S_,BR2+20C_,BR2+20L_,BR2-10E_,BR2-10S_,BR2-20L_,BR2_____,C2H4F2__,C2H5F___,C2HF5___,CCL2F2__,CCL3F___,CH2CL2__,CH2CLF__,CH2FCF3_,CH3BR___,CH3CL___,CHCL2F__,CHCL3___,CHCLF2__,CJTHF___,CL2FCCH3,CLF2CCH3,CPENTA++,CPENTA--,CPENTAN+,CPENTAN-,CPENTANE,ISOBUTAN,LJAR____,LJME____,LJXE____,NEGTHF__,PROPANE_ > $@
#EPM2CO2_
test2.csv: Makefile
	python system2.py -T 183,353,3 -L C2H4F2__,C2H5F___,C2HF5___,CCL3F___,CH2CL2__,CH2CLF__,CH2FCF3_,CH3BR___,CH3CL___,CHCL2F__,CHCL3___,CHCLF2__,CL2FCCH3,CLF2CCH3,ISOBUTAN,PROPANE_,PROPANHT,LJPR____,LJAR____,LJME____,LJXE____,NEGTHF__,BR2_____,BR2-10E_,BR2-10S_,BR2-20L_,BR2+10C_,BR2+10E_,BR2+10S_,BR2+10SE,BR2+20C_,BR2+20L_,CJTHF___,CPENTA--,CPENTA++,CPENTAN-,CPENTAN+,CPENTANE > $@
#,EPM2CO2_
test3.csv: #Makefile
	basename -s .16hedra.histo `ls data/*.16*` > $@.mol
	python system3.py -T 183,353,3 -L `tr '\n' ',' < $@.mol` > $@
#,EPM2CO2_
%.labels: %.csv Makefile
	awk '{print "set label \""$$1"-"$$2"\" at "$$3","$$4" center offset 0,-0.4 font \"Helvetica,12\""}' $< > $@
