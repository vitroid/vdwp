BINPATH=~/Dropbox/gitbox/CageIntegral
me:
	for mol in LJMEs2__ LJME____ LJMEs3__ LJMEs4__ LJMEs5__ LJMEe0__ LJMEe2__ \
	LJAR____ LJXE____ LJET____ LJBR2___ LJCS2___ LJKR____ LJCO2___ EKABR2__ \
	LJCF4___ \
	CH2CL2__ CHCL3___ CH3CL___ CH3BR___ C2H4F2__ 	C2H5F___ 	C2HF5___ \
		CCL2F2__	CCL3F___	CH2CLF__	CH2FCF3_	CH3CL___ \
		CHCL2F__	CHCLF2__	CL2FCCH3	CLF2CCH3	PROPANE_ \
		ISOBUTAN; do\
		for cage in 12 14 15 16; do\
			echo $$mol.ice$${cage}.histo; \
		done; \
	done | xargs make -j 8 -k 


all:
	for mol in TYTHF___ CH2CL2__ CHCL3___ CH3CL___ CH3BR___ BR2_____ \
	BR2+20C_ BR2+10C_ BR2+20L_ BR2-20L_ BR2+10S_ BR2+10SE BR2-10S_ \
	BR2+10E_ BR2-10E_ CJTHF___ NEGTHF__ CPENTANE LJME____ CPENTAN- \
	CPENTAN+ CPENTA-- CPENTA++ C2H4F2__ 	C2H5F___ 	C2HF5___ \
		CCL2F2__	CCL3F___	CH2CLF__	CH2FCF3_	CH3CL___ \
		CHCL2F__	CHCLF2__	CL2FCCH3	CLF2CCH3	PROPANE_ \
		ISOBUTAN	; do\
		for cage in 12 14 15 16; do\
			echo $$mol.$${cage}hedra.histo; \
		done; \
	done | xargs make -j 8 -k 

%.12hedra.histo: 12hedra.nx4a
	cat DEFR $< | time $(BINPATH)/histogram -i "$*" > $@
%.14hedra.histo: 14hedra.nx4a
	cat DEFR $< | time $(BINPATH)/histogram -i "$*" > $@
%.15hedra.histo: 15hedra.nx4a
	cat DEFR $< | time $(BINPATH)/histogram -i "$*" > $@
%.16hedra.histo: 16hedra.nx4a
	cat DEFR $< | time $(BINPATH)/histogram -i "$*" > $@

%.ice12.histo: 12hedra.nx4a
	cat DEFR $< | time $(BINPATH)/histogram -i "$*" -w "TIP4PICE" > $@
%.ice14.histo: 14hedra.nx4a
	cat DEFR $< | time $(BINPATH)/histogram -i "$*" -w "TIP4PICE"  > $@
%.ice15.histo: 15hedra.nx4a
	cat DEFR $< | time $(BINPATH)/histogram -i "$*" -w "TIP4PICE"  > $@
%.ice16.histo: 16hedra.nx4a
	cat DEFR $< | time $(BINPATH)/histogram -i "$*" -w "TIP4PICE"  > $@
