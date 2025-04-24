HISTOGRAM=../histogram/histogram
GUESTS=CH2CL2__ CHCL3___ CH3CL___ CH3BR___ BR2_____ BR2+20C_ BR2+10C_ BR2+20L_ BR2-20L_ BR2+10S_ BR2+10SE BR2-10S_ BR2+10E_ BR2-10E_ CJTHF___ NEGTHF__ CPENTANE LJME____ CPENTAN- CPENTAN+ CPENTA-- CPENTA++ C2H4F2__ 	C2H5F___ 	C2HF5___ 	CCL2F2__	CCL3F___	CH2CLF__	CH2FCF3_	CH3CL___	CHCL2F__	CHCLF2__	CL2FCCH3	CLF2CCH3	PROPANE_	ISOBUTAN
test: LJME____.12hedra.histo LJME____.14hedra.histo LJME____.16hedra.histo
%.12hedra.histo: $(HISTOGRAM)
	cat DEFR 12hedra.nx4a | time $(HISTOGRAM) -i $* > $@
%.14hedra.histo: $(HISTOGRAM)
	cat DEFR 14hedra.nx4a | time $(HISTOGRAM) -i $* > $@
%.16hedra.histo: $(HISTOGRAM)
	cat DEFR 16hedra.nx4a | time $(HISTOGRAM) -i $* > $@
# %.fvalues: vdwp/physconst.py vdwp/chempot.py fvalue.py %.12hedra.histo %.14hedra.histo %.16hedra.histo
# 	python fvalue.py $* | tee $@

all: $(patsubst %, %.12hedra.histo, $(GUESTS)) $(patsubst %, %.14hedra.histo, $(GUESTS)) $(patsubst %, %.16hedra.histo, $(GUESTS))

distclean:
	-rm *.fvalues *.histo
