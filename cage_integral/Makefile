HISTOGRAM=../histogram/histogram

test: LJME____.12hedra.histo LJME____.14hedra.histo LJME____.16hedra.histo
%.12hedra.histo: $(HISTOGRAM)
	cat DEFR 12hedra.nx4a | time $(HISTOGRAM) -i $* > $@
%.14hedra.histo: $(HISTOGRAM)
	cat DEFR 14hedra.nx4a | time $(HISTOGRAM) -i $* > $@
%.16hedra.histo: $(HISTOGRAM)
	cat DEFR 16hedra.nx4a | time $(HISTOGRAM) -i $* > $@
# %.fvalues: vdwp/physconst.py vdwp/chempot.py fvalue.py %.12hedra.histo %.14hedra.histo %.16hedra.histo
# 	python fvalue.py $* | tee $@

distclean:
	-rm *.fvalues *.histo
