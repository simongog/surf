include ../external/sdsl-lite/Make.helper

CXX_FLAGS = $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(INC_DIR) -L$(LIB_DIR)  
LIBS = -lsdsl -ldivsufsort -ldivsufsort64

is_int_col=$(shell if [[ "$1" == *INT || "$1" == "gov2" || "$1" == "trec8" ]]; then echo "_int"; else echo ""; fi)
is_int_idx=$(shell if [[ "$1" == *_INT ]]; then echo "1"; else echo "0"; fi)

# Experiment 1
# Fix collection, pattern length, and index
# Vary k, 
EXP1IDX:=$(call config_column,experiment1.config,1)
EXP1COL:=$(call config_column,experiment1.config,2)
EXP1LEN:=$(call config_column,experiment1.config,3)
EXP1K:=$(call config_column,experiment1.K.config,1)
EXP1RES:=$(foreach K,$(EXP1K),results/experiment1.$(K).txt)

experiment1: results/experiment1.txt

experiment2: results/experiment2.txt

experiment3: results/experiment3.txt

results/experiment1.txt: exp1.txt ../build/surf_query-$(EXP1IDX) $(EXP1RES) 
	@cat $(EXP1RES) > $@

results/experiment1.%.txt: pattern/$(EXP1COL).$(EXP1LEN).byte.txt
	$(eval ONLYMATCH:=$(call config_select,experiment1.K.config,$*,2)) 
	../build/surf_query-$(EXP1IDX) -c ../collections/$(EXP1COL) \
		                -q pattern/$(EXP1COL).$(EXP1LEN).byte.txt -k $* \
						-m $(ONLYMATCH) > $@
	../build/surf_query-$(EXP1IDX) -c ../collections/$(EXP1COL) \
		                -q pattern/$(EXP1COL).$(EXP1LEN).byte.txt -k $* \
						$(ONLYMATCH) >> $@

exp1.txt: ../build/surf_index-$(EXP1IDX) ../collections/$(EXP1COL)
	../build/surf_index-$(EXP1IDX) -c ../collections/$(EXP1COL) > $@ 
	
pattern/%.byte.txt: ../build/gen_patterns
	$(eval COL:=$(call dim,1,$*)) 
	$(eval LEN:=$(call dim,2,$*)) 
	../build/gen_patterns$(SUF) -c ../collections/$(COL) -m $(LEN) -x 40000 > $@

pattern/%.int.txt: ../build/gen_patterns_int
	$(eval COL:=$(call dim,1,$*)) 
	$(eval LEN:=$(call dim,2,$*)) 
	../build/gen_patterns_int -c ../collections/$(COL) -m $(LEN) -x 40000 > $@



# Experiment 2
# Fix k, pattern length, 
# Vary index (IDX_NN, IDX_NNX), collections
EXP2IDX:=$(call config_column,experiment2.IDX.config,1)
EXP2IDXINT:=$(call config_column,experiment2.IDX_INT.config,1)
EXP2COL:=$(call config_column,experiment2.COL.config,1)
EXP2COLINT:=$(call config_column,experiment2.COL_INT.config,1)
EXP2K:=$(call config_column,experiment2.config,1)
EXP2LEN:=$(call config_column,experiment2.config,2)
EXP2LENINT:=$(call config_column,experiment2.config,3)
EXP2RES:=$(foreach IDX,$(EXP2IDX),\
			$(foreach COL,$(EXP2COL),results/experiment2.$(IDX).$(COL).byte.txt))
EXP2RESINT:=$(foreach IDX,$(EXP2IDXINT),\
			$(foreach COL,$(EXP2COLINT),results/experiment2.$(IDX).$(COL).int.txt))

EXP2:=$(foreach IDX,$(EXP2IDX),\
			$(foreach COL,$(EXP2COL),exp2.build.$(IDX).$(COL).byte))
EXP2INT:=$(foreach IDX,$(EXP2IDXINT),\
			$(foreach COL,$(EXP2COLINT),exp2.build.$(IDX).$(COL).int))

EXP2PAT:=$(foreach COL,$(EXP2COL),pattern/$(COL).$(EXP2LEN).byte.txt)
EXP2PATINT:=$(foreach COL,$(EXP2COLINT),pattern/$(COL).$(EXP2LENINT).int.txt)

results/experiment2.%.byte.txt: 
	$(eval IDX:=$(call dim,1,$*))
	$(eval COL:=$(call dim,2,$*))
	../build/surf_query-$(IDX) -c ../collections/$(COL) \
		-q pattern/$(COL).$(EXP2LEN).byte.txt -k $(EXP2K) \
		> $@

results/experiment2.%.int.txt: 
	$(eval IDX:=$(call dim,1,$*))
	$(eval COL:=$(call dim,2,$*))
	../build/surf_query-$(IDX) -c ../collections/$(COL) \
		-q pattern/$(COL).$(EXP2LENINT).int.txt -k $(EXP2K) \
		> $@

results/experiment2.txt: exp2.txt $(EXP2PAT) $(EXP2RES) $(EXP2PATINT) $(EXP2RESINT)
	@cat $(EXP2RES) $(EXP2RESINT) > $@

exp2.txt: $(EXP2) $(EXP2INT)

exp2.build.%.byte:
	rm -f exp2.txt
	$(eval IDX:=$(call dim,1,$*))
	$(eval COL:=$(call dim,2,$*))
	../build/surf_index-$(IDX) -c ../collections/$(COL) > exp2.txt

exp2.build.%.int:
	rm -f exp2.txt
	$(eval IDX:=$(call dim,1,$*))
	$(eval COL:=$(call dim,2,$*))
	../build/surf_index-$(IDX) -c ../collections/$(COL) > exp2.txt


# Experiment 3
# As in Figure 2+3 in Plug&Play paper.
# Fix k=10, index=IDX_NNX
# Vary index, pattern length, collection

EXP3IDX:=$(call config_column,experiment3.config,1)
EXP3IDXINT:=$(call config_column,experiment3.config,2)
EXP3K:=$(call config_column,experiment3.config,3)
EXP3LEN:=$(call config_column,experiment3.PAT.config,1)
EXP3COL:=$(call config_column,experiment3.COL.config,1)
EXP3COLINT:=$(call config_column,experiment3.COL_INT.config,1)
EXP3RES:=$(foreach COL,$(EXP3COL),\
			$(foreach LEN,$(EXP3LEN),\
				results/experiment3.$(EXP3IDX).$(COL).$(LEN).byte.txt))
EXP3RESINT:=$(foreach COL,$(EXP3COLINT),\
			$(foreach LEN,$(EXP3LEN),\
				results/experiment3.$(EXP3IDXINT).$(COL).$(LEN).int.txt))
		
EXP3:=$(foreach COL,$(EXP3COL),exp3.build.$(EXP3IDX).$(COL).byte)
EXP3INT:=$(foreach COL,$(EXP3COLINT),exp3.build.$(EXP3IDXINT).$(COL).int)
EXP3PAT:=$(foreach COL,$(EXP3COL),\
			$(foreach LEN,$(EXP3LEN),pattern/$(COL).$(LEN).byte.txt))
EXP3PATINT:=$(foreach COL,$(EXP3COLINT),\
			$(foreach LEN,$(EXP3LEN),pattern/$(COL).$(LEN).int.txt))

results/experiment3.%.txt:
	$(eval IDX:=$(call dim,1,$*))
	$(eval COL:=$(call dim,2,$*))
	$(eval LEN:=$(call dim,3,$*))
	$(eval ALPH:=$(call dim,4,$*))
	../build/surf_query-$(IDX) -c ../collections/$(COL) \
		-q pattern/$(COL).$(LEN).$(ALPH).txt -k $(EXP3K) \
		> $@

results/experiment3.txt: exp3.txt $(EXP3PAT) $(EXP3RES) $(EXP3PATINT) $(EXP3RESINT)
	@cat $(EXP3RES) $(EXP3RESINT) > $@

exp3.txt: $(EXP3) $(EXP3INT)

exp3.build.%.byte:
	rm -f exp3.txt
	$(eval IDX:=$(call dim,1,$*))
	$(eval COL:=$(call dim,2,$*))
	../build/surf_index-$(IDX) -c ../collections/$(COL) > exp3.txt

exp3.build.%.int:
	rm -f exp3.txt
	$(eval IDX:=$(call dim,1,$*))
	$(eval COL:=$(call dim,2,$*))
	../build/surf_index-$(IDX) -c ../collections/$(COL) > exp3.txt
