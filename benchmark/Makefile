include ../external/sdsl-lite/Make.helper

CXX_FLAGS = $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(INC_DIR) -L$(LIB_DIR)  
LIBS = -lsdsl -ldivsufsort -ldivsufsort64

# Experiment 1
# Fix collection, pattern length, and index
# Vary k, 
EXP1IDX:=$(call config_column,experiment1.config,1)
EXP1COL:=$(call config_column,experiment1.config,2)
EXP1LEN:=$(call config_column,experiment1.config,3)
EXP1K:=$(call config_column,experiment1.K.config,1)
EXP1RES:=$(foreach K,$(EXP1K),results/experiment1.$(K).txt)

results/experiment1.txt: exp1.txt ../build/surf_query-$(EXP1IDX) $(EXP1RES) 
	@cat $(EXP1RES) > $@

results/experiment1.%.txt: pattern/$(EXP1COL).$(EXP1LEN).txt
	$(eval ONLYMATCH:=$(call config_select,experiment1.K.config,$*,2)) 
	../build/surf_query-$(EXP1IDX) -c ../collections/$(EXP1COL) \
		                -q pattern/$(EXP1COL).$(EXP1LEN).txt -k $* \
						-m $(ONLYMATCH) > $@
	../build/surf_query-$(EXP1IDX) -c ../collections/$(EXP1COL) \
		                -q pattern/$(EXP1COL).$(EXP1LEN).txt -k $* \
						$(ONLYMATCH) >> $@

exp1.txt: ../build/surf_index-$(EXP1IDX) ../collections/$(EXP1COL)
	../build/surf_index-$(EXP1IDX) -c ../collections/$(EXP1COL) > $@ 
	
pattern/%: ../build/gen_patterns
	$(eval COL:=$(call dim,1,$*)) 
	$(eval LEN:=$(call dim,2,$*)) 
	../build/gen_patterns -c ../collections/$(COL) -m $(LEN) -x 40000 > $@

# Experiment 2
# Fix k, pattern length, 
# Vary index (IDX_NN, IDX_NNX), collections
EXP2IDX:=$(call config_column,experiment2.IDX.config,1)
EXP2COL:=$(call config_column,experiment2.COL.config,1)
EXP2K:=$(call config_column,experiment2.config,1)
EXP2LEN:=$(call config_column,experiment2.config,2)
EXP2RES:=$(foreach IDX,$(EXP2IDX),\
			$(foreach COL,$(EXP2COL),results/experiment2.$(IDX).$(COL).txt))
EXP2:=$(foreach IDX,$(EXP2IDX),\
			$(foreach COL,$(EXP2COL),exp2.build.$(IDX).$(COL)))
EXP2PAT:=$(foreach COL,$(EXP2COL),pattern/$(COL).$(EXP2LEN).txt)

results/experiment2.%.txt: 
	$(eval IDX:=$(call dim,1,$*))
	$(eval COL:=$(call dim,2,$*))
	../build/surf_query-$(IDX) -c ../collections/$(COL) \
		-q pattern/$(COL).$(EXP2LEN).txt -k $(EXP2K) \
		> $@

results/experiment2.txt: exp2.txt $(EXP2PAT) $(EXP2RES)
	@cat $(EXP2RES) > $@

exp2.txt: $(EXP2)

exp2.build.%:
	rm exp2.txt
	$(eval IDX:=$(call dim,1,$*))
	$(eval COL:=$(call dim,2,$*))
	../build/surf_index-$(IDX) -c ../collections/$(COL) > exp2.txt


# Experiment 3
# As in Figure 2+3 in Plug&Play paper.
# Fix k=10
# Vary index (GREEDY, SORT, IDX_NNX), pattern length


