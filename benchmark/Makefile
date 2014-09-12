include ../external/sdsl-lite/Make.helper

CXX_FLAGS = $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(INC_DIR) -L$(LIB_DIR)  
LIBS = -lsdsl -ldivsufsort -ldivsufsort64

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
	../build/gen_patterns -c ../collections/$(COL) -m $(LEN) -x 4000 > $@

