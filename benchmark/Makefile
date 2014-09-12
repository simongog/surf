include ../external/sdsl-lite/Make.helper

CXX_FLAGS = $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(INC_DIR) -L$(LIB_DIR)  
LIBS = -lsdsl -ldivsufsort -ldivsufsort64

EXP1IDX:=$(call config_column,experiment1.config,1)
EXP1COL:=$(call config_column,experiment1.config,2)
EXP1LEN:=$(call config_column,experiment1.config,3)
EXP1K:=$(call config_column,experiment1.K.config,1)
EXP1RES:=$(foreach K,$(EXP1K),results/experiment1.$(K).txt)

results/experiment1.txt: ../build/surf_query-$(EXP1IDX) $(EXP1RES) 
	@cat $(EXP1RES) > $@

results/experiment1.%.txt: pattern/$(EXP1COL).$(EXP1LEN).txt
	$(eval ONLYMATCH:=$(call config_select,experiment1.K.config,$*,2)) 
	../build/surf_query-$(EXP1IDX) -c ../collections/$(EXP1COL) \
		                -q pattern/$(EXP1COL).$(EXP1LEN) -k $* \
						-m $(ONLYMATCH) > $@
	../build/surf_query-$(EXP1IDX) -c ../collections/$(EXP1COL) \
		                -q pattern/$(EXP1COL).$(EXP1LEN) -k $* \
						$(ONLYMATCH) >> $@
	

#RESULT_FILE=results/all.txt
#RESULT_FILE_INT=results/all_int.txt
#
#COLS:=$(call config_column,collections.config,1)
#IDXS:=$(call config_column,index.config,1)
#LENS:=$(call config_column,pattern_len.config,1)
#KS:=$(call config_column,K.config,1)
#
#PATS:=$(foreach COL,$(COLS),\
#    	$(foreach LEN,$(LENS),pattern/$(COL).$(LEN)))
#RESS:=$(foreach COL,$(COLS),\
#    	$(foreach LEN,$(LENS),\
#			$(foreach K,$(KS),\
#				$(foreach IDX,$(IDXS),results/$(COL).$(LEN).$(K).$(IDX)))))
#
#
#experiment1: $(PATS) $(RESS)
#	cat $(RESS) > results/experiment1.txt
#
#pattern/[COL].[LEN]
pattern/%: ../build/gen_patterns
	$(eval COL:=$(call dim,1,$*)) 
	$(eval LEN:=$(call dim,2,$*)) 
	@../build/gen_patterns -c ../collections/$(COL) -m $(LEN) -x 4000 > $@

##pattern/[COL].[LEN].[IDX]
#results/%:
#	$(eval COL:=$(call dim,1,$*)) 
#	$(eval LEN:=$(call dim,2,$*)) 
#	$(eval K:=$(call dim,3,$*)) 
#	$(eval IDX:=$(call dim,4,$*)) 
#	../build/surf_query-$(IDX) -c ../collections/$(COL) -q pattern/$(COL).$(LEN) -k $(K) -m > $@
