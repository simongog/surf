#!/bin/bash

COLLECTION="../collections/wikishort/"

RANKERS="BM25 LMDS TFIDF"

INDEXES="IDX_D IDX_DR IDX_D1R1 INVIDX_E"

SURF_PATH="../"
QRYBIN=$SURF_PATH/build/surf_query
QRYFILE="wikishort.qry"

for col in $COLLECTION
do
    for rank in $RANKERS
    do
        OUTPUT_FILES_OR=""
        OUTPUT_FILES_AND=""
        for idx in $INDEXES
        do
            DAEMONBIN="$SURF_PATH/build/surf_daemon-${idx}_$rank"
            $DAEMONBIN -c $col > /dev/null 2>&1 &
            $QRYBIN -q $QRYFILE -k 10 -r 1 -R 1> OUT_OR_${idx}_$rank 2>/dev/null
            $QRYBIN -q $QRYFILE -k 10 -r 1 -R -a 2>/dev/null 1> OUT_AND_${idx}_$rank 
            OUTPUT_FILES_OR="$OUTPUT_FILES_OR OUT_OR_${idx}_$rank"
            OUTPUT_FILES_AND="$OUTPUT_FILES_AND OUT_AND_${idx}_$rank"
            $QRYBIN -q $QRYFILE -s > /dev/null 2>&1 1> /dev/null
        done

        # cmp for equality now...
        diff $OUTPUT_FILES_OR > /dev/null 2>&1
        if [ $? -eq 1 ]
        then
            echo "output for OR NOT EQUAL for $rank $col !!!"
        else
            echo "all good for OR and $rank"
        fi
        diff $OUTPUT_FILES_AND > /dev/null 2>&1
        if [ $? -eq 1 ]
        then
            echo "output for AND NOT EQUAL for $rank $col !!!"
        else
            echo "all good for AND and $rank"
        fi
    done
done


# cleanup
#rm -f OUT_*
