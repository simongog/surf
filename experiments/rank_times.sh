#!/bin/bash
CUR_DIR=`pwd`
MY_DIR="$( cd "$( dirname "$0" )" && pwd )" # gets the directory where the script is located in
cd "${MY_DIR}"
MY_DIR=`pwd`
#SURF_PATH=$MY_DIR/..
SURF_PATH=/scratch/VR0052/ESA2014/surf

COLLECTIONS="$SURF_PATH/collections/gov2"
EXP_DIR="$SURF_PATH/experiments"
PORT=12345

INDEXES="IDX_DR IDX_D IDX_D1R1 INVIDX_E"
RANKERS="BM25 TFIDF LMDS"

echo "qryid;collection;ranker;index;qrymode;k;qrylen;res_size;qry_time;search_time;nodes_evaluated;nodes_total;postings_evaluated;postings_total;client_time" > $EXP_DIR/nodes_evaluated.csv

for col in $COLLECTIONS
do
    for idx in $INDEXES
    do
        for rank in $RANKERS
        do
            $SURF_PATH/build/surf_daemon-${idx}_$rank -c $col -p $PORT &
            for k in 100
            do
                $SURF_PATH/build/surf_query -h localhost:$PORT -q $SURF_PATH/queries/trec2005-efficiency-1000.qry -k $k -r 1 -a -p >> $EXP_DIR/ranker_times_2005.csv
                $SURF_PATH/build/surf_query -h localhost:$PORT -q $SURF_PATH/queries/trec2006-efficiency-1000.qry -k $k -r 1 -a -p >> $EXP_DIR/ranker_times_2006.csv
		    done
            # shut down daemon
            $SURF_PATH/build/surf_query -h localhost:$PORT -q $SURF_PATH/queries/wiki.q -k 1 -s > /dev/null
        done
    done
done

cd "${CUR_DIR}"
