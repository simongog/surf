#!/bin/bash
CUR_DIR=`pwd`
MY_DIR="$( cd "$( dirname "$0" )" && pwd )" # gets the directory where the script is located in
cd "${MY_DIR}"
MY_DIR=`pwd`
SURF_PATH=$MY_DIR/..
#SURF_PATH=/scratch/VR0052/ESA2014/surf

#COLLECTIONS="$SURF_PATH/collections/gov2 $SURF_PATH/collections/cluewebB"
COLLECTIONS="/devhome3/sgog/ESA2014/surf/collections/gov2"
EXP_DIR="$SURF_PATH/experiments"
PORT=12345

INDEXES="IDX_D"

echo "qryid;collection;index;qrymode;k;qrylen;res_size;qry_time;search_time;nodes_evaluated;nodes_total;postings_evaluated;postings_total;client_time" > $EXP_DIR/phrase_time_2005-2.csv

for col in $COLLECTIONS
do
    for idx in $INDEXES
    do
        $SURF_PATH/build/surf_daemon-$idx -c $col -p $PORT &
        for k in 10 100
        do
            $SURF_PATH/build/surf_query -h localhost:$PORT -q $SURF_PATH/queries/trec2006-efficiency-1000.qry -k $k -r 1 -P 10 -p >> $EXP_DIR/phrase_time_2006-2.csv
		done
        # shut down daemon
        $SURF_PATH/build_turpin/surf_query -h localhost:$PORT -q $SURF_PATH/queries/wiki.q -k 1 -s > /dev/null
    done
done

cd "${CUR_DIR}"
