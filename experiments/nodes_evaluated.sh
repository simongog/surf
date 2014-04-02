#!/bin/bash
CUR_DIR=`pwd`
MY_DIR="$( cd "$( dirname "$0" )" && pwd )" # gets the directory where the script is located in
cd "${MY_DIR}"
MY_DIR=`pwd`
SURF_PATH=$1

COLLECTIONS="$SURF_PATH/collections/gov2/"
PORT=12345

INDEXES="IDX_DR IDX_DR_SANSLEN IDX_D IDX_D_SANSLEN"

echo "qryid;collection;index;qrymode;k;qrylen;res_size;qry_time;search_time;nodes_evaluated;nodes_total;postings_evaluated;postings_total;client_time" > nodes_evaluated.csv

for col in $COLLECTIONS
do
    for idx in $INDEXES
    do
        $SURF_PATH/build/surf_daemon-$idx -c $col -p $PORT &
        for k in 10 100 1000 
        do
            $SURF_PATH/build/surf_query -h localhost:$PORT -q $SURF_PATH/queries/trec2005-efficiency-1000.qry -k $k -r 1 -p >> nodes_evaluated.csv
		done
        # shut down daemon
        $SURF_PATH/build/surf_query -h localhost:$PORT -q $SURF_PATH/queries/wiki.q -k 1 -s > /dev/null
    done
done

cd "${CUR_DIR}"
