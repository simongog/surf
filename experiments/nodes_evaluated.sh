#!/bin/bash
CUR_DIR=`pwd`
MY_DIR="$( cd "$( dirname "$0" )" && pwd )" # gets the directory where the script is located in
cd "${MY_DIR}"
MY_DIR=`pwd`
SURF_PATH="$MY_DIR/.."

COLLECTIONS="/devhome3/sgog/ESA2014/surf/collections/gov2/"
PORT=12345

INDEXES="IDX_D IDX_D_SANSLEN IDX_DR IDX_DR_SANSLEN"

for col in $COLLECTIONS
do
    for idx in $INDEXES
    do
        echo "$SURF_PATH/build/surf_daemon_$idx -c $col -p $PORT &"
        for k in 10 100 1000 
        do
            echo "$SURF_PATH/build/surf_query -h localhost:$PORT -q $SURF_PATH/queries/trec2005-efficiency-100.qry -k $k -r 1 -p >> nodes_evaluated.csv"
		done
        # shut down daemon
        echo "$SURF_PATH/build/surf_query -h localhost:$PORT -q $SURF_PATH/queries/wiki.q -k 1 -s > /dev/null"
    done
done

cd "${CUR_DIR}"
