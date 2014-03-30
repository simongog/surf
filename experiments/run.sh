#!/bin/bash
CUR_DIR=`pwd`
MY_DIR="$( cd "$( dirname "$0" )" && pwd )" # gets the directory where the script is located in
cd "${MY_DIR}"
MY_DIR=`pwd`
SURF_PATH="$MY_DIR/.."

COLLECTIONS="$SURF_PATH/collections/gov2/"
PORT=12345

for col in $COLLECTIONS
do
    for idx in $SURF_PATH/build/surf_daemon-IDX_SAWIT2
    do
        IDXNAME=$(echo $idx | sed 's/.*-\(.*\)/\1/g')
        $idx -c $col -p $PORT &
        for k in 10 100 1000 
        do
            $SURF_PATH/build/surf_query -h localhost:$PORT -q $SURF_PATH/queries/trec2005-efficiency-100.qry -k $k -r 1 -p >> trec-2005.csv
            $SURF_PATH/build/surf_query -h localhost:$PORT -q $SURF_PATH/queries/trec2006-efficiency-100.qry -k $k -r 1 -p >> trec-2006.csv
            $SURF_PATH/build/surf_query -h localhost:$PORT -q $SURF_PATH/queries/trec2005-efficiency-100.qry -k $k -r 1 -p -a >> trec-2005.csv
            $SURF_PATH/build/surf_query -h localhost:$PORT -q $SURF_PATH/queries/trec2006-efficiency-100.qry -k $k -r 1 -p -a >> trec-2006.csv
        done
        # shut down daemon
        $SURF_PATH/build/surf_query -h localhost:$PORT -q $SURF_PATH/queries/wiki.q -k 1 -s > /dev/null
    done
done

cd "${CUR_DIR}"
