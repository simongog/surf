
COLLECTIONS="/devhome3/sgog/ESA2014/surf/collections/gov2/"
SURF_PATH="/devhome3/sgog/ESA2014/surf"
PORT=12345

for col in $COLLECTIONS
do
    for idx in $SURF_PATH/build/surf_daemon-*
    do
        IDXNAME=$(echo $idx | sed 's/.*-\(.*\)/\1/g')
        $idx -c $col -p $PORT &
        for k in 10 100 1000 1000000000
        do
            $SURF_PATH/build/surf_query -h localhost:$PORT -q $SURF_PATH/queries/trec2005-efficiency-1000.qry -k $k >> trec-2005-time-$IDXNAME.csv
            $SURF_PATH/build/surf_query -h localhost:$PORT -q $SURF_PATH/queries/trec2006-efficiency-1000.qry -k $k >> trec-2006-time-$IDXNAME.csv
            $SURF_PATH/build/surf_query -h localhost:$PORT -q $SURF_PATH/queries/trec2005-efficiency-1000.qry -k $k -p >> trec-2005-profile-$IDXNAME.csv
            $SURF_PATH/build/surf_query -h localhost:$PORT -q $SURF_PATH/queries/trec2006-efficiency-1000.qry -k $k -p >> trec-2006-profile-$IDXNAME.csv
            $SURF_PATH/build/surf_query -h localhost:$PORT -q $SURF_PATH/queries/trec2005-efficiency-1000.qry -k $k -a >> trec-2005-and-time-$IDXNAME.csv
            $SURF_PATH/build/surf_query -h localhost:$PORT -q $SURF_PATH/queries/trec2006-efficiency-1000.qry -k $k -a >> trec-2006-and-time-$IDXNAME.csv
            $SURF_PATH/build/surf_query -h localhost:$PORT -q $SURF_PATH/queries/trec2005-efficiency-1000.qry -k $k -a -p >> trec-2005-and-profile-$IDXNAME.csv
            $SURF_PATH/build/surf_query -h localhost:$PORT -q $SURF_PATH/queries/trec2006-efficiency-1000.qry -k $k -a -p >> trec-2006-and-profile-$IDXNAME.csv
        done
        # shut down daemon
        $SURF_PATH/build/surf_query -h localhost:$PORT -q $SURF_PATH/queries/wiki.q -k 1 -s > /dev/null
    done
done
