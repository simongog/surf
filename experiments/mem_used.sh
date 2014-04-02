#!/bin/bash
CUR_DIR=`pwd`
MY_DIR="$( cd "$( dirname "$0" )" && pwd )" # gets the directory where the script is located in
cd "${MY_DIR}"
MY_DIR=`pwd`
SURF_PATH=/scratch/VR0052/ESA2014/surf

COLLECTIONS="$SURF_PATH/collections/gov2 $SURF_PATH/collections/cluewebB"
EXP_DIR="$SURF_PATH/experiments"
PORT=12345

INDEXES="IDX_D IDX_DR"

for col in $COLLECTIONS
do
    echo $col
    for idx in $INDEXES
    do
	echo $idx
        $SURF_PATH/build_turpin/surf_index-$idx -c $col 
    done
done

cd "${CUR_DIR}"
