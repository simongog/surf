#!/bin/bash
info=`find ../collections/*/index/space_*.html`

res=results/experiment0.txt

rm -f ${res}

for f in ${info}; do
    suf=${f#../collections/}
    col=${suf%/index/*}
    col_dir=${f%/index/*}
    col_size=`wc -c < ${col_dir}/text*_SURF.sdsl | tr -d ' '`
    size_info=`../build/size_info < $f`
    for line in ${size_info}; do
        echo "${line};${col};${col_size}" >> ${res}
    done
done
