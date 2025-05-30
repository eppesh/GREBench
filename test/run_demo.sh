#!/bin/bash

# Set binary path based on mode
PROGRAM="../build/release/microbench"
# Set the thread number here
THREAD_NUM=64

# Determine index list based on thread number
if [[ "$THREAD_NUM" -eq 1 ]]; then
    INDEX_LIST="libox"
else
    INDEX_LIST="alexol,btreeolc,libox,xindex"
fi

if [[ "$1" == "debug" ]]; then
    PROGRAM="gdb --args ../build/debug/microbench"
    INDEX_LIST="libox"
fi


# Run the benchmark
${PROGRAM} \
  --keys_file=/home/shuaihua/traces/libox/text_format/osm.csv \
  --keys_file_type=text \
  --config_file=/home/shuaihua/traces/libox/segments/segments_osm.csv \
  --read=0.0 \
  --insert=0.05 \
  --operations_num=1000000 \
  --scan=0.95 --scan_num=100 \
  --output_path=../result/out_demo.csv \
  --table_size=-1 \
  --init_table_ratio=0.5 \
  --thread_num=${THREAD_NUM} \
  --index=${INDEX_LIST}
