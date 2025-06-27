#!/bin/bash

# Set binary path based on mode
PROGRAM="../build/release/microbench"
# Set the thread number here
THREAD_NUM=1

# Determine index list based on thread number
if [[ "$THREAD_NUM" -eq 1 ]]; then
    INDEX_LIST="hope,alex,lipp,lisa"
else
    INDEX_LIST="alexol,btreeolc,libox,xindex"
fi

if [[ "$1" == "debug" ]]; then
    PROGRAM="gdb --args ../build/debug/microbench"
    INDEX_LIST="hope"
fi


# Run the benchmark
${PROGRAM} \
  --keys_file=/mnt/shared_traces/libox/text_format/w048_sort_unique.csv \
  --keys_file_type=text \
  --config_file=/home/shuaihua/traces/libox/segments/segments_osm.csv \
  --read=1 \
  --insert=0 \
  --operations_num=1000000 \
  --output_path=../result/out_demo.csv \
  --table_size=-1 \
  --init_table_ratio=1 \
  --node_capacity=1000 \
  --thread_num=${THREAD_NUM} \
  --index=${INDEX_LIST}
