#!/bin/bash

# Set binary path based on mode
PROGRAM="../build/release/microbench"
if [[ "$1" == "debug" ]]; then
    PROGRAM="gdb --args ../build/debug/microbench"
fi

# Set the thread number here
THREAD_NUM=1

# Determine index list based on thread number
if [[ "$THREAD_NUM" -eq 1 ]]; then
    INDEX_LIST="alex,alexol,lipp,lippol,btreeolc,pgm,libox,xindex"
else
    INDEX_LIST="alexol,lippol,btreeolc,libox,xindex"
fi

# Run the benchmark
${PROGRAM} \
  --keys_file=/mnt/shared_traces/vmware_u32/text_format/w106.csv \
  --keys_file_type=text \
  --config_file=/mnt/shared_traces/vmware_u32/segments/optimized_segments_w106.csv \
  --read=1.0 \
  --insert=0.0 \
  --operations_num=200000000 \
  --output_path=../result/out_demo.csv \
  --table_size=-1 \
  --init_table_ratio=1 \
  --thread_num=${THREAD_NUM} \
  --index=${INDEX_LIST}
