#!/bin/bash

# Set binary path based on mode
PROGRAM="../build/release/microbench"
# Set the thread number here
THREAD_NUM=1

# Determine index list based on thread number
if [[ "$THREAD_NUM" -eq 1 ]]; then
    # INDEX_LIST="alex,lipp,hope,lisa"
    INDEX_LIST="lisa,hope,alex,lipp"
else
    INDEX_LIST="alexol,btreeolc,libox,xindex"
fi

if [[ "$1" == "debug" ]]; then
    PROGRAM="gdb --args ../build/debug/microbench"
    INDEX_LIST="hope"
fi


# Run the benchmark
${PROGRAM} \
  --keys_file=/mnt/shared_traces/hope/w048.csv \
  --keys_file_type=text \
  --sample_distribution=uniform \
  --config_file=w087 \
  --read=1 \
  --insert=0 \
  --operations_num=10000000 \
  --output_path=../result/out_demo.csv \
  --table_size=-1 \
  --init_table_ratio=1 \
  --node_capacity=100 \
  --top_k=0.05 \
  --density_factor=4 \
  --temp_node_cap=20 \
  --max_error_rs=32 \
  --min_line_len=5 \
  --thread_num=${THREAD_NUM} \
  --index=${INDEX_LIST}
