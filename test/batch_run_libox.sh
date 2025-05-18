#!/bin/bash

TRACE_DIR="/mnt/shared_traces/libox/text_format"
SEGMENT_DIR="/mnt/shared_traces/libox/segments/longitudes"
RESULT_DIR="../result/libox"
mkdir -p "$RESULT_DIR" "$RESULT_DIR/logs"

OUTPUT_FILE="$RESULT_DIR/out_libox_all_0516.csv"  # Single CSV file for all results

PROGRAM_BASE="../build"
MODE=${1:-release}  # Default mode is 'release'

THREAD_COUNTS=(1 2 4 6 8 10 12)

# Set binary path
if [[ "$MODE" == "debug" ]]; then
    PROGRAM="gdb --args $PROGRAM_BASE/debug/microbench"
else
    PROGRAM="$PROGRAM_BASE/release/microbench"
fi

# Define workload configurations
declare -A WORKLOADS
WORKLOADS[readonly]="--read=1.0 --insert=0.0 --init_table_ratio=1 --operations_num=200000000"
#WORKLOADS[readintensive]="--read=0.8 --insert=0.2 --init_table_ratio=0.5 --operations_num=200000000"
WORKLOADS[balanced]="--read=0.5 --insert=0.5 --init_table_ratio=0.5 --operations_num=200000000"
#WORKLOADS[writeintensive]="--read=0.2 --insert=0.8 --init_table_ratio=0.5 --operations_num=200000000"
WORKLOADS[writeonly]="--read=0.0 --insert=1.0 --init_table_ratio=0.5 --operations_num=200000000"

# Loop through segment config files
for SEGMENT_FILE in "$SEGMENT_DIR"/segments_*.csv; do
    BASENAME=$(basename "$SEGMENT_FILE")
    TRACE_NAME="${BASENAME#segments_}"
    TRACE_NAME="${TRACE_NAME%.csv}"

    # Find corresponding trace file
    TRACE_FILE=$(find "$TRACE_DIR" -maxdepth 1 -type f | grep -E "/${TRACE_NAME}(\.csv|\.txt)?$" | head -n 1)
    if [[ -z "$TRACE_FILE" ]]; then
        echo "No matching trace found for config $SEGMENT_FILE"
        continue
    fi

    echo "Using Trace: $TRACE_FILE"
    echo "Segment: $SEGMENT_FILE"
    echo "--------------------------------------------------------------------------------"

    LOG_FILE="$RESULT_DIR/logs/log_${TRACE_NAME}_0516.log"
    for WORKLOAD_NAME in "${!WORKLOADS[@]}"; do
        WORKLOAD_PARAMS="${WORKLOADS[$WORKLOAD_NAME]}"

        for THREAD_NUM in "${THREAD_COUNTS[@]}"; do
            echo "Running [$WORKLOAD_NAME] with $THREAD_NUM threads on [$TRACE_NAME]"

            if [[ "$THREAD_NUM" -eq 1 ]]; then
                INDEX_LIST="alexol,lippol,btreeolc,xindex,artolc,libox"
            else
                INDEX_LIST="alexol,lippol,btreeolc,xindex,artolc,libox"
            fi

            {
                echo "===== Test: $TRACE_NAME | Workload: $WORKLOAD_NAME | Threads: $THREAD_NUM ====="
                echo "Start time: $(date '+%Y-%m-%d %H:%M:%S')"
                START_TIME=$(date +%s)

                $PROGRAM \
                --keys_file="$TRACE_FILE" \
                --keys_file_type=text \
                --config_file="$SEGMENT_FILE" \
                $WORKLOAD_PARAMS \
                --output_path="$OUTPUT_FILE" \
                --table_size=-1 \
                --memory \
                --thread_num="$THREAD_NUM" \
                --index="$INDEX_LIST"

                END_TIME=$(date +%s)
                echo "End time: $(date '+%Y-%m-%d %H:%M:%S')"
                echo "Duration: $((END_TIME - START_TIME)) seconds"
                echo "==================================================================================="
            } 2>&1 | tee -a "$LOG_FILE"
        done
    done
done
