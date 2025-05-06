#!/bin/bash

TRACE_DIR="/mnt/shared_traces/total_text_traces"
SEGMENT_DIR="$TRACE_DIR/segments"
RESULT_DIR="../result"
PROGRAM_BASE="../build"
MODE=${1:-release}  # default to 'release', use 'debug' for debug mode

THREAD_COUNTS=(1 2 4 6 8 10 12)

# Set binary path
if [[ "$MODE" == "debug" ]]; then
    PROGRAM="gdb --args $PROGRAM_BASE/debug/microbench"
else
    PROGRAM="$PROGRAM_BASE/release/microbench"
fi

mkdir -p "$RESULT_DIR"

for SEGMENT_FILE in "$SEGMENT_DIR"/optimized_segments_*.csv; do
    # Extract trace base name from segment file
    echo "seg file: ${SEGMENT_FILE}"
    BASENAME=$(basename "$SEGMENT_FILE")
    echo "base file: ${BASENAME}"
    TRACE_NAME="${BASENAME#optimized_segments_}"
    echo "trace file: ${TRACE_NAME}"
    TRACE_NAME="${TRACE_NAME%.csv}"
    echo "Final trace file: ${TRACE_NAME}"

    # Try to find the corresponding trace file (first match)
    # Try to find a trace file that matches exactly (with or without .csv suffix)
    TRACE_FILE=$(find "$TRACE_DIR" -maxdepth 1 -type f | grep -E "/${TRACE_NAME}(\.csv)?$" | head -n 1)


    if [[ -z "$TRACE_FILE" ]]; then
        echo "No matching trace found for config $SEGMENT_FILE"
        continue
    fi
    echo "Trace: ${TRACE_FILE}"
    echo "Segment: ${SEGMENT_FILE}"
    echo "---------------------------------------------------------------------------------"
    
    LOG_FILE="$RESULT_DIR/log_${TRACE_NAME}.log"
    for THREAD_NUM in "${THREAD_COUNTS[@]}"; do
        echo "Running trace [$TRACE_FILE] with $THREAD_NUM threads"

        if [[ "$THREAD_NUM" -eq 1 ]]; then
            INDEX_LIST="alex,alexol,lipp,lippol,btreeolc,pgm,libox,xindex"
        else
            INDEX_LIST="alexol,lippol,btreeolc,libox,xindex"
        fi

        OUTPUT_FILE="$RESULT_DIR/out_readonly.csv"

        {
            echo "===== Test: $TRACE_NAME | Threads: $THREAD_NUM ====="
            echo "Start time: $(date '+%Y-%m-%d %H:%M:%S')"
            START_TIME=$(date +%s)

            $PROGRAM \
            --keys_file="$TRACE_FILE" \
            --keys_file_type=text \
            --config_file="$SEGMENT_FILE" \
            --read=1.0 \
            --insert=0.0 \
            --operations_num=200000000 \
            --output_path="$OUTPUT_FILE" \
            --table_size=-1 \
            --init_table_ratio=1 \
            --thread_num="$THREAD_NUM" \
            --index="$INDEX_LIST"

            END_TIME=$(date +%s)
            echo "End time: $(date '+%Y-%m-%d %H:%M:%S')"
            echo "Duration: $((END_TIME - START_TIME)) seconds"
            echo "==============================================="
        } 2>&1 | tee "$LOG_FILE"

    done
done
