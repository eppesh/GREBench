#!/bin/bash
# ==============================================================================
# Flexible performance testing script for LiBox and other indexes
# 
# Usage examples:
#   ./run_tests.sh                       # Run all tests with defaults
#   ./run_tests.sh -t osm,fb             # Test only osm and fb traces
#   ./run_tests.sh -w readonly,balanced  # Test only readonly and balanced workloads
#   ./run_tests.sh -n 16,32,64           # Test only with 16, 32, and 64 threads
#   ./run_tests.sh -i alexol,loft        # Test only alexol and loft indexes
#   ./run_tests.sh -t osm -w writeonly -i libox,loft  # Test osm with writeonly workload using only libox and loft
#   ./run_tests.sh -o custom_results     # Output to custom filename
# ==============================================================================

# Default configuration
#TRACE_DIR="/home/shuaihua/traces/libox/alpha_beta"
#SEGMENT_DIR="/home/shuaihua/traces/libox/alpha_beta/segments"
TRACE_DIR="/mnt/shared_traces/hope"
SEGMENT_DIR="/home/shuaihua/traces/libox/segments"
# TRACE_DIR="/home/shuaihua/traces/libox/memory_test"
# SEGMENT_DIR="/home/shuaihua/traces/libox/memory_test/segments"
RESULT_DIR="../result/hope"
PROGRAM_BASE="../build"
DATE_TAG=$(date +"%m%d")
OUTPUT_FILE="$RESULT_DIR/out_hope_$DATE_TAG.csv"

# Available configurations 
#ALL_TRACES=(
    # "fb_10_10.csv"
    # "fb_20_10.csv"
    # "fb_30_10.csv"
    # "fb_40_10.csv"
    # "fb_50_5.csv"
    # "fb_50_10.csv"
    # "fb_50_20.csv"
    # "fb_50_30.csv"
    # "fb_50_40.csv"
    # "fb_50_50.csv"
    # "fb_60_10.csv"
    # "fb_70_10.csv"
    # "fb_80_10.csv"
    # "osm_10_10.csv"
    # "osm_20_10.csv"
    # "osm_30_10.csv"
    # "osm_40_10.csv"
    # "osm_50_5.csv"
    # "osm_50_10.csv"
    # "osm_50_20.csv"
    # "osm_50_30.csv"
    # "osm_50_40.csv"
    # "osm_50_50.csv"
    # "osm_60_10.csv"
    # "osm_70_10.csv"
    # "osm_80_10.csv"
    # "genome_10_10.csv"
    # "genome_20_10.csv"
    # "genome_30_10.csv"
    # "genome_40_10.csv"
    # "genome_50_5.csv"
    # "genome_50_10.csv"
    # "genome_50_20.csv"
    # "genome_50_30.csv"
    # "genome_50_40.csv"
    # "genome_50_50.csv"
    # "genome_60_10.csv"
    # "genome_70_10.csv"
    # "genome_80_10.csv"
    # "w048_10_10.csv"
    # "w048_20_10.csv"
    # "w048_30_10.csv"
    # "w048_40_10.csv"
    # "w048_50_5.csv"
    # "w048_50_10.csv"
    # "w048_50_20.csv"
    # "w048_50_30.csv"
    # "w048_50_40.csv"
    # "w048_50_50.csv"
    # "w048_60_10.csv"
    # "w048_70_10.csv"
    # "w048_80_10.csv"
    # "planet_10_10.csv"
    # "planet_20_10.csv"
    # "planet_30_10.csv"
    # "planet_40_10.csv"
    # "planet_50_5.csv"
    # "planet_50_10.csv"
    # "planet_50_20.csv"
    # "planet_50_30.csv"
    # "planet_50_40.csv"
    # "planet_50_50.csv"
    # "planet_60_10.csv"
    # "planet_70_10.csv"
    # "planet_80_10.csv"
    # "longitudes_10_10.csv"
    # "longitudes_20_10.csv"
    # "longitudes_30_10.csv"
    # "longitudes_40_10.csv"
    # "longitudes_50_5.csv"
    # "longitudes_50_10.csv"
    # "longitudes_50_20.csv"
    # "longitudes_50_30.csv"
    # "longitudes_50_40.csv"
    # "longitudes_50_50.csv"
    # "longitudes_60_10.csv"
    # "longitudes_70_10.csv"
    # "longitudes_80_10.csv"
 #)
ALL_TRACES=()
# Add w001.csv to w106.csv
# for i in $(seq -w 1 106); do
#     ALL_TRACES+=("w${i}.csv")
# done
# ALL_TRACES=("umass_financial1"  "umass_financial2"  "umass_websearch1"  "umass_websearch2"  "umass_websearch3")
#ALL_TRACES=("msr_proj.csv"  "msr_src1.csv"  "msr_src2.csv"  "msr_web.csv")
ALL_TRACES=("msr_proj.csv" "msr_src1.csv" "msr_src2.csv" "msr_web.csv" "umass_financial1.csv" "w025.csv" "w046.csv" "w048.csv" "w054.csv" "w056.csv" "w057.csv" "w087.csv")

# echo "${ALL_TRACES[@]}"  # Correct way to print all elements
# exit

# ALL_TRACES=(
#     "fb_100.csv"  "fb_60.csv"       "genome_20.csv"  "genome_80.csv"            "longitudes-200M_40.csv"  "msr_web_100.csv"  "msr_web_60.csv"  "osm_20.csv"  "osm_80.csv"    "w048_20.csv"  "w048_80.csv"
#     "fb_20.csv"   "fb_80.csv"       "genome_40.csv"  "longitudes-200M_100.csv"  "longitudes-200M_60.csv"  "msr_web_20.csv"   "msr_web_80.csv"  "osm_40.csv"  "w048_40.csv"
#     "fb_40.csv"   "genome_100.csv"  "genome_60.csv"  "longitudes-200M_20.csv"   "longitudes-200M_80.csv"  "msr_web_40.csv"   "osm_100.csv"     "osm_60.csv"  "w048_100.csv"  "w048_60.csv"
# )
# ALL_SEGMENTS=("segments_osm.csv" "segments_genome.csv" "segments_fb.csv" "segments_w027.csv" "segments_w045.csv" "segments_longitudes-200M.csv" "segments_msr_web.csv" "segments_w048.csv" "segments_umass_financial1.csv")
ALL_SEGMENTS=(
    "segments_osm.csv",
    "segments_genome.csv",
    "segments_fb.csv",
    "segments_longitudes-200M.csv",
    "segments_msr_web.csv",
    "segments_w048.csv"
)
# ALL_SEGMENTS=(
#     "segments_fb_100.csv"  "segments_fb_80.csv"       "segments_genome_60.csv"            "segments_longitudes-200M_40.csv"  "segments_msr_web_20.csv"  "segments_osm_100.csv"  "segments_osm_80.csv"    "segments_w048_60.csv"
#     "segments_fb_20.csv"   "segments_genome_100.csv"  "segments_genome_80.csv"            "segments_longitudes-200M_60.csv"  "segments_msr_web_40.csv"  "segments_osm_20.csv"   "segments_w048_100.csv"  "segments_w048_80.csv"
#     "segments_fb_40.csv"   "segments_genome_20.csv"   "segments_longitudes-200M_100.csv"  "segments_longitudes-200M_80.csv"  "segments_msr_web_60.csv"  "segments_osm_40.csv"   "segments_w048_20.csv"
#     "segments_fb_60.csv"   "segments_genome_40.csv"   "segments_longitudes-200M_20.csv"   "segments_msr_web_100.csv"         "segments_msr_web_80.csv"  "segments_osm_60.csv"   "segments_w048_40.csv"
# )
#ALL_SEGMENTS=(
    # "segments_fb_10_10.csv",
    # "segments_fb_20_10.csv",
    # "segments_fb_30_10.csv",
    # "segments_fb_40_10.csv",
    # "segments_fb_50_5.csv",
    # "segments_fb_50_10.csv",
    # "segments_fb_50_20.csv",
    # "segments_fb_50_30.csv",
    # "segments_fb_50_40.csv",
    # "segments_fb_50_50.csv",
    # "segments_fb_60_10.csv",
    # "segments_fb_70_10.csv",
    # "segments_fb_80_10.csv",    
    # "segments_osm_10_10.csv",
    # "segments_osm_20_10.csv",
    # "segments_osm_30_10.csv",
    # "segments_osm_40_10.csv",
    # "segments_osm_50_5.csv",
    # "segments_osm_50_10.csv",
    # "segments_osm_50_20.csv",
    # "segments_osm_50_30.csv",
    # "segments_osm_50_40.csv",
    # "segments_osm_50_50.csv",
    # "segments_osm_60_10.csv",
    # "segments_osm_70_10.csv",
    # "segments_osm_80_10.csv",
    # "segments_genome_10_10.csv",
    # "segments_genome_20_10.csv",
    # "segments_genome_30_10.csv",
    # "segments_genome_40_10.csv",
    # "segments_genome_50_5.csv",
    # "segments_genome_50_10.csv",
    # "segments_genome_50_20.csv",
    # "segments_genome_50_30.csv",
    # "segments_genome_50_40.csv",
    # "segments_genome_50_50.csv",
    # "segments_genome_60_10.csv",
    # "segments_genome_70_10.csv",
    # "segments_genome_80_10.csv",
    # "segments_w048_10_10.csv",
    # "segments_w048_20_10.csv",
    # "segments_w048_30_10.csv",
    # "segments_w048_40_10.csv",
    # "segments_w048_50_5.csv",
    # "segments_w048_50_10.csv",
    # "segments_w048_50_20.csv",
    # "segments_w048_50_30.csv",
    # "segments_w048_50_40.csv",
    # "segments_w048_50_50.csv",
    # "segments_w048_60_10.csv",
    # "segments_w048_70_10.csv",
    # "segments_w048_80_10.csv",
    # "segments_planet_10_10.csv",
    # "segments_planet_20_10.csv",
    # "segments_planet_30_10.csv",
    # "segments_planet_40_10.csv",
    # "segments_planet_50_5.csv",
    # "segments_planet_50_10.csv",
    # "segments_planet_50_20.csv",
    # "segments_planet_50_30.csv",
    # "segments_planet_50_40.csv",
    # "segments_planet_50_50.csv",
    # "segments_planet_60_10.csv",
    # "segments_planet_70_10.csv",
    # "segments_planet_80_10.csv",
    # "segments_longitudes_10_10.csv",
    # "segments_longitudes_20_10.csv",
    # "segments_longitudes_30_10.csv",
    # "segments_longitudes_40_10.csv",
    # "segments_longitudes_50_5.csv",
    # "segments_longitudes_50_10.csv",
    # "segments_longitudes_50_20.csv",
    # "segments_longitudes_50_30.csv",
    # "segments_longitudes_50_40.csv",
    # "segments_longitudes_50_50.csv",
    # "segments_longitudes_60_10.csv",
    # "segments_longitudes_70_10.csv",
    # "segments_longitudes_80_10.csv",
#)
ALL_WORKLOADS=("readonly" "balanced" "writeonly" "read20" "read40" "read60" "read80" "init20" "init40" "init60" "init80" "init100" "scanonly" "ycsbe" "scan10" "scan100" "scan1000" "scan10000" "scan100000" "scan500" "scan1500" "scan2000" "scan2500")
ALL_THREAD_COUNTS=(1 8 16 24 32 40 48 56 64 72 80 84)
DEFAULT_INDEX_LIST="alexol,lippol,btreeolc,artolc,libox,xindex"
ALL_INDEXES=("alex" "alexol" "lipp" "lippol" "btree" "btreeolc" "artunsync" "artolc" "libox" "xindex" "loft" "pgm" "rs" "hope" "lisa")

# Define workload configurations
declare -A WORKLOAD_PARAMS
WORKLOAD_PARAMS[readonly]="--read=1.0 --insert=0.0 --init_table_ratio=1 --operations_num=100000000"
WORKLOAD_PARAMS[read20]="--read=0.2 --insert=0.8 --init_table_ratio=0.5 --operations_num=200000000"
WORKLOAD_PARAMS[read40]="--read=0.4 --insert=0.6 --init_table_ratio=0.5 --operations_num=200000000"
WORKLOAD_PARAMS[balanced]="--read=0.5 --insert=0.5 --init_table_ratio=0.5 --operations_num=100000000"
WORKLOAD_PARAMS[read60]="--read=0.6 --insert=0.4 --init_table_ratio=0.5 --operations_num=200000000"
WORKLOAD_PARAMS[read80]="--read=0.8 --insert=0.2 --init_table_ratio=0.5 --operations_num=200000000"
WORKLOAD_PARAMS[writeonly]="--read=0.0 --insert=1.0 --init_table_ratio=0.5 --operations_num=100000000"
# workloads for memory consumption test (different bulkloading ratio)
WORKLOAD_PARAMS[init20]="--read=0.0 --insert=1.0 --init_table_ratio=0.2 --operations_num=0"
WORKLOAD_PARAMS[init40]="--read=0.0 --insert=1.0 --init_table_ratio=0.4 --operations_num=0"
WORKLOAD_PARAMS[init60]="--read=0.0 --insert=1.0 --init_table_ratio=0.6 --operations_num=0"
WORKLOAD_PARAMS[init80]="--read=0.0 --insert=1.0 --init_table_ratio=0.8 --operations_num=0"
WORKLOAD_PARAMS[init100]="--read=0.0 --insert=1.0 --init_table_ratio=1.0 --operations_num=0"
# workloads for range search
WORKLOAD_PARAMS[scanonly]="--read=0.0 --insert=0.0 --init_table_ratio=1 --operations_num=200000000 --scan=1 --scan_num=100"
WORKLOAD_PARAMS[ycsbe]="--read=0.0 --insert=0.05 --init_table_ratio=0.5 --operations_num=1000000 --scan=0.95 --scan_num=100"
WORKLOAD_PARAMS[scan10]="--read=0.0 --insert=0.0 --init_table_ratio=1 --operations_num=20000000 --scan=1 --scan_num=10"
WORKLOAD_PARAMS[scan100]="--read=0.0 --insert=0.0 --init_table_ratio=1 --operations_num=20000000 --scan=1 --scan_num=100"
WORKLOAD_PARAMS[scan500]="--read=0.0 --insert=0.0 --init_table_ratio=1 --operations_num=20000000 --scan=1 --scan_num=500"
WORKLOAD_PARAMS[scan1000]="--read=0.0 --insert=0.0 --init_table_ratio=1 --operations_num=20000000 --scan=1 --scan_num=1000"
WORKLOAD_PARAMS[scan1500]="--read=0.0 --insert=0.0 --init_table_ratio=1 --operations_num=20000000 --scan=1 --scan_num=1500"
WORKLOAD_PARAMS[scan2000]="--read=0.0 --insert=0.0 --init_table_ratio=1 --operations_num=20000000 --scan=1 --scan_num=2000"
WORKLOAD_PARAMS[scan2500]="--read=0.0 --insert=0.0 --init_table_ratio=1 --operations_num=20000000 --scan=1 --scan_num=2500"
WORKLOAD_PARAMS[scan10000]="--read=0.0 --insert=0.0 --init_table_ratio=1 --operations_num=20000000 --scan=1 --scan_num=10000"
WORKLOAD_PARAMS[scan100000]="--read=0.0 --insert=0.0 --init_table_ratio=1 --operations_num=20000000 --scan=1 --scan_num=100000"


# Parse command line arguments
SELECTED_TRACES=()
SELECTED_WORKLOADS=()
SELECTED_THREADS=()
SELECTED_INDEXES=""

function print_usage {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -t, --traces TRACES    Comma-separated list of traces to test (default: all)"
    echo "                         Available: osm, genome, fb, covid, longitudes-200M"
    echo "  -w, --workloads WLDS   Comma-separated list of workloads to test (default: all)"
    echo "                         Available: readonly, balanced, writeonly, read20, read40, read60, read80"
    echo "  -n, --threads THREADS  Comma-separated list of thread counts to test (default: all)"
    echo "  -i, --indexes INDEXES  Comma-separated list of indexes to test (default: all)"
    echo "                         Available: alexol, lippol, btreeolc, artolc, libox, loft, xindex"
    echo "  -o, --output FILENAME  Custom output filename (default: out_libox_MMDD.csv)"
    echo "  -h, --help             Show this help message"
    echo
    echo "Examples:"
    echo "  $0                          # Run all tests with defaults"
    echo "  $0 -t osm,fb                # Test only osm and fb traces"
    echo "  $0 -w readonly,balanced     # Test only readonly and balanced workloads"
    echo "  $0 -n 16,32,64              # Test only with 16, 32, and 64 threads"
    echo "  $0 -i libox,loft            # Test only libox and loft indexes"
    echo "  $0 -t osm -w writeonly -i libox  # Test osm trace with writeonly workload on libox only"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -t|--traces)
            IFS=',' read -r -a SELECTED_TRACES <<< "$2"
            shift 2
            ;;
        -w|--workloads)
            IFS=',' read -r -a SELECTED_WORKLOADS <<< "$2"
            shift 2
            ;;
        -n|--threads)
            IFS=',' read -r -a SELECTED_THREADS <<< "$2"
            shift 2
            ;;
        -i|--indexes)
            SELECTED_INDEXES="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_FILE="$RESULT_DIR/$2"
            shift 2
            ;;
        -h|--help)
            print_usage
            ;;
        *)
            echo "Unknown option: $1"
            print_usage
            ;;
    esac
done

# Use defaults if nothing is specified
if [ ${#SELECTED_TRACES[@]} -eq 0 ]; then
    # Remove .csv extension from trace names
    for trace in "${ALL_TRACES[@]}"; do
        # SELECTED_TRACES+=("${trace%}") # for umass
        SELECTED_TRACES+=("${trace%.csv}")
    done
else
    # Validate selected traces
    for trace in "${SELECTED_TRACES[@]}"; do
        found=0
        for valid_trace in "${ALL_TRACES[@]}"; do
            if [[ "${valid_trace%.csv}" == "$trace" ]]; then
                found=1
                break
            fi
        done
        if [[ $found -eq 0 ]]; then
            echo "Error: Invalid trace '$trace'. Available traces: ${ALL_TRACES[*]%.csv}"
            exit 1
        fi
    done
fi

if [ ${#SELECTED_WORKLOADS[@]} -eq 0 ]; then
    SELECTED_WORKLOADS=("${ALL_WORKLOADS[@]}")
else
    # Validate selected workloads
    for workload in "${SELECTED_WORKLOADS[@]}"; do
        if [[ ! " ${ALL_WORKLOADS[*]} " =~ " ${workload} " ]]; then
            echo "Error: Invalid workload '$workload'. Available workloads: ${ALL_WORKLOADS[*]}"
            exit 1
        fi
    done
fi

if [ ${#SELECTED_THREADS[@]} -eq 0 ]; then
    SELECTED_THREADS=("${ALL_THREAD_COUNTS[@]}")
else
    # Validate thread counts
    for thread in "${SELECTED_THREADS[@]}"; do
        if ! [[ "$thread" =~ ^[0-9]+$ ]]; then
            echo "Error: Invalid thread count '$thread'. Must be a positive integer."
            exit 1
        fi
    done
fi

# Process index list
INDEX_LIST="$DEFAULT_INDEX_LIST"
if [[ -n "$SELECTED_INDEXES" ]]; then
    # Validate selected indexes
    IFS=',' read -r -a INDEX_ARRAY <<< "$SELECTED_INDEXES"
    for idx in "${INDEX_ARRAY[@]}"; do
        if [[ ! " ${ALL_INDEXES[*]} " =~ " ${idx} " ]]; then
            echo "Error: Invalid index '$idx'. Available indexes: ${ALL_INDEXES[*]}"
            exit 1
        fi
    done
    INDEX_LIST="$SELECTED_INDEXES"
fi

# Create necessary directories
mkdir -p "$RESULT_DIR" "$RESULT_DIR/logs"

# Show test configuration
echo "==== Test Configuration ===="
echo "Traces:   ${SELECTED_TRACES[*]}"
echo "Workloads: ${SELECTED_WORKLOADS[*]}"
echo "Threads:  ${SELECTED_THREADS[*]}"
echo "Indexes:  $INDEX_LIST"
echo "Output file: $OUTPUT_FILE"
echo "========================="

# Run the tests
for trace_name in "${SELECTED_TRACES[@]}"; do
    # Add .csv extension for file lookups
    trace_file="$trace_name.csv"
    # trace_file="$trace_name" # for umass
    segment_file="segments_$trace_name.csv"
    
    # Find the trace and segment files
    TRACE_FILE="$TRACE_DIR/$trace_file"
    SEGMENT_FILE="$SEGMENT_DIR/$segment_file"
    
    # Check if files exist
    if [[ ! -f "$TRACE_FILE" ]]; then
        echo "Error: Trace file not found: $TRACE_FILE"
        continue
    fi
    
    # if [[ ! -f "$SEGMENT_FILE" ]]; then
    #     echo "Error: Segment file not found: $SEGMENT_FILE"
    #     continue
    # fi
    
    # Set binary path based on trace name
    if [[ "$trace_name" == "longitudes-200M" ]]; then
    # if [[ "$trace_name" == *"longitudes"* ]]; then
        PROGRAM="$PROGRAM_BASE/release/microbench_longitudes"
    else
        PROGRAM="$PROGRAM_BASE/release/microbench"
    fi
    
    # Check if program exists
    if [[ ! -f "$PROGRAM" ]]; then
        echo "Error: Program not found: $PROGRAM"
        echo "Did you build the project correctly?"
        continue
    fi
    
    echo "Starting tests for trace: $trace_name"
    echo "Trace file: $TRACE_FILE"
    echo "Segment file: $SEGMENT_FILE"
    echo "Program: $PROGRAM"
    echo "--------------------------------------------------------------------------------"
    
    LOG_FILE="$RESULT_DIR/logs/log_${trace_name}_$DATE_TAG.log"
    
    for workload_name in "${SELECTED_WORKLOADS[@]}"; do
        echo "Workload: $workload_name (${WORKLOAD_PARAMS[$workload_name]})"
        
        for thread_num in "${SELECTED_THREADS[@]}"; do
            echo "Running [$workload_name] with $thread_num threads on [$trace_name] using indexes [$INDEX_LIST]"
            
            {
                echo "===== Test: $trace_name | Workload: $workload_name | Threads: $thread_num | Indexes: $INDEX_LIST ====="
                echo "Start time: $(date '+%Y-%m-%d %H:%M:%S')"
                START_TIME=$(date +%s)
                
                $PROGRAM \
                --keys_file="$TRACE_FILE" \
                --keys_file_type=text \
                --config_file="$SEGMENT_FILE" \
                ${WORKLOAD_PARAMS[$workload_name]} \
                --output_path="$OUTPUT_FILE" \
                --table_size=-1 \
                --node_capacity=100000 \
                --top_k=0.05 \
                --density_factor=4 \
                --temp_node_cap=20 \
                --max_error_rs=32 \
                --min_line_len=10 \
                --thread_num="$thread_num" \
                --index="$INDEX_LIST"
                
                EXIT_CODE=$?
                END_TIME=$(date +%s)
                DURATION=$((END_TIME - START_TIME))
                
                echo "End time: $(date '+%Y-%m-%d %H:%M:%S')"
                echo "Duration: $DURATION seconds"
                
                if [[ $EXIT_CODE -ne 0 ]]; then
                    echo "WARNING: Test exited with code $EXIT_CODE"
                fi
                
                # Sleep to allow system to stabilize
                echo "Sleeping for 5 seconds before the next test..."
                sleep 5
                
                echo "==================================================================================="
            } 2>&1 | tee -a "$LOG_FILE"
        done
    done
done

echo "All tests completed. Results saved to: $OUTPUT_FILE"
echo "Logs directory: $RESULT_DIR/logs/"