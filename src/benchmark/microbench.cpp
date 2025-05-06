#include "./benchmark.h"

int main(int argc, char **argv) {
    Benchmark <int64_t, int64_t> bench;
    bench.parse_args(argc, argv);
    bench.run_benchmark();
}

