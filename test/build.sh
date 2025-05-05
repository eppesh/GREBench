#!/bin/bash

# Determine build type from argument
BUILD_TYPE="Release"
BUILD_DIR="../build/release"

if [[ "$1" == "debug" ]]; then
    BUILD_TYPE="Debug"
    BUILD_DIR="../build/debug"
fi

# Create build directory if it doesn't exist
mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

# Run cmake and make
echo "Building in ${BUILD_TYPE} mode..."
cmake -DCMAKE_BUILD_TYPE=${BUILD_TYPE} ../..
make -j$(nproc)
