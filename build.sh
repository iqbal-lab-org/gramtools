#!/usr/bin/env bash
set -eu -o pipefail

usage(){
    echo "Builds the gramtools backend"
    echo "usage: $0 [build_target] [build_mode]"
    echo "\tbuild_target: $TARGETS (default: gram)"
    echo "\tbuild_mode: $BUILD_TYPES (default: REL_WITH_ASSERTS)"
    exit 0
}

TARGETS="gram|all|test_main"
BUILD_TYPES="DEBUG|REL_WITH_ASSERTS"
TARGET=${1:-gram}
BUILD_TYPE=${2:-REL_WITH_ASSERTS}
[[ "$BUILD_TYPE" =~ $BUILD_TYPES ]] || usage
[[ "$TARGET" =~ $TARGETS ]] || usage

BASE_DIR=$(realpath $(dirname $0))
CUR_DIR=$(pwd)

EXIT_CODE=0
conan >/dev/null 2>&1 || EXIT_CODE=$?
if [[ "$EXIT_CODE" != 0 ]]; then 
    echo "Missing conan: first run 'pip install conan'"
    exit 1
fi


BUILD_DIR="${BASE_DIR}/build"
STDOUT_FILE="${BUILD_DIR}/build_stdout.txt"
echo "Building in: ${BUILD_DIR}" >&1
echo "Writing stdout to: ${STDOUT_FILE}" >&1
mkdir -p "${BUILD_DIR}" && cd "${BUILD_DIR}"

conan install .. -s compiler.libcxx=libstdc++11 --build=missing > "$STDOUT_FILE"
CC=gcc CXX=g++ cmake -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" .. >> "$STDOUT_FILE"
make -j 4 "${TARGET}" | tee -a "$STDOUT_FILE"

cd "$CUR_DIR"
