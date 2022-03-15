#!/usr/bin/env bash
set -eu -o pipefail

TARGETS="gram|all|test_main"
BUILD_TYPES="DEBUG|REL_WITH_ASSERTS"

usage(){
    echo "Builds the gramtools backend"
    echo "usage: $0 [--target ($TARGETS)] [--build_type ($BUILD_TYPES)] [--static]"
    echo -e "\tDefaults: --target gram --build_type REL_WITH_ASSERTS"
    echo -e "\tPass --static to make a standalone executable"
    exit 0
}

while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        --target)
        TARGET="$2"
        [[ "$TARGET" =~ $TARGETS ]] || usage
        shift
        shift
        ;;
        --build_type)
        BUILD_TYPE="$2"
        [[ "$BUILD_TYPE" =~ $BUILD_TYPES ]] || usage
        shift
        shift
        ;;
        --static)
        STATIC="-static"
        shift
        ;;
        *)    # unknown option
        usage
        ;;
    esac
done

TARGET=${TARGET:-gram}
BUILD_TYPE=${BUILD_TYPE:-REL_WITH_ASSERTS}
STATIC=${STATIC:- }
[[ "$STATIC" =~ -static| ]] || usage

BASE_DIR=$(realpath $(dirname $0))
CUR_DIR=$(pwd)

EXIT_CODE=0
conan >/dev/null 2>&1 || EXIT_CODE=$?
if [[ "$EXIT_CODE" != 0 ]]; then 
    echo "Missing conan: first run 'pip install conan'"
    exit 1
fi


BUILD_DIR="${BASE_DIR}/cmake-build"
STDOUT_FILE="${BUILD_DIR}/build_stdout.txt"
echo "Building in: ${BUILD_DIR}" >&1
echo "Writing stdout to: ${STDOUT_FILE}" >&1
mkdir -p "${BUILD_DIR}" && cd "${BUILD_DIR}"

CMAKE_OPTS="-DCMAKE_EXE_LINKER_FLAGS=${STATIC} -DCMAKE_BUILD_TYPE=${BUILD_TYPE}"

conan install .. -s compiler.libcxx=libstdc++11 --build=missing > "$STDOUT_FILE"
CC=gcc CXX=g++ cmake "$CMAKE_OPTS" .. | tee -a "$STDOUT_FILE"
make -j 4 "${TARGET}" | tee -a "$STDOUT_FILE"

[[ "${TARGET}" == "test_main" ]] && ctest -V

cd "$CUR_DIR"
