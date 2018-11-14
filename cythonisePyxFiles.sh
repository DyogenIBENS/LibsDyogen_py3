#!/bin/bash
# cythonise the extractDiags.pyx file into extractDiags.c
# the first argument should be the path to the root of LibsDyogen

set -eu

PATH_LIBSDYOGEN=${1:-.}
# Path to the most recent Python3 development headers (apt install python3-dev)
PATH_PYTHON3HEADER=$(find /usr/include/ -maxdepth 1 -name "python3*" | sort -rV | head -1)

echo "* Working directory: '${PATH_LIBSDYOGEN}'
* Path to Python 3 development files: '${PATH_PYTHON3HEADER}'" >&2

cython3 ${PATH_LIBSDYOGEN}/extractDiags.pyx
# C compilation of extractDiags.c into extractDiags.so that can be used as a python package
gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I${PATH_PYTHON3HEADER} -o ${PATH_LIBSDYOGEN}/extractDiags.so ${PATH_LIBSDYOGEN}/extractDiags.c
