#!/bin/bash
# cythonise the extractDiags.pyx file into extractDiags.c
# the first argument should be the path to the root of LibsDyogen

set -eu

PATH_LIBSDYOGEN=${1}
cython3 ${PATH_LIBSDYOGEN}/extractDiags.pyx
# C compilation of extractDiags.c into extractDiags.so that can be used as a python package
gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/usr/include/python2.7 -o ${PATH_LIBSDYOGEN}/extractDiags.so ${PATH_LIBSDYOGEN}/extractDiags.c
