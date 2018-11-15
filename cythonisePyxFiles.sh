#!/bin/bash
# cythonise the extractDiags.pyx file into extractDiags*.so

set -eu

cythonize -i extractDiags.pyx
