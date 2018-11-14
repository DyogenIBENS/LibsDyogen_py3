#!/bin/bash

set -euo pipefail

# First we need to know the directory of this file, which is the package directory.
thisfiledir="$( cd "$(dirname "$0")" ; pwd -P )"
echo $thisfiledir

PYTHONPATH="$(dirname $thisfiledir)":${PYTHONPATH:-} python3 -m LibsDyogen._myTests.all

echo "Done. Return code: $?" >&2
