#!/usr/bin/env python3


"""Test all imports assuming the PATH has an entry **containing** the `LibsDyogen` folder.
"""


import sys

if sys.version_info[0] < 3:
    raise RuntimeError("This version is designed for Python version >= 3")

# Remove possibility to import from files located in the current directory
# (implicit relative import)
sys.path.remove('')

from .. import \
               myCondor, \
               myCyntenator, \
               myFile, \
               myGenomes, \
               myGraph, \
               myIntervals, \
               myKaryoDrawer, \
               myMapping, \
               myMaths, \
               myMultiprocess, \
               myPhylTree, \
               myProbas, \
               myProteinTree, \
               myPsOutput, \
               mySvgDrawer, \
               myTools, \
               walktrap, \
               myADHoRe, \
               myDiags, \
               myLightGenomes, \
               myGeneTeams, \
               myGenomesDrawer
