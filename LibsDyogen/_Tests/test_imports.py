#!/usr/bin/env python3


"""Test all imports assuming the PATH has an entry **containing** the `LibsDyogen` folder.
"""


import sys

def test_python_is_version_3():
    assert sys.version_info[0] >= 3, "This version is designed for Python version >= 3"

# Remove possibility to import from files located in the current directory
# (implicit relative import)
#sys.path.remove('')


def test_relative_imports():
    from .. import \
                   myFile, \
                   myTools, \
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
                   myCondor, \
                   myLightGenomes, \
                   walktrap, \
                   myDiags, \
                   myADHoRe, \
                   myCyntenator, \
                   myGeneTeams, \
                   myGenomesDrawer


if __name__ == '__main__':
    test_python_is_version_3()
    test_relative_imports()

