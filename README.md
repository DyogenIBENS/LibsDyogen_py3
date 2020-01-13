     _     _ _         ____
    | |   (_) |__  ___|  _ \ _   _  ___   __ _  ___ _ __
    | |   | | '_ \/ __| | | | | | |/ _ \ / _` |/ _ \ '_ \
    | |___| | |_) \__ \ |_| | |_| | (_) | (_| |  __/ | | |
    |_____|_|_.__/|___/____/ \__, |\___/ \__, |\___|_| |_|
                             |___/       |___/

**LibsDyogen_py3** version 0

It is a port for Python 3 of the original [LibsDyogen](https://github.com/DyogenIBENS/LibsDyogen) (python2.7) v.1 (27/01/2017).

Python v3.2 at least is needed.

This code may be freely distributed and modified under the terms of the GNU
General Public License version 3 (GPL v3) and the CeCILL licence version 2 of
the CNRS. These licences are contained in the files:

1. LICENSE-GPL.txt (or <http://www.gnu.org/licenses/gpl-3.0-standalone.html>)
2. LICENCE-CeCILL.txt (or <http://www.cecill.info/licences/Licence_CeCILL_V2-en.html>)

Copyright for this code is held jointly by the Dyogen (DYnamic and Organisation
of GENomes) team of the Institut de Biologie de l'École Normale Supérieure
(IBENS) 46 rue d'Ulm Paris, and the individual authors.

Copyright © 2015 IBENS/Dyogen and authors

Authors' contributions:
-----------------------

* Matthieu MUFFATO: started the library and created: myGenomes, myTools,
  myPsOutput, myProteinTree, myPhylTree, myKaryoDrawer, myFile, myGenomes,
  myGraph, myMultiprocess, myMaths
* Joseph LUCAS: Updated the library and organised it, created: myDiags,
  myLightGenomes, myCondor, myGeneTeams, myGenomesDrawer, myIntervals,
  myMapping, myProbas, mySvgDrawer and updated: myTools
* Nga THI THUY NGUYEN: Updated some functions of myDiags and myGraph
* Lucas TITTMANN: Updated some functions of myDiags
* Guillaume LOUVEL: conversion to Python 3 and various updates in myFile,
    myTools, myPhylTree, myProteinTree.
* Hugues ROEST CROLLIUS: supervisor


Mail: dyogen_git[at]biologie.ens.fr

About Licences:
---------------

## Justify the use of snippets from stack overflow

<http://programmers.stackexchange.com/questions/12171/how-does-fair-use-apply-to-code-snippets>

Installation:
-------------

## For Ubuntu and Debian-like Linux distributions

### Requirements

```bash
# Install python3
sudo apt-get install python3

# Install cython (not mandatory)
# http://docs.cython.org/src/quickstart/install.html
sudo apt-get install python3-dev cython3

# Install git
sudo apt-get install git

# Clone this repository
git clone https://github.com/DyogenIBENS/LibsDyogen_py3.git

# Install
cd LibsDyogen_py3/
python3 setup.py install --user
# Or install in developer mode (the current folder is just symlinked)
python3 setup.py develop --user
```

#### Install optional plugable softwares (related to PhylDiag):

Homology teams:

```bash
HOMOLOGYTEAMS_DIR="$HOME/Libs"
cd ${HOMOLOGYTEAMS_DIR}
wget http://euler.slu.edu/~goldwasser/homologyteams/homologyteams-1.1.zip
unzip homologyteams-1.1.zip
cd homologyteams-1.1/src
make

# Make sure the executable is in your path, e.g:
echo "export PATH=\"\$PATH:${HOMOLOGYTEAMS_DIR}/homologyteams-1.1/src\"" >> ~/.bashrc
```

