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
* Hugues ROEST CROLLIUS: supervisor

Mail: hrc[at]ens.fr or jlucas[at]ens.fr

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

# Install git
sudo apt-get install git

# Install cython (not mandatory)
# http://docs.cython.org/src/quickstart/install.html
sudo apt-get install cython3
```

### Setup

#### Recommended installation path

```bash
# This path is automatically used by Python.
# Change the Python version number according to your version.
LIBSDIR=$HOME/.local/lib/python3.5/site-packages/
```

#### Custom installation path

```bash
# Choose a path for the installation (here it is ${HOME}/Libs)
LIBSDIR="${HOME}/Libs"
```

Add LIBSDIR to the `PYTHONPATH` environment variable, e.g in the `~/.bashrc` file:

```bash
export PYTHONPATH="${LIBSDIR}:${PYTHONPATH}"
```

#### Clone the LibsDyogen library

```bash
# Create the root folder of LibsDyogen
mkdir -p ${LIBSDIR}
cd ${LIBSDIR}
git clone <url to this repository> LibsDyogen
#(example: https://gitlab.com/DyogenIBENS/LibsDyogen_py3.git)
```

#### Compile the extractDiags.pyx file

```bash
cd LibsDyogen
cythonize -i extractDiags.pyx
```

#### Install optional plugable softwares:

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

