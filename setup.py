import setuptools
from Cython.Build import cythonize


with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name="LibsDyogen_py3",
    version="0.0.1",
    author="Alexandra LOUIS, Matthieu MUFFATO, Joseph LUCAS, Guillaume LOUVEL, Hugues ROEST CROLLIUS",
    author_email="dyogen_git@biologie.ens.fr",
    description="Python 3 version of the Dyogen team library for comparative genomics. http://www.ibens.ens.fr/?rubrique43&lang=en",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/DyogenIBENS/LibsDyogen_py3",
    packages=setuptools.find_packages(),
    python_requires='>=3',
    setup_requires=['Cython'],
    install_requires=['enum34;python_version<"3.4"'],
    extras_require={'plotstats': ['numpy', 'scipy', 'matplotlib']},  # Only in bin/statsOnGenesInGenomes.py,
    #optional dependencies (handled in try-except)
                      #'java.lang',
                      #'psyco',  # development stopped in 2011
    ext_modules=cythonize('LibsDyogen/extractDiags.pyx'),
    scripts = ['scripts/ancGenesFromGeneTrees.py',
               'scripts/indexGenesUsingExtremityThatStartsTranscription.py',
               'scripts/newickSpeciesTree2phylTreeSpeciesTree.py',
               'scripts/nhxGeneTrees2phylTreeGeneTrees.py',
               'scripts/statsOnGenesInGenomes.py'
               ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        #"Operating System :: Microsoft :: Windows",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    include_package_data=True,
    zip_safe=False
)
