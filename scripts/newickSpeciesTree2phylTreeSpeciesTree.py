#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# from PhylDiag v1.02
# python 3
# Copyright © 2013 IBENS/Dyogen : Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

""" convert species tree (newick to phylTree) or (phylTree to newick) """

from LibsDyogen import myFile
from LibsDyogen import myTools
from LibsDyogen import myPhylTree

arguments = myTools.checkArgs([("phylTree.conf",file)], [("fromNewick",bool,True)], __doc__)

phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])


if arguments["fromNewick"]:

    # print into phylTree format, with tabulations
    def do(node, indent):
        node = node.replace(myPhylTree.SYMBOL6X, "")
        node = node.replace(myPhylTree.SYMBOL2X, "")
        names = myFile.myTSV.printLine([node] + [x for x in phylTree.commonNames.get(node,"") if isinstance(x, str) and (x != node)], delim="|")
        if node in phylTree.listSpecies :
            print(("\t" * indent) + str(names))
        elif node in phylTree.items:
            print(("\t" * indent) + str(names) + "\t" + str(int(phylTree.ages[node])))
            for (f,_) in phylTree.items[node]:
                do(f, indent+1)
    do(phylTree.root, 0)

else:
    # return the tree into the newick tree
    def convertToFlatFile(anc):
        a = phylTree.fileName[anc]
        if anc in phylTree.listSpecies:
            return a
        else:
            return "(" + ",".join([convertToFlatFile(e) + ":" + str(l) for (e,l) in phylTree.items[anc]]) + ")%s" % a #" |%d" % (a,phylTree.ages[anc])
    print(convertToFlatFile(phylTree.root), ";")
