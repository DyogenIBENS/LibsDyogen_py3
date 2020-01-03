#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# from PhylDiag v1.02
# python 3
# Copyright © 2015 IBENS/Dyogen : Joseph LUCAS, and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France
import collections

"""
        Rank genes in chromosomes by increasing coordinate of the 5' extremity
        (the TSS) and assert that the beg (2nd column) and end (3rd column) of
        genes are consistent with their orientation.
"""

import sys, itertools
from LibsDyogen import myTools, myFile, myLightGenomes

def intOrNone(string):
    return int(string) if string != 'None' else None

def readerDependingOnFileWithDebAndEnd(fileName):
        flb = myFile.firstLineBuffer(myFile.openFile(fileName, 'r'))
        c = flb.firstLine.split("\t")
        if len(c) == 6:
            print("(c, beg, end, s, gName, transcriptName) -> (c, s, gName)", file=sys.stderr)
            # c, beg, end, s,  gName, transcriptName
            reader = myFile.myTSV.readTabular(fileName, [str, str, str, str, str, str])
            reader = ((c, intOrNone(beg), intOrNone(end), intOrNone(strand), gName) for (c, beg, end, strand, gName, tName) in reader)
        elif len(c) == 5:
            print("(c, beg, end, s, gName) -> (c, s, gName)", file=sys.stderr)
            # c, beg, end, s,  gName
            tmpReader = myFile.myTSV.readTabular(fileName, [str, str, str, str, str])
            # check, with the first line, if there are several gene names (the format genome of Matthieu contains several gene names)
            (c, beg, end, strand, gNames) = next(tmpReader)
            severalNames = True if len(gNames.split(' ')) > 0 else False
            reader = itertools.chain([(c, beg, end, strand, gNames)], tmpReader)
            if severalNames:
                # if gNames contains more than one gene name, only take the first gene name
                reader = ((c, intOrNone(beg), intOrNone(end), intOrNone(strand), gNames.split(' ')[0]) for (c, beg, end, strand, gNames) in reader)
            else:
                reader = ((c, intOrNone(beg), intOrNone(end), intOrNone(strand), gName) for (c, beg, end, strand, gName) in reader)
        else:
            raise ValueError("%s file is badly formatted" % fileName)
        return reader

# arguments
arguments = myTools.checkArgs(
    [
        ("genome", myTools.file),
    ],
    [
        ("orderGenesByIncreasingTranscriptionStart", bool, True),
        ("orderChromosomesBy", bool, 'decreasingNbOfGenes'),
        ("out:genome", str, 'genes.outPutGenome.list.bz2'),
        ("removeUnofficialChrNames", bool, False)
    ],
    __doc__)
assert arguments['orderChromosomesBy'] in {'decreasingNbOfGenes', 'names'}

# 1) Load the genome and ensure the consistency of gene strand ans gene coordinates
unofficialChromRemoved = False
iniGenomeLength = 0
genomeListByChr = myTools.DefaultOrderedDict(list)
reader = readerDependingOnFileWithDebAndEnd(arguments['genome'])
cptNbEditedGenes = 0
for (chr, x1, x2, s, gNames) in reader:
    if arguments['removeUnofficialChrNames'] and myLightGenomes.contigType(chr) != myLightGenomes.ContigType.Chromosome:
        unofficialChromRemoved = True
        continue
    assert isinstance(x1, int)
    assert isinstance(x2, int)
    if x1 == x2:
        raise ValueError('The current gene (%s) has a null length' % gNames)
    if s in {+1, None}:
        if x1 < x2:
            (beg, end) = (x1, x2)
        else:
            assert x2 < x1
            (beg, end) = (x2, x1)
            cptNbEditedGenes += 1
    else:
        assert s == -1
        if x2 < x1:
            (beg, end) = (x1, x2)
        else:
            assert x1 < x2
            (beg, end) = (x2, x1)
            cptNbEditedGenes += 1
    genomeListByChr[chr].append((beg, end, s, gNames))
    iniGenomeLength += 1
if cptNbEditedGenes > 0:
    print('%s genes have been edited because gene strand was not coherent with gene.beg and gene.end columns' % cptNbEditedGenes, file=sys.stderr)
if unofficialChromRemoved:
    print('At least one unofficial chromosome (int(chr) > 100 or discarded name in myLightGenome.contigType()) has been removed', file=sys.stderr)

rankOfChrHasChanged = False
if arguments['orderChromosomesBy'] == 'decreasingNbOfGenes':
    if not myTools.isSorted(list(genomeListByChr.items()), key=lambda x: len(x[1]), stricly=False, increasingOrder=False):
        genomeListByChr = collections.OrderedDict(sorted(list(genomeListByChr.items()), key=lambda x: len(x[1]), reverse=True))
        rankOfChrHasChanged = True
elif arguments['orderChromosomesBy'] == 'names':
    if not myTools.isSorted(list(genomeListByChr.items()), key=lambda x: myTools.keyNaturalSort(x[0])):
        genomeListByChr = collections.OrderedDict(sorted(list(genomeListByChr.items()), key=lambda x: myTools.keyNaturalSort(x[0]), reverse=True))
        rankOfChrHasChanged = True
if rankOfChrHasChanged:
    print('The rank of at least one chromosome has changed while sorting chrNames using the length of chromosomes', file=sys.stderr)

# 2) If necessary rank genes by increasing beg coordinates
geneRankHasChanged = False
for chr, chrom in  genomeListByChr.items():
    if not myTools.isSorted(chrom, key=lambda x: x[0], stricly=False, increasingOrder=True):
        chrom.sort(key=lambda x: x[0])
        geneRankHasChanged = True
if geneRankHasChanged:
    print('The rank of at least one gene has changed while sorting using the 5\' extremities', file=sys.stderr)
assert sum(len(chrom) for chrom in list(genomeListByChr.values())) == iniGenomeLength

# 3) Print the genome
f = myFile.openFile(arguments['out:genome'], 'w')
for chr, chrom in genomeListByChr.items():
    for (beg, end, s, gNames) in chrom:
        print(myFile.myTSV.printLine([chr, beg, end, s, gNames]), file=f)
f.close()

