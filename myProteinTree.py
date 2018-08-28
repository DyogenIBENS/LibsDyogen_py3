# -*- coding: utf-8 -*-

# LibsDyogen version 1.0 (6/11/2015)
# python v2.7 at least is needed
# Copyright © 2015 IBENS/Dyogen : Matthieu MUFFATTO, Joseph LUCAS and Hugues ROEST CROLLIUS
# mail : jlucas@ens.fr
# Licences GLP v3 and CeCILL v2

from __future__ import print_function


"""Manage gene trees in the form of (data,info)"""


import sys
import collections

from . import myTools
from . import myFile


def hasLowScore_alwaysTrue(tree, node):
    """Default function `hasLowScore` for `ProteinTree.rebuildTree`.
    Always return True."""
    return True


def getIntermediateAnc(phylTree, previousAnc, lastWrittenAnc, newAnc, isDuplication):
    """return the ancestor names since the previousAnc"""

    if lastWrittenAnc: # cruise mode
        # genes of the species between lastWritten and newAnc are recorded
        # exepted if it is a duplication node (see after)
        toWrite =list(phylTree.dicLinks[lastWrittenAnc][newAnc][1:])

    # in the first species of the gene tree or at the first node outside the
    # first species
    elif not lastWrittenAnc:
        if previousAnc == None: # at the root
            # the gene is recorded except if it is a duplication (see after)
            toWrite =list([newAnc])

        elif previousAnc: # not at the root

            if previousAnc == newAnc: # still in the first species
                # the gene is recorded except if it is a duplication (see after)
                toWrite=list([newAnc])

            elif previousAnc != newAnc: # in a child species
                # the genes of the species between previousAnc and newAnc are
                # recorded and newAnc is removed if the node is a duplication
                # node (see after)
                toWrite=list(phylTree.dicLinks[previousAnc][newAnc])

    # if the current node is a duplication newAnc is not recorded in the list
    # to be written
    if isDuplication:
        toWrite.remove(newAnc)

    root=False # root refers to terminal genes of the first species
    # in the first species or at the first node outside the first species
    if not lastWrittenAnc:
        if not isDuplication: # if it is a speciation node
            root = True
        # if at the first node outside the first species and if this node is a
        # duplication node
        if previousAnc != newAnc and isDuplication:
            root = True

    if len(toWrite) == 0:
        newLastWritten = lastWrittenAnc
    else:
        newLastWritten = toWrite[-1]

    return (toWrite,newLastWritten,root)


class ProteinTree:
    """class for managing a forest of gene trees"""

    def __init__(self, data=None, info=None, root=None):
        self.data = {} if data is None else data
        self.info = {} if info is None else info
        self.root = root


    def doBackup(self):
        self.backRoot = self.root
        self.backData = dict((node,values[:]) for (node,values) in self.data.items())
        self.backInfo = dict((node,values.copy()) for (node,values) in self.info.items())


    def printTree(self, f, node=None):
        """print the tree into the phylTree format (with tabulations)"""
            
        def rec(n, node):

            indent = "\t" * n
            # id of the node
            print("%sid\t%d" % (indent, node), file=f)
            # informations
            print("%sinfo\t{%s}" % (indent, ", ".join(repr(key) + ": " + repr(value) for (key, value) in sorted(self.info[node].items()))), file=f)
            # children
            for (g,d) in self.data.get(node,[]):
                print("%s\tlen\t%f" % (indent, d), file=f)
                rec(n+1, g)

        rec(0, self.root if node is None else node)
        try:
            f.flush()
        except AttributeError:
            pass

    def printNewick(self, f, root=None, withDist=True, withTags=False, withAncSpeciesNames=False, withAncGenesNames=False, withID=False):
        """print the tree into the Newick format (with parentheses)"""
        NHX = {"Duplication": "D",
               "Bootstrap": "B",
               "taxon_name": "S",
               "duplication_confidence_score": "SIS",
               "dubious_duplication": "DD"}
        def rec(node):
            if node in self.data:
                return "(" + ",".join(
                        rec(g)
                        + ((":%f" % l) if withDist else "")
                        + ("[&&NHX:"
                            + (":I=%d" % node if withID else "")
                            + ":".join(
                                ("%s=%s" % (
                                    (NHX[tag],self.info[g][tag]) 
                                     if tag!="Duplication" 
                                     else (NHX[tag],"N" if self.info[g][tag]== 0 else "Y")
                                           )
                                ).replace(" ", ".") for tag in NHX if tag in self.info[g]
                                      ) + "]" if withTags else ""
                          )
                          for (g,l) in self.data[node]
                        ) + ")" + (
                        self.info[node]["taxon_name"].replace(' ', '.') if withAncSpeciesNames and ("taxon_name" in self.info[node]) else ''
                        )+(
                         self.info[node].get('family_name', '').split("/")[0] if withAncGenesNames and ("taxon_name" in self.info[node]) else ''
                        )
            else:
                return self.info[node].get('gene_name', 
                           (self.info[node]['taxon_name'].replace(' ', '.') if withAncGenesNames and ("taxon_name" in self.info[node]) else '') +
                           (self.info[node].get('family_name', '').split("/")[0] if withAncGenesNames and ("taxon_name" in self.info[node]) else '')
                           )

        if root is None:
            root = self.root
        print(rec(root) + ("[&&NHX:" + ":".join(("%s=%s" % ((NHX[tag],self.info[root][tag]) if tag!="Duplication" else (NHX[tag],"N" if self.info[root][tag]== 0 else "Y"))).replace(" ", ".") for tag in NHX if tag in self.info[root]) + "]" if withTags else "") + ";", file=f)
        try:
            f.flush()
        except AttributeError:
            pass


    #FIXME
    def printNewickTree(self, f, node=None):
        """print tree into the Newick format"""
        genes = []
        def rec(node):
            if node not in self.data:
                genes.append(self.info[node]['gene_name'])
                return self.info[node]['gene_name']
            else:
                return "(" + ",".join([rec(x) + ":" + str(l) for (x,l)  in self.data[node]]) + ") " + self.info[node]['family_name']
        tr = rec(self.root if node is None else node)
        print(" ".join(genes), file=f)
        print(tr, ";", file=f)


    def compactTree(self, phylTree, node=None):
        """Compact a tree by removing intermediary nodes that have only one
        child"""

        def do(node):
            # end of the process on one leaf
            if node not in self.data:
                return False

            flag = False
            # recursive calls
            for (gg,_) in self.data[node]:
                flag |= do(gg)

            if len(self.data[node]) > 1:
                return flag

            # edition of the current node
            (g,l) = self.data[node][0]
            if g in self.data:
                self.data[node] = [(gg,ll+l) for (gg,ll) in self.data[g]]
                del self.data[g]
            else:
                del self.data[node]
            self.info[node] = self.info[g]
            del self.info[g]
            return True

        return do(self.root if node is None else node)


    def renameTree(self, phylTree, node=None):
        """rename all nodes in order to make them match with the common
        ancestors of their descendants"""

        def do(node):
            # end of the process on one leaf
            if node not in self.data:
                return False

            flag = False
            # recursive calls
            for (g,_) in self.data[node]:
                flag |= do(g)

            # rename the current node
            newName = phylTree.lastCommonAncestor( [self.info[g]['taxon_name'] for (g,_) in self.data[node]] )
            flag |= (self.info[node]['taxon_name'] != newName)
            self.info[node]['taxon_name'] = newName

            return flag

        return do(self.root if node is None else node)


    def flattenTree(self, phylTree, rec, node=None):
        """Flatten a node and direct children if they represent the same taxon
        and if their is no duplication For instance:
            ((Eutheria1,Eutheria2)XA,(Eutheria3,Eutheria4)XB)XC is transformed
            into (Eutheria1,Eutheria2,Eutheria3,Eutheria4) only if XA, XB et XC
            are speciation nodes"""

        def do(node):
            # end of the process on one leaf
            if node not in self.data:
                return False
            assert len(self.data[node]) > 0

            flag = False
            # recursive calls
            if rec:
                for (g,_) in self.data[node]:
                    flag |= do(g)

            self.info[node]['taxon_name'] = phylTree.lastCommonAncestor( [self.info[g]['taxon_name'] for (g,_) in self.data[node]] )

            # if it is a true duplication, their is nothing more to do
            if self.info[node]['Duplication'] >= 2:
                return flag

            newData = []
            taxonName = self.info[node]['taxon_name']
            for (g,d) in self.data[node]:
                inf = self.info[g]
                # 2x the same taxon and no duplication

                if (inf['taxon_name'] == taxonName) and (inf['Duplication'] < 2) and (g in self.data):

                    newData.extend([(g2,d+d2) for (g2,d2) in self.data[g]])
                    del self.data[g]

                    self.info[node].update(self.info[g])
                    del self.info[g]
                    flag = True
                else:
                    newData.append( (g,d) )

            assert len(newData) == len(set(g for (g,_) in newData)), newData
            assert len(newData) > 0, (node,self.data[node],newData)

            self.info[node]['Duplication'] = 0
            self.data[node] = newData

            if len(self.data[node]) > 0:
                assert self.info[node]['taxon_name'] == phylTree.lastCommonAncestor( [self.info[g]['taxon_name'] for (g,_) in self.data[node]] )

            return flag

        return do(self.root if node is None else node)


    def rebuildTree(self, phylTree, hasLowScore=hasLowScore_alwaysTrue,
                    node=None):
        """give back the expected topology to the tree (to match the species tree)
          gather equivalent nodes under the same child.
          
          hasLowScore: function that takes (tree, node) and return True/False.
                       example: check that the duplication score of a node is
                       above a given threshold.
          """

        def do(node):

            # end of the process on one leaf
            if node not in self.data:
                return False

            flag = False

            # only if it is not a true duplication the children are changed
            if self.info[node]['Duplication'] < 2:

                # the node will be a speciation except special case below
                self.info[node]['Duplication'] = 0

                # the children are redefined for *our* phylogenetic tree by
                # organising the children into packs
                children = collections.defaultdict(list)  # *gene* children, classified by taxon lineage
                anc = self.info[node]['taxon_name']

                lchildren = phylTree.items.get(anc, [])   # taxon children
                for (g,d) in self.data[node]:
                    gname = self.info[g]['taxon_name']
                    # Place back the gene child if it has a taxon that fits.
                    for (a,_) in lchildren:
                        # Check that the taxon of child genes are
                        # **descendants** of the expected child taxa.
                        if phylTree.isChildOf(gname, a):
                            children[a].append( (g,d) )
                            break
                    # There was no 'break' in the above loop: no gene child fits the expected taxa:
                    else:
                        # Special case:
                        # we can be here only if g is in the same ancestor than
                        # node and if g is a duplication,
                        # this usually entails Duplication >= 2, except for the new created nodes
                        assert (gname == anc), "ERROR: name!=anc [%s / %s / %s]" % (node, anc, gname)
                        # false: g is not always a duplication !
                        # assert self.info[g]['Duplication'] >= 2
                        # the current node will thus be a duplication node
                        self.info[node]['Duplication'] = 3
                        children[anc].append( (g,d) )

                #print >> sys.stderr, "NEW CHILD", child
                # len(child):
                #  1 -> only anc
                #  2 or 3 among C1/C2/anc
                assert (len(children) != 1) or (anc in children), "ERROR: 1=anc [%s / %s]" % (node, children)
                assert (len(children) <= (1+len(lchildren))), "ERROR: len>(1+nbChildren) [%s / %s]" % (node, children)

                todo = []  # flattened nodes

                if len(children) > 1:
                    if anc in children:
                        lst1 = children.pop(anc)
                        lst2 = []
                        for tmp in children.values():
                            lst2.extend(tmp)
                        items = [(anc,lst1), (anc,lst2)]
                    else:
                        items = list(children.items())

                    newData = set()
                    for (anc,lst) in items:
                        if len(lst) == 1:
                            newData.add( lst[0] )
                        elif len(lst) > 1:
                            for (g,l) in self.data[node]:
                                if (g in self.data) and (self.data[g] == lst):
                                    newData.add( (g,l) )
                                    break
                                if g in self.data:
                                    assert sorted(self.data[g]) != sorted(lst)
                            else:
                                # Need to insert the correct common ancestor
                                global nextNodeID
                                nextNodeID += 1
                                # Put the common anc in the middle of the shortest branch
                                length = min([d for (_,d) in lst]) / 2
                                self.data[nextNodeID] = [(g,d-length) for (g,d) in lst]
                                anc = phylTree.lastCommonAncestor([self.info[g]['taxon_name'] for (g,_) in lst])
                                self.info[nextNodeID] = {'taxon_name':anc}
                                self.info[nextNodeID]["Duplication"] = 1 if hasLowScore(self, nextNodeID) else 3
                                todo.append(nextNodeID)
                                newData.add( (nextNodeID,length) )
                                self.flattenTree(phylTree, False,  nextNodeID)
                                flag = True
                    assert len(newData) == len(set(g for (g,_) in newData)), newData
                    # assign newData (reordered)
                    self.data[node] = [x for x in self.data[node] if x in newData] + list(newData.difference(self.data[node]))
                    for x in todo:
                        if hasLowScore(self, x):
                            self.info[x]["Duplication"] = 0
                            self.flattenTree(phylTree, False, x)
            # recursive calls
            for (g,_) in self.data[node]:
                flag |= do(g)
            return flag
        return do(self.root if node is None else node)

nextNodeID = -1

@myTools.deprecated
def printTree(ft, data, info, root):
    ProteinTree(data, info, root).printTree(ft)



def getDupSuffix(n, upper=False):
    """return the suffix associated to a duplication rank n.
    ('a' -> 'z', 'aa' -> 'az', 'ba' ...)"""
    base = 64 if upper else 96
    assert 1 <= n
    s = "."
    while n > 26:
        s = s + chr(base + (n - 1) % 26)
        n = 1 + (n - 1)/26
    return s + chr(base + n)


def loadTree(name):
    """Load a LibsDyogen ProteinTree from a file"""

    ns = myTools.Namespace()

    # read the next line of the file (and bufferise the next one)
    def nextLine():
        # Uses variables `f` and `ns` from the outside scope
        old = ns.curr
        try:
            l = ""
            while (l == "") or l.startswith("#"):
                # the final '\n' is removed and we cut owing to the '\t'
                l = next(f).replace('\n', '')
            l = l.split('\t')
            # the triplet (indentation,key,value) is recorded
            try:
                ns.curr = (len(l)-2, l[-2], l[-1])
            except IndexError:
                print('Badly formatted line:\n%r' % '\t'.join(l), file=sys.stderr)
                raise
        except StopIteration:
            ns.curr = None
        return old

    # the analysing process of the lines of the file
    def recLoad(tree, indent):
        # Uses variable `ns` from the outside scope

        # id of the point
        currID = int(nextLine()[2])
        # associated informations
        tree.info[currID] = eval(nextLine()[2])

        # children ?
        child = []
        while (ns.curr != None) and (ns.curr[0] == indent+1):
            length = float(nextLine()[2])
            child.append((recLoad(tree, indent+1), length))
        if len(child) > 0:
            tree.data[currID] = child

        return currID


    print("Loading the forest of gene trees %s ..." % name, end=' ', file=sys.stderr)
    f = myFile.openFile(name, "r") if isinstance(name, str) else name
    ns.curr = None
    nextLine()
    n = (0,0,0)
    while True:
        tree = ProteinTree()
        tree.root = recLoad(tree, 0)
        yield tree
        n = (n[0]+1, n[1]+len(tree.data), n[2]+len(tree.info)-len(tree.data))
        if ns.curr == None:
            break
    print("%d roots, %d internal nodes, %d leaves: OK" % n, file=sys.stderr)

    f.close()
