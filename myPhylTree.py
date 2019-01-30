# -*- coding: utf-8 -*-

# LibsDyogen_py3 version 0 (2018/11/16)
# Copyright © 2015 IBENS/Dyogen : Matthieu MUFFATTO, Joseph LUCAS and Hugues ROEST CROLLIUS
# mail : jlucas@ens.fr
# Licenses GPL v3 and CeCILL v2

from __future__ import print_function


import os
import sys
import itertools
import collections

from . import myFile

SYMBOL6X = '.'
SYMBOL2X = '*'

GeneSpeciesPosition = collections.namedtuple("GeneSpeciesPosition", ['species', 'chromosome', 'index'])

class PhylogeneticTree:
    """class mother of all types of trees"""

    ParentItem = collections.namedtuple("ParentItem", ['name', 'distance'])

    def reinitTree(self, stream=open(os.devnull, 'w')):
        # object initialisations
        self.parent = self.newCommonNamesMapperInstance()
        """self.parent = {..., son: parentItem, ...}"""

        self.species = self.newCommonNamesMapperInstance()
        """self.species = {..., nodeI: [extant nodes that are descendant from nodeI], ...}
        interior nodes are not included in the list of descendants"""

        self.fileName = self.newCommonNamesMapperInstance()
        """self.fileName = {..., node: "fileName", ...} with 'fileName' the name
        given in genes or ancGenes files (eg Homo.sapiens for node 'Homo
        sapiens'"""

        self.dicLinks = self.newCommonNamesMapperInstance()
        """self.dicLinks = {..., node1: {..., node2: [list of nodes on the
        evolutive path between node1 and node 2], ...}, ...}
        Usage:
        self.dicLinks[node1][node2] = [list of node on the
        evolutive path between node1 and node 2]"""

        self.dicParents = self.newCommonNamesMapperInstance()
        """self.dicParents = {..., node1: {..., node2: [list of nodes that are
        descendants of the LCA(node1, node2)], ...}, ...}"""
        # NOT AT ALL. It contains the LCA(node1, node2)

        self.allDescendants = self.newCommonNamesMapperInstance()
        """self.allDescendants = {..., nodeI: [nodes that are descendant from nodeI], ...}
        interior nodes are included in the list of descendants"""
        
        self.tmpS = []
        self.tmpA = []

        def recInitialize(node, stream=open(os.devnull, 'w')):
            """analysing process of the tree"""
            print(".", file=stream)
            self.dicLinks.setdefault(node, self.newCommonNamesMapperInstance())
            # each link between a node and itself is a list containing the node itself
            self.dicLinks.get(node).setdefault(node, [node])
            self.dicParents.setdefault(node, self.newCommonNamesMapperInstance())
            self.dicParents.get(node).setdefault(node, node)
            family = [node]  # All descendants of this node.
            if node in self.items:
                # interior node (not a leaf of the tree)
                self.fileName.setdefault(node, str(node).replace(' ', '_'))
                # List of descendant species
                s = []
                # list of subFamilies
                lsf = []
                self.tmpA.append(node)
                for (son, bLength) in self.items.get(node):
                    self.parent.setdefault(son, PhylogeneticTree.ParentItem(node, bLength))
                    subFamily = recInitialize(son, stream)
                    s.extend(self.species.get(son))
                    # we go back up
                    family.extend(subFamily)
                    lsf.append(subFamily)
                    for e in subFamily:
                        self.dicParents.get(e).setdefault(node, node)
                        self.dicParents.get(node).setdefault(e, node)
                        self.dicLinks.get(e).setdefault(node, self.dicLinks.get(e).get(son) + [node])
                        self.dicLinks.get(node).setdefault(e, [node] + self.dicLinks.get(son).get(e))
                self.species.setdefault(node, frozenset(s))
                # phylogenetic relationships
                for (sf1, sf2) in itertools.combinations(lsf, 2):
                    for e1 in sf1:
                        for e2 in sf2:
                            # node is the LCA of e1, e2
                            self.dicParents.get(e1).setdefault(e2, node)
                            self.dicParents.get(e2).setdefault(e1, node)
                            self.dicLinks.get(e1).setdefault(e2, self.dicLinks.get(e1).get(node) + self.dicLinks.get(node).get(e2)[1:])
                            self.dicLinks.get(e2).setdefault(e1, self.dicLinks.get(e2).get(node) + self.dicLinks.get(node).get(e1)[1:])
            else:
                # leaf of the tree
                self.fileName.setdefault(node, str(node).replace(' ', '.'))
                self.species.setdefault(node, frozenset([node]))
                self.tmpS.append(node)
            self.allDescendants.setdefault(node, frozenset(family))
            return family

        print("Data analysis ", end=' ', file=stream)
        # filling
        recInitialize(self.root, stream=stream)

        self.listSpecies = frozenset(self.species).difference(self.items)
        """list of extant species"""

        assert frozenset(self.tmpS) == self.listSpecies

        self.listAncestr = frozenset(self.items)
        """list of ancestral species"""
        
        assert frozenset(self.tmpA) == self.listAncestr

        # post-analysis
        self.outgroupSpecies = self.newCommonNamesMapperInstance()
        # FIXME: bad name, this is a dict that take a name as keys and return
        # the set of outgroup species names.
        self.indNames = self.newCommonNamesMapperInstance()
        """This is a dict that takes a name as key and return the
        branch ID."""
        # setdefault: This method returns the key value available in the
        # dictionary and if given key is not available then it will return
        # provided default value.
        self.indBranches = self.newCommonNamesMapperInstance()
        # This loop defines the order of iterations in the tree
        self.allNames = self.tmpA + self.tmpS
        for (i, e) in enumerate(self.allNames):
            self.outgroupSpecies.setdefault(e, self.listSpecies.difference(self.species.get(e)))
            self.indNames.setdefault(e, i)
            if e == self.root:
                assert i == 0
            elif e != self.root:
                # each branch correspond to the son species, the root species
                # has no associated branch
                self.indBranches.setdefault(e, i-1)
            self.officialName[i] = e
        # (e)numerateItems, all nodes excepted the root using indNames
        self.numItems = [[]] * len(self.allNames)
        for a in self.listAncestr:
            self.numItems[self.indNames.get(a)] = [(self.indNames.get(son), l) for (son, l) in self.items.get(a)]
        # (e)numerateParents, all nodes excepted the leaves using indNames
        self.numParent = [None] * len(self.allNames)
        for (e, (par, l)) in self.parent.items():
            self.numParent[self.indNames.get(e)] = (self.indNames.get(par), l)
        # self.numParent[self.indNames.get('self.root')] = None
        assert set(self.indBranches.values()) == set(range(len(self.numParent)-1))
        # official names / common names
        tmp = collections.defaultdict(list)
        for (curr, off) in self.officialName.items():
            tmp[off].append(curr)
        self.commonNames = self.newCommonNamesMapperInstance()
        self.commonNames.update(tmp)
        self.dicGenes = {}
        self.dicGenomes = self.newCommonNamesMapperInstance()

        print(" OK", file=stream)

    def __init__(self, file, skipInit=False, stream=open(os.devnull, 'w')):
        if type(file) == tuple:
            print("Creation of the phylogenetic tree ...", end=' ', file=stream)
            (self.items, self.root, self.officialName) = file
        else:
            print("Loading phylogenetic tree %s ..." % file, end=' ', file=stream)
            self.officialName = {}
            self.items = self.newCommonNamesMapperInstance()
            self.ages = self.newCommonNamesMapperInstance()
            self.lstEsp2X = set()
            self.lstEsp6X = set()
            self.lstEspFull = set()
            # name and instance of file
            f = myFile.openFile(file, 'r')
            try:
                self.name = f.name
            except AttributeError:
                self.name = file
            f = myFile.firstLineBuffer(f)
            if (';' in f.firstLine) or ('(' in f.firstLine):
                self.__loadFromNewick__(' '.join(f).replace('\n', '') + " ;")
            else:
                self.__loadFromMyFormat__(f)
            f.close()
            if not skipInit:
                self.reinitTree()
            else:
                print("OK", file=stream)

    def __loadFromMyFormat__(self, f):
        """load a phylogenetic tree into the phylTree format (with tabulations)"""

        def loadFile():
            """loading process"""
            lines = []
            for line in f:
                # the final "\n" is removed and we cut owing to the "\t"
                l = line.split('\t')
                # indentation level
                indent = 0
                for x in l:
                    if x != '':
                        break
                    indent += 1
                # the triplet (indentation,names,age) is recorded
                assert l[indent] == x
                names = [xx.strip() for xx in x.split('|')]
                if len(l) > (indent+1):
                    try:
                        age = int(l[indent+1])
                    except ValueError:
                        age = float(l[indent+1])
                else:
                    # NOT A NICE IDEA. TODO. change to None?
                    age = 0
                lines.append((indent, names, age))
            lines.reverse()
            return lines

        def recLoad(indent, lines):
            """lines analysis process
            lines = [..., (indent, names, age), ..., (0, namesOfRoot, age)]"""
            # we continue only if the next line is indented as expected
            if len(lines) == 0 or lines[-1][0] != indent:
                return None
            # the line is loaded
            currLine = lines.pop()
            # all the children sub-trees are loaded, as long as it is possible
            children = []
            while True:
                tmp = recLoad(indent + 1, lines)
                if tmp is None:
                    break
                # record (name_of_child, evolution_time)
                children.append((tmp, currLine[2] - self.ages.get(tmp)))
            # current node
            (currIndent, currNames, currAge) = currLine
            
            #def parse_names(self, currNames):
            name = currNames[0]
            if len(children) == 0:
                if name[0] == SYMBOL6X:
                    currLine[1][0] = name = currNames[0][1:]
                    self.lstEsp6X.add(name)
                elif name[0] == SYMBOL2X:
                    currLine[1][0] = name = currNames[0][1:]
                    self.lstEsp2X.add(name)
                else:
                    name = currNames[0]
                    self.lstEspFull.add(name)
            for s in currLine[1]:
                self.officialName[s] = name

            # FIXE: also return anc node if only one extant species
            if len(children) == 1:
                self.items.setdefault(name, children)
            # several children, they are recorded
            elif len(children) > 1:
                self.items.setdefault(name, children)
            # standard informations
            self.ages.setdefault(name, currLine[2])
            return name

        lines = loadFile()
        self.root = recLoad(0, lines)

    def __loadFromNewick__(self, s):
        """tree in the Newick format (between parentheses)"""

        s = " ".join(s.split()) # remove any consecutive whitespaces
        
        #FIXME
        def readStr(nb):
            """read the nb next characters of the tree"""
            ret = s[self.pos:self.pos+nb]
            self.pos += nb
            return ret

        def keepWhile(car):
            """remove forbidden characters"""
            x = self.pos
            while s[x] in car:
                x += 1
            return readStr(x-self.pos)

        def keepUntil(car):
            """keep the authorised characters"""
            x = self.pos
            while s[x] not in car:
                x += 1
            return readStr(x-self.pos)

        def readTree():
            """read the tree in the form of text informations"""

            keepWhile(' \t')

            # Guillaume 2017-10-03: allow spaces in labels
            # --------------------------------------------
            # On the contrary to the `skbio` newick parser [implementation](
            # http://scikit-bio.org/docs/latest/generated/skbio.io.format.newick.html),
            # here you cannot escape node labels with single quotes (meaning 
            # the single quotes become part of the node label. Similarly to
            # `Ete3`).
            # Therefore it wasn't originally possible to have spaces in labels,
            # but I allowed it here, on one *condition*: multiple consecutive
            # whitespaces in labels are **converted into a single space**.
            
            if s[self.pos] == '(':
                children = []
                # '(' the first time, then some ',' untill the final ')'
                while readStr(1) != ')':
                    children.append(readTree())
                    keepWhile(' \t')
                keepWhile(' \t')
                # the result is the list of children + the name
                elt = (children, keepUntil("),:;[\t"))
            else:
                # the result is a name, for leaves, no children
                elt = ([], keepUntil("),:;[\t"))
            keepWhile(' \t')

            # possibly a non-null branch length
            if s[self.pos] == ':':
                # ":"
                readStr(1)
                length = keepWhile("0123456789.eE-")
                try:
                    length = int(length)
                except ValueError:
                    length = float(length)
            else:
                length = 0

            keepWhile(' \t')

            # possibly informations between brackets
            if s[self.pos] == '[':
                # "["
                readStr(1)
                info = keepUntil("]")
                info = dict(x.split("=") for x in info.split(":") if "=" in x)
                # "]"
                readStr(1)
            else:
                info = {}

            keepWhile(' \t')

            return (elt, length, info)

        # Write integer in base 26.
        def format_uniq_suffix(n):
            r = [n]
            while r[0] >= 26:
                r.insert(1, r[0] % 26)
                r[0] //= 26
            return ''.join(chr(65+k) for k in r)

        def make_uniq_name(name):  # Could be a class method.
            if not name: name = 'UNNAMED'
            n = 0
            uniqname = name
            while uniqname in self.officialName:
                uniqname = '%s_%s' %(name, format_uniq_suffix(n))
                n += 1
            return uniqname

        def storeTree(data):
            """fill the variables of a GenericTree"""
            ((children, name), length, info) = data
            
            #name = self.parse_names(currNames)
            currNames = [xx.strip() for xx in name.split('|')]
            name = currNames[0]
            if len(children) == 0:
                if name[0] == SYMBOL6X:
                    name = currNames[0][1:]
                    self.lstEsp6X.add(name)
                elif name[0] == SYMBOL2X:
                    name = currNames[0][1:]
                    self.lstEsp2X.add(name)
                else:
                    name = currNames[0]
                    self.lstEspFull.add(name)

            name = make_uniq_name(name)
            currNames[0] = name
            for s in currNames:
                self.officialName[make_uniq_name(s)] = name
            
            self.info[name] = info

            tmp_ages = []
            if len(children) > 0:
                items = []
                for arbre in children:
                    child, dist = storeTree(arbre)
                    items.append((child, dist))
                    tmp_ages.append(dist + self.ages[child])

                self.items.setdefault(name, items)
            self.ages.setdefault(name, max(tmp_ages, default=0))

            self.root = name
            return (name, length)

        #FIXME add the age of ancestors for an easier conversion into myFormat
        self.pos = 0
        data = readTree()
        self.pos = 0
        self.info = {}
        storeTree(data)

    def to_ete3(self, nosinglechild=False):
        """Convert to Ete3 tree object"""
        import ete3
        # TODO: do not import ete3 here, just try to use it and raise error
        tree = ete3.Tree(name=self.root)
        tree.add_features(treename=self.name)
        current_nodes = [tree]
        while current_nodes:
            current = current_nodes.pop()
            current.add_features(age=self.ages.get(current.name),
                                 commonNames=self.commonNames.get(current.name),
                                 fileName=self.fileName.get(current.name),
                                 indBranch=self.indBranches.get(current.name),
                                 indName=self.indNames.get(current.name),
                                 )
                                 #officialNames=self.officialNames.get(current.name),
            childlist = self.items.get(current.name, [])
            if not childlist:
                current.add_features(Esp2X=(current.name in self.lstEsp2X),
                                     Esp6X=(current.name in self.lstEsp6X),
                                     EspFull=(current.name in self.lstEspFull))
                # dicGenes
                # dicGenomes
            while nosinglechild and (len(childlist) == 1):
                current.name = childlist[0][0]
                current.dist += childlist[0][1]
                childlist = self.items.get(current.name, [])

            for child, dist in childlist:
                current_nodes.append(current.add_child(name=child, dist=dist))
        return tree

    def draw(self, images=None):
        """Draw tree in a gui window using Ete3"""
        import ete3
        ptree = self.to_ete3()
        ts = ete3.TreeStyle()
        ts.scale = 2
        ts.show_branch_length = True
        nameface = ete3.faces.AttrFace("name", fsize=12)
        def mylayout(node):
            if node.is_leaf():
                # add species image
                pass
            else:
                ete3.faces.add_face_to_node(nameface, node, column=0,
                                            position='branch-top')
        ptree.show(layout=mylayout, tree_style=ts, name='Species tree')


    def lastCommonAncestor(self, species):
        """return the name of the last common ancestor of several species"""
        anc = species[0]
        for e in species[1:]:
            anc = self.dicParents[anc][e]
            if anc == self.root:
                return self.root
        return anc

    def isChildOf(self, child, parent):
        """assess if 'child' is a **descendant** of 'parent'"""
        return self.dicParents[child][parent] == self.officialName[parent]

    #FIXME
    def compact(self, maxlength=1e-4):
        def do(node):
            if node not in self.items:
                return
            newitems = []
            for x in self.items[node]:
                do(x[0])
                if (x[1] >= maxlength) or (x[0] not in self.items):
                    newitems.append(x)
                else:
                    print("Removing node", x[0], file=sys.stderr)
                    newitems.extend(self.items.pop(x[0]))
            self.items[node] = newitems
        do(self.root)
        self.reinitTree()

    def getTargets(self, s):
        """return the list of species to be compared
        ancestors must be prefixed with '.' to be taken into account as genomes
        otherwise, that is their descendant species that are used"""

        # are we exclusively on the common ancestors
        allIntermediates = (s[-1] == "+")
        if allIntermediates:
            s = s[:-1]

        # list of extant species
        listSpecies = set()
        listAncestors = set()
        for x in s.split(','):
            if x[0] != '.':
                listSpecies.update(self.species[x])
            else:
                listSpecies.add(x[1:])
                listAncestors.add(x[1:])

        # list of ancestors
        for (e1, e2) in itertools.combinations(listSpecies, 2):
            if allIntermediates:
                listAncestors.update(self.dicLinks[e1][e2][1:-1])
            else:
                listAncestors.add(self.dicParents[e1][e2])
            # if e1 or e2 is already an ancestor

        return (listSpecies, listAncestors)

    def getTargetsSpec(self, target):
        """return the list of species pointed out by target
          =esp
          +anc [children]
          /anc [outgroup]
          _** in order to remove instead of adding"""
        lesp = set()
        for x in target.split(","):
            if x.startswith("_"):
                if x.startswith("_/"):
                    lesp.difference_update(self.outgroupSpecies[x[2:]])
                elif x.startswith("_="):
                    lesp.discard(x[2:])
                else:
                    lesp.difference_update(self.species[x[1:]])
            else:
                if x.startswith("/"):
                    lesp.update(self.outgroupSpecies[x[1:]])
                elif x.startswith("="):
                    lesp.add(x[1:])
                else:
                    lesp.update(self.species[x])
        return lesp

    def getTargetsAnc(self, target):
        """return the list of ancestors pointed out by target
          =anc
          +anc [children]
          -anc [parents]
          /anc [parents+outgroups]
          _** in order to remove instead of adding"""
        lanc = set()
        for x in target.split(","):
            if x.startswith("_"):
                if x.startswith("_/"):
                    lanc.difference_update(e for e in self.listAncestr if not self.isChildOf(e, x[2:]))
                elif x.startswith("_="):
                    lanc.discard(x[2:])
                elif x.startswith("_\\"):
                    lanc.difference_update(e for e in self.listAncestr if self.isChildOf(x[2:], e))
                else:
                    lanc.difference_update(e for e in self.listAncestr if self.isChildOf(e, x[1:]))
            else:
                if x.startswith("/"):
                    lanc.update(e for e in self.listAncestr if not self.isChildOf(e, x[1:]))
                elif x.startswith("="):
                    lanc.add(x[1:])
                elif x.startswith("\\"):
                    lanc.update(e for e in self.listAncestr if self.isChildOf(x[1:], e))
                else:
                    lanc.update(e for e in self.listAncestr if self.isChildOf(e, x))
        return lanc

    def getSubTree(self, goodSpecies, rootAnc=None):
        """return the structure of the sub-tree in which are exclusively the
        chosen species"""
        goodAnc = set(self.dicParents[e1][e2] for (e1, e2) in itertools.combinations(goodSpecies, 2))
        newtree = collections.defaultdict(list)
        from . import myMaths

        def do(node):
            if node in goodAnc:
                for (x, _) in self.items[node]:
                    newtree[node].extend(do(x))
                return [node]
            elif node in self.items: # i.e, not an extant species.
                # TODO: itertools.chain?
                return myMaths.flatten([do(x) for (x, _) in self.items[node]])
            elif node in goodSpecies:
                return [node]
            else:
                return []
        root = do(self.root if rootAnc is None else rootAnc)
        assert len(root) == 1
        return (root[0], newtree)

    def getSubAncTree(self, goodNodes, rootAnc=None):
        """return the structure of the sub-tree in which are exclusively the
        chosen nodes.
        
        It's different from `getSubTree` because if it will not omit the most
        recent nodes given, where not descendants leaves were listed.
        For example:
            
            >>> phyltree.getSubTree(['Lagomorpha', 'Rodentia'])
            ('Glires', defaultdict(list, {'Glires': []}))

        but

            >>> phyltree.getSubAncTree(['Lagomorpha', 'Rodentia'])
            ('Glires', defaultdict(list, {'Glires': ['Lagomorpha', 'Rodentia']}))
        """
        goodAnc = set(self.dicParents[e1][e2] for (e1, e2) in itertools.combinations(goodNodes, 2))
        goodLeaves = set(goodNodes) - goodAnc
        newtree = collections.defaultdict(list)
        from . import myMaths

        def do(node):
            if node in goodAnc:
                for (x, _) in self.items[node]:
                    newtree[node].extend(do(x))
                return [node]
            elif node in goodLeaves:
                newtree[node] = []
                return [node]
            elif node in self.items: # i.e, not an extant species.
                # TODO: itertools.chain?
                return myMaths.flatten([do(x) for (x, _) in self.items[node]])
            else:
                return []
        root = do(self.root if rootAnc is None else rootAnc)
        assert len(root) == 1
        return (root[0], newtree)

    def calcWeightedValue(self, values, notdefined, resultNode):
        """calculate the values for the nodes of the tree
         - 'values' represents some defined values (nodes or leaves)
         - 'notdefined' is the value to be returned if no result
         - 'resultNode' indicates at which ancestral node we want to find the result"""
        import numpy
        n = len(self.allNames)
        # by default the results will be wil be "notdefined"
        matriceA = numpy.identity(n)
        matriceB = numpy.empty((n,))
        matriceB.fill(notdefined)

        def recBuild(ianc, parent, length):
            """recursively build the matrix"""
            anc = self.allNames[ianc]
            #FIXME recursive calls: children are selected OK (supposing that we are it ourself)
            items = [(e, p) for (e, p) in self.numItems[ianc] if recBuild(e, ianc, p)]
            # parent if available
            if parent is not None:
                items.append((parent, length))

            if anc in values:
                # if we have a value, "x = val" is assigned
                matriceB[ianc] = values.get(anc)
                return True

            elif len(items) >= 2:
                # if it has enough neighbours, the equation is written
                s = 0.
                for (e, p) in items:
                    p = 1./max(p, 0.00001)
                    matriceA[ianc][e] = p
                    s += p
                matriceA[ianc][ianc] = -s
                matriceB[ianc] = 0
                return True
            else:
                return False

        # matrix building
        if len(values) == 0:
            return matriceB

        rootNode = self.indNames[self.lastCommonAncestor(list(values.keys()))]
        recBuild(rootNode, None, 0)
        # equation solving
        res = numpy.linalg.solve(matriceA, matriceB)

        def searchBestValue(node, origin, tripLength):
            """search the neirest value around 'node' knowing that we come from
            'origin'"""
            # Does node have a defined value
            val = res[node]
            if val != notdefined:
                return (tripLength, node, val)

            # Paths to go through
            test = self.numItems[node][:]
            par = self.numParent[node]
            if par is not None:
                test.append(par)
            # Iteration
            best = (1e300, node, notdefined)
            for (e, l) in test:
                # In order to not loop
                if e == origin:
                    continue
                # Break
                if best[0] < tripLength+l:
                    continue
                tmp = searchBestValue(e, node, tripLength+l)
                if tmp[0] < best[0]:
                    best = tmp
            return best

        if resultNode is None:
            return res
        resultNode = self.indNames[resultNode]
        best = searchBestValue(resultNode, None, 0)
        return (best[0], self.allNames[best[1]], best[2])

    def newCommonNamesMapperInstance(self):
        """return a dict that uses internally official names of taxons
         but that can be used with common names"""
        dsi = dict.__setitem__
        dgi = dict.__getitem__

        class commonNamesMapper(dict):

            def __getitem__(d, name):
                if name in self.officialName:
                    return dgi(d, self.officialName[name])
                else:
                    return dgi(d, name)

            def __setitem__(d, name, value):
                if name in self.officialName:
                    dsi(d, self.officialName[name], value)
                else:
                    dsi(d, name, value)

            # Because it's recursive, use a non-ambiguous name for this method,
            # since 'to_dict' for instance, is already a method of
            # pandas.DataFrames.
            def common_names_mapper_2_dict(self):
                """Recursively converts the commonNamesMapper instance to a
                dictionary. (to allow pickling for example)"""
                as_dict = {}
                for key, value in self.items():
                    if hasattr(value, 'common_names_mapper_2_dict'):
                        value = value.common_names_mapper_2_dict()
                    as_dict[key] = value
                return as_dict

        return commonNamesMapper()

    def loadAllSpeciesSince(self, ancestr, template, **args):
        """load all the species that come from an ancestor"""
        if ancestr in self.species:
            l = self.species[ancestr]
        else:
            l = self.listSpecies
        self.loadSpeciesFromList(l, template, **args)

    def loadAllSpeciesBefore(self, ancestr, template, **args):
        """load all the outgroup species of an ancestor"""
        if ancestr in self.outgroupSpecies:
            l = self.outgroupSpecies[ancestr]
        else:
            l = self.listSpecies
        self.loadSpeciesFromList(l, template, **args)

    def loadSpeciesFromList(self, lst, template, storeGenomes=True):
        """load all the species of a list"""

        from . import myGenomes

        for esp in lst:
            esp = self.officialName[esp]
            g = myGenomes.Genome(template % self.fileName[esp])
            if storeGenomes:
                self.dicGenomes[esp] = g
            for (x, (c, i)) in g.dicGenes.items():
                #self.dicGenes[x] = (esp, c, i)
                self.dicGenes[x] = GeneSpeciesPosition(esp, c, i)

    def findFamilyComposition(self, fam):
        """return the number of genes in each species for a given family"""
        score = self.newCommonNamesMapperInstance()
        score.update(dict.fromkeys(self.officialName, []))
        for g in fam:
            if g in self.dicGenes:
                (e, _, _) = self.dicGenes[g]
                score[e].append(g)
        return score

    def reroot(self, node, useOutgroups, newname=0):
        """topple over ("basculer" in french) the tree around a node"""
        self.tmpItems = self.items.copy()
        if useOutgroups:
            anc = node
            while anc in self.parent:
                (par, l) = self.parent[anc]
                self.tmpItems[anc] = self.tmpItems[anc] + [(par, l)]
                self.tmpItems[par] = [x for x in self.tmpItems[par] if x[0] != anc]
                anc = par
        self.tmpItems[newname] = self.tmpItems[node]

    def printBranchName(self, child, stream=sys.stderr):
        # To print branch name to sys.stderr
        # Please be aware that changing this breaks the analysis of the stderr by myScore.effectiveValuesFromSimulationLogErr()
        (parent, bLength) = self.parent[child]
        foo = "# Branch %s -> %s (%s My) #" % (parent, child, bLength)
        bar = "#" * len(foo)
        print(bar, file=stream)
        print(foo, file=stream)
        print(bar, file=stream)

    # FIXME, ensure that no two species has the same acronym!!
    def speciesAcronym(self, node):
        if node in self.listSpecies:
            # First char of each word with a upper char for the first
            speciesAcronym = [word for word in node.split(' ')]
            speciesAcronym = [speciesAcronym[0][0].upper()] + [char[0].lower() for char in speciesAcronym[1:]]
        elif node in self.listAncestr:
            # The two first chars
            speciesAcronym = node[:2]
        return ''.join(speciesAcronym)
