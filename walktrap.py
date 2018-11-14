
"""
Fonctions communes de traitement des communautes
Formats d'entree des donnees:
 - Fonction de score
 - Dictionnaire d'aretes
 - Fichier "node1 node2 weight"
Format de sortie des donnees:
 - Liste pour chaque composante connexe
   - Les noeuds de la composante
   - Les coupures interessantes (score,alpha)
   - Le dendogramme (liste des (alphas, fils, pere))
   - Le dendogramme (instance de WalktrapDendogram)
"""

import sys


def doWalktrap(edges, **args):
	"""Lancer un walktrap a partir d'un graphe.

	Permet de detecter les composantes connexes.
	"""

	from . import myTools
	from . import _walktrap

	print("Computing connected components ...", end=' ', file=sys.stderr)
	# Les composantes connexes
	combin = myTools.myCombinator()
	superdelete = []
	
	for (x,l) in edges.items():
		delete = [y for y in l if l[y] <= 0.]
		for y in delete:
			del l[y]
		l.pop(x, None)
		if len(l) > 0:
			combin.addLink(list(l.keys()) + [x])
		else:
			superdelete.append(x)
	
	res = []
	for x in superdelete:
		del edges[x]
		res.append( ([x], [], None, None) )

	print("Launching walktrap ", end=' ', file=sys.stderr)
	n = len(edges)
	for nodes in combin:
		# Reindexation des noeuds
		indNodes = {}
		for (i,node) in enumerate(nodes):
			indNodes[node] = i
	
		(relevantCuts,dend) = _walktrap.doWalktrap(indNodes, edges, **args)

		# On doit revenir aux noms de noeuds originels
		def translate(x):
			if x < len(nodes):
				return nodes[x]
			else:
				return (x,)
		dend = [(cut,tuple(translate(f) for f in fils),translate(pere)) for (cut,fils,pere) in dend]
		res.append( (nodes, relevantCuts, dend, WalktrapDendogram(dend, nodes)) )
		sys.stderr.write('.')
	print(" OK", file=sys.stderr)

	return res


def clusterWithNb(l, funcScore, funcBestChoice, putLonelyinNone, **kwargs):
	"""Lance un walktrap sur la liste de blocs l, en utilisant func pour le
	calcul du score."""

	(edges,none) = funcScore(l)
	print("input", len(edges), len(none), len(l), file=sys.stderr)

	#print >> sys.stderr, edges
	s = len(edges)
	res = doWalktrap(edges, **kwargs)

	# Chaque composante connexe
	chrOrder = []
	for (nodes,cuts,_,dend) in res:
		if len(cuts) == 0:
			chrOrder.append(nodes)
		else:
			(alpha,score,(clust,lonely)) = funcBestChoice([(alpha,score,dend.cut(alpha)) for (alpha,score) in cuts])
			assert sum(len(x) for x in clust) + len(lonely) == len(nodes)
			chrOrder.extend(clust)
			if putLonelyinNone:
				none.update(lonely)
			else:
				chrOrder.append(lonely)

	print("output", len(edges), len(none), [len(x) for x in chrOrder], file=sys.stderr)
	assert len(l) == (len(none) + sum(len(x) for x in chrOrder)), (len(l), len(none), sum(len(x) for x in chrOrder), [len(x) for x in chrOrder])
	return (chrOrder,none)


def applyMultipleClust(ll, lfuncScore, lfuncBestChoice, putLonelyinNone, **kwargs):
	"""Applique succesivement plusieurs clusterWithNb sur chaque liste de blocs
	dans ll."""
	res = []
	none = []
	for l in ll:
		for (funcScore, funcBestChoice) in zip(lfuncScore, lfuncBestChoice):
			(l,n) = clusterWithNb(l, funcScore, funcBestChoice, putLonelyinNone, **kwargs)
			none.extend(n)
		res.extend(l)
	assert sum(len(l) for l in ll) == (len(none) + sum(len(l) for l in res))
	return (res,none)
	return res + [[x] for x in none]



def loadWalktrapOutput(f):
	"""Chargement d'un fichier de resultat de walktrap"""
	
	# On charge les fusions
	allMerges = []
	lstFils = {}
	for line in f:
		if line == "\n":
			break
		
		l = line.split(':')
		scale = float(l[0])
		l = l[1].split('-->')
		
		allMerges.append( (scale,tuple([int(x) for x in l[0].split('+')]),int(l[1])) )

	allMerges.sort( reverse = True )

	lstCoup = []
	for line in f:
		try:
			# On extrait les lignes "alpha relevance"
			c = line.split()
			lstCoup.append( (float(c[0]),float(c[1])) )
		except ValueError:
			pass
	
	return (lstCoup, allMerges)


class WalktrapDendogram:
	"""Le dendogramme resultat, que l'on peut couper a un niveau pour recuperer les classes"""

	def __init__(self, lstMerges, lstNodes):

		self.allMerges = lstMerges
		self.allMerges.sort(reverse = True)
		
		self.lstFils = {}
		for (_,fils,pere) in self.allMerges:
			self.lstFils[pere] = fils
		
		self.lstAll = lstNodes

	def cut(self, scale):
		"""Renvoie (la liste des clusters, les noeuds en dehors des clusters)"""

		# On extrait les communautes correspondantes
		lstClusters = []
		fathersAlreadySeen = set()
		nodesNotSeen = set(self.lstAll)
		for (s,_,pere) in self.allMerges:
			if (s < scale) and (pere not in fathersAlreadySeen):
				cluster = []
				todo = [pere]
				while len(todo) > 0:
					father = todo.pop()
					if father in self.lstFils:
						fathersAlreadySeen.add(father)
						todo.extend(self.lstFils[father])
					else:
						nodesNotSeen.discard(father)
						cluster.append(father)
				lstClusters.append( cluster )
		return (lstClusters, list(nodesNotSeen))


def askPartitionChoice(dend, cuts):
	"""Demande a l'utilisateur quelle partition choisir"""

	from . import myMaths

	def mystr(xxx_todo_changeme ):
		(alpha,relevance,(clusters,lonely)) = xxx_todo_changeme
		return "alpha=%f relevance=%f clusters=%d size=%d lonely=%d sizes={%s}" % \
			(alpha,relevance,len(clusters),sum([len(c) for c in clusters]),len(lonely),myMaths.myStats.txtSummary([len(c) for c in clusters]))

	res = [(alpha,relevance,dend.cut(alpha)) for (alpha,relevance) in cuts]
	# Le choix par defaut
	if len(res) == 1:
		print("1 seule possibilite", file=sys.stderr)
		x = 0
	else:
		# Si on peut, on propose a l'utilisateur de choisir
		for x in res:
			print("> " + mystr(x), file=sys.stderr)
		import os
		if os.isatty(sys.stdin.fileno()):
			print("Aucune entree utilisateur", file=sys.stderr)
			x = 0
		else:
			while True:
				try:
					print("Choix ? ", end=' ', file=sys.stderr)
					x = int(input())
					break
				except ValueError:
					pass
				except IOError:
					x = 0
					break
				except EOFError:
					x = 0
					break
	print("Choix de " + mystr(res[x]), file=sys.stderr)
	return res[x]

