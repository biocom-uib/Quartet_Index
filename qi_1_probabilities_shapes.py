# qi_1_probabilities_shapes.py
# (c) 2018 Tomas Martinez, Arnau Mir, Francesc Rossello, Gabriel Valiente
# 
# input:	---
# output:	files "P2.txt" ... "P8.txt"

from ete3 import Tree
from sympy import *

a, g = symbols("a g") # alpha, gamma

# Zhang - 2008 - Discovering frequent agreement subtrees from phylogenetic data
def canonical_newick_string(s):
	T = Tree(s, format=9)
	for node in T.traverse("postorder"):
		# taxa = sorted([leaf.name for leaf in node.get_leaves()])
		if node.is_leaf():
			taxa = node.name
		else:
			taxa = min([child.taxa for child in node.children])
		node.add_features(taxa = taxa)
	for node in T.traverse("postorder"):
		if not node.is_leaf():
			node.children.sort(key=lambda x: x.taxa)
	return T.write(format=9)

flatten = lambda l: [item for sublist in l for item in sublist]

def tree_to_code(t):
	for node in t.traverse("postorder"):
		if node.is_leaf():
			code = [1]
		else:
			size = len([desc for desc in node.traverse()])
			code = [size] + flatten(sorted([child.code for child in node.children]))
		node.add_features(code = code)
	return t.code

def newick_to_code(s):
	t = Tree(s, format=9)
	for node in t.traverse("postorder"):
		if node.is_leaf():
			code = [1]
		else:
			size = len([desc for desc in node.traverse()])
			code = [size] + flatten(sorted([child.code for child in node.children]))
		node.add_features(code = code)
	return t.code

def split(tail):
	chunks = []
	chunk = []
	last = -1
	for i in range(len(tail)):
		if tail[i] >= last: # start new chunk
			chunks.append(chunk)
			last = tail[i]
			chunk = [last]
		else: # extend current chunk
			chunk.append(tail[i])
	chunks.append(chunk)
	return chunks[1:]

def code_to_tree(code):
	head, tail = code[0], code[1:]
	t = Tree()
	for chunk in split(tail):
		t.add_child(code_to_tree(chunk))
	return t

def shape(s):
	return code_to_tree(newick_to_code(s)).write(format=100)

def extend(t, P, Q):
	n = len(t)
	s = canonical_newick_string(t.write(format=9))
	# 1: add the leaf to an external arc
	# P[m+1] = factor(((1-a)/(m-a))*P[m])
	leaves = [node for node in t]
	for node in leaves:
		x = node.add_child(name=node.name)
		y = node.add_child(name=n+1)
		z = t.copy()
		# print "case 1", s, z.write(format=9)
		Q[canonical_newick_string(z.write(format=9))] = factor(((1-a)/(n-a))*P[s])
		x.detach()
		y.detach()
	# 2: add the leaf to an internal arc
	# P[m+1] = factor((g/(m-a))*P[m])
	internal = [node for node in t.iter_descendants() if not node.is_leaf()]
	for node in internal:
		node.add_features(name="internal")
		x = t.copy("deepcopy")
		node.del_feature("name")
		y = x.search_nodes(name="internal")[0]
		parent = y.up
		old = y.detach()
		new = Tree()
		new.add_child(old)
		new.add_child(name=n+1)
		parent.add_child(new)
		# print "case 2", s, x.write(format=9)
		Q[canonical_newick_string(x.write(format=9))] = factor((g/(n-a))*P[s])
	# 3: add the leaf to a new root
	# P[m+1] = factor(g/(m-a)*P[m])
	x = t.copy()
	y = Tree()
	y.add_child(x)
	y.add_child(name=n+1)
	# print "case 3", s, y.write(format=9)
	Q[canonical_newick_string(y.write(format=9))] = factor(g/(n-a)*P[s])
	# 4: add the leaf to the root
	# P[m+1] = factor(((len(t.children)-1)*a-g)/(m-a)*P[m])
	if not t.is_leaf():
		x = t.copy()
		x.add_child(name=n+1)
		# print "case 4r", s, x.write(format=9)
		Q[canonical_newick_string(x.write(format=9))] = factor(((len(t.children)-1)*a-g)/(n-a)*P[s])
	# 4: add the leaf to an internal node
	# P[m+1] = factor(((len(node.children)-1)*a-g)/(m-a)*P[m])
	for node in t.iter_descendants():
		if not node.is_leaf():
			y = node.add_child(name=n+1)
			# print "case 4i", s, t.write(format=9)
			Q[canonical_newick_string(t.write(format=9))] = factor(((len(node.children)-2)*a-g)/(n-a)*P[s])
			y.detach()

def store(D, s, f):
	ts = shape(s)
	if ts in D:
		D[ts] = factor(D[ts] + f)
	else:
		D[ts] = f

def fill_taxa(s):
	s = s.replace("(,", "(*,")
	s = s.replace(",,", ",*,")
	s = s.replace(",,", ",*,")
	return s.replace(",)", ",*)")

P2 = dict()
P3 = dict()
P4 = dict()
P5 = dict()
P6 = dict()
P7 = dict()
P8 = dict()

s2 = "(1,2);"
P2[s2] = 1

PS2 = dict()
for s in P2.keys():
	store(PS2, s, P2[s])

out = open("P2.txt", "w")
for s in sorted(PS2.keys()):
	print >>out, fill_taxa(s), PS2[s]

for s in P2.keys():
	t = Tree(s, format=9)
	extend(t, P2, P3)

PS3 = dict()
for s in P3.keys():
	store(PS3, s, P3[s])

out = open("P3.txt", "w")
for s in sorted(PS3.keys()):
	print >>out, fill_taxa(s), PS3[s]

for s in P3.keys():
	t = Tree(s, format=9)
	extend(t, P3, P4)

PS4 = dict()
for s in P4.keys():
	store(PS4, s, P4[s])

out = open("P4.txt", "w")
for s in sorted(PS4.keys()):
	print >>out, fill_taxa(s), PS4[s]

for s in P4.keys():
	t = Tree(s, format=9)
	extend(t, P4, P5)

PS5 = dict()
for s in P5.keys():
	store(PS5, s, P5[s])

out = open("P5.txt", "w")
for s in sorted(PS5.keys()):
	print >>out, fill_taxa(s), PS5[s]

for s in P5.keys():
	t = Tree(s, format=9)
	extend(t, P5, P6)

PS6 = dict()
for s in P6.keys():
	store(PS6, s, P6[s])

out = open("P6.txt", "w")
for s in sorted(PS6.keys()):
	print >>out, fill_taxa(s), PS6[s]

for s in P6.keys():
	t = Tree(s, format=9)
	extend(t, P6, P7)

PS7 = dict()
for s in P7.keys():
	store(PS7, s, P7[s])

out = open("P7.txt", "w")
for s in sorted(PS7.keys()):
	print >>out, fill_taxa(s), PS7[s]

for s in P7.keys():
	t = Tree(s, format=9)
	extend(t, P7, P8)

PS8 = dict()
for s in P8.keys():
	store(PS8, s, P8[s])

out = open("P8.txt", "w")
for s in sorted(PS8.keys()):
	print >>out, fill_taxa(s), PS8[s]

