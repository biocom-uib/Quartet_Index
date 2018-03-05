# qi_6_value_fast.py
# (c) 2018 Tomas M. Coronado, Arnau Mir, Francesc Rossello, Gabriel Valiente
# 
# input:	phylogenetic tree T in Newick format
# output:	value of QI(T)

from ete3 import Tree
from itertools import combinations
import sys

in_file = sys.argv[1]
T = Tree(in_file, format=1)

def combin2(n):
	return n*(n - 1)/2

[q0, q1, q2, q3, q4] = [0, 1, 2, 3, 4]

def qi(t):
	for node in t.traverse("postorder"):
		node.add_features(k=0)
		if node.is_leaf():
			node.k = 1
		for child in node.children:
			node.k += child.k
		node.add_features(E3=0)
		for [v1, v2, v3] in combinations(node.children, 3):
			node.E3 += v1.k*v2.k*v3.k
		node.add_features(E4=0)
		for [v1, v2, v3, v4] in combinations(node.children, 4):
			node.E4 += v1.k*v2.k*v3.k*v4.k
	for node in t.traverse("postorder"):
		node.add_features(G=node.E3)
		for child in node.children:
			if not child.is_leaf():
				node.G += child.E3
	for node in t.traverse("postorder"):
		node.add_features(F1=0)
		for [v1, v2, v3] in combinations(node.children, 3):
			node.F1 += combin2(v1.k)*v2.k*v3.k
			node.F1 += combin2(v2.k)*v1.k*v3.k
			node.F1 += combin2(v3.k)*v1.k*v2.k
		node.add_features(F2=0)
		for [v1, v2] in combinations(node.children, 2):
			node.F2 += v1.k*v2.G + v2.k*v1.G
		node.add_features(F3=0)
		for [v1, v2] in combinations(node.children, 2):
			node.F3 += combin2(v1.k)*combin2(v2.k)
	x = 0
	for node in t.traverse("postorder"):
		if not node.is_leaf():
			x += q1*node.F1 + q2*node.F2 + q3*node.F3 + q4*node.E4
	return x

print qi(T)

