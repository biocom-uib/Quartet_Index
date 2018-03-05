# qi_6_value.py
# (c) 2018 Tomas M. Coronado, Arnau Mir, Francesc Rossello, Gabriel Valiente
# 
# input:	phylogenetic tree T in Newick format
# output:	value of QI(T)

from ete3 import Tree
from itertools import combinations
import sys

in_file = sys.argv[1]
T = Tree(in_file, format=1)

s4 = 4
b4 = 3
q31 = 2
q211 = 1
k4 = 0

def qi_value(t, a, b, c, d):
    ab = t.get_common_ancestor(str(a), str(b))
    ac = t.get_common_ancestor(str(a), str(c))
    ad = t.get_common_ancestor(str(a), str(d))
    bc = t.get_common_ancestor(str(b), str(c))
    bd = t.get_common_ancestor(str(b), str(d))
    cd = t.get_common_ancestor(str(c), str(d))
    d = dict()
    for w in [ab, ac, ad, bc, bd, cd]:
    	if w in d:
    		d[w] += 1
    	else:
    		d[w] = 1
    s = str(sorted(d.values()))
    if s == '[6]': # S_4
    	return s4
    elif s == '[1, 5]': # Q_{2,1,1}
    	return q211
    elif s == '[3, 3]': # Q_{3,1}
    	return q31
    elif s == '[1, 1, 4]': # B_4
    	return b4
    elif s == '[1, 2, 3]': # K_4
    	return k4
    else:
    	print "error", s

def qi(t):
	sum = 0
	taxa = [node.name for node in t.get_leaves()]
	for [a, b, c, d] in combinations(taxa, 4):
		sum += qi_value(t, a, b, c, d)
	return sum

print qi(T)
