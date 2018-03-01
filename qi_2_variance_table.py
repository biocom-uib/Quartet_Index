# qi_2_variance_table.py
# (c) 2018 Tomas Martinez, Arnau Mir, Francesc Rossello, Gabriel Valiente
# 
# input:	files "P2.txt" ... "P8.txt"
# output:	file "variance_table.txt"

from sympy import *
from itertools import combinations, permutations
from math import factorial
from ete3 import Tree

a, g = symbols("a g") # alpha, gamma

q0 = 0; q1, q2, q3, q4 = symbols("q1 q2 q3 q4")
q = [q0, q1, q2, q3, q4]

def binom(n, k):
    if n < k:
        return 0
    else:
        return factorial(n)/factorial(k)/factorial(n - k)

def read_data(D, f):
	for line in open(f):
		line = line.rstrip().split(' ', 1)
		shape, formula = line[0], line[1]
		D[shape] = sympify(formula)
	print "read", len(D), "probabilities"

P = dict()
P[2] = dict(); read_data(P[2], "P2.txt")
P[3] = dict(); read_data(P[3], "P3.txt")
P[4] = dict(); read_data(P[4], "P4.txt")
P[5] = dict(); read_data(P[5], "P5.txt")
P[6] = dict(); read_data(P[6], "P6.txt")
P[7] = dict(); read_data(P[7], "P7.txt")
P[8] = dict(); read_data(P[8], "P8.txt")

S4 = "(*,*,*,*);"
B4 = "((*,*),(*,*));"
Q31 = "(*,(*,*,*));"
Q211 = "(*,*,(*,*));"
K4 = "(*,(*,(*,*)));"
T4 = [K4, Q211, Q31, B4, S4]

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

def shape(t):
	tt = code_to_tree(tree_to_code(t))
	for leaf in tt:
		leaf.name = "*"
	return tt.write(format=9)

out = open("variance_table.txt", "w+")
S = dict()
for i in range(1, 5):
	S[i] = dict()
	for j in range(1, 5):
		S[i][j] = dict()
		for k in range(5, 9):
			sum = 0
			for s in P[k].keys():
				t = Tree(s, format=9)
				num = 0
				for leaf in t:
					leaf.name = str(num)
					num += 1
				count = 0
				for Q in combinations(range(k), 4):
					Q = list(map(str, Q))
					TQ = t.copy()
					TQ.prune(Q)
					for Qp in combinations(range(k), 4):
						Qp = list(map(str, Qp))
						TQp = t.copy()
						TQp.prune(Qp)
						if len(set(Q).intersection(set(Qp))) == 8 - k and shape(TQ) == T4[i] and shape(TQp) == T4[j]:
							count += 1
				sum += simplify(q[i]*q[j]*count*P[k][s])
			S[i][j][k] = simplify(sum)
			print "round", i + j + (k - 4) - 2, "of", 4*4*4
			print >>out, i, j, k, S[i][j][k]

