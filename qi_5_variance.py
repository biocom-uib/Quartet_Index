# qi_5_variance.py
# (c) 2018 Tomas Martinez, Arnau Mir, Francesc Rossello, Gabriel Valiente
# 
# input:	file "Var_n.txt", real alpha, real gamma, integer n
# output:	value of Var_{M_{alpha,gamma}}(QI_n)

from sympy import *

a, g, n = symbols("a g n") # alpha, gamma, n

import sys

aa = sys.argv[1]
gg = sys.argv[2]
nn = sys.argv[3]

with open("Var_n.txt", "r") as f:
	v = f.read().rstrip()
v = sympify(v)

print v.subs(a, aa).subs(g, gg).subs(n, nn)

