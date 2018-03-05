# qi_3_variance_alpha_gamma.py
# (c) 2018 Tomas M. Coronado, Arnau Mir, Francesc Rossello, Gabriel Valiente
# 
# input:	file "variance_table.txt"
# output:	file  "Var_2_64.txt" with n, Var_{M_{alpha,gamma}}(QI_n) for n = 2 ... 64
#           as a function of alpha, gamma

from sympy import *
from math import factorial

a, g = symbols("a g") # alpha, gamma

q0 = 0
q1, q2, q3, q4 = symbols("q1 q2 q3 q4")
q = [q0, q1, q2, q3, q4]

def binom(n, k):
    if n < k:
        return 0
    else:
        return factorial(n)/factorial(k)/factorial(n - k)

S4 = "(*,*,*,*);"
B4 = "((*,*),(*,*));"
Q31 = "(*,(*,*,*));"
Q211 = "(*,*,(*,*));"
K4 = "(*,(*,(*,*)));"

T4 = [K4, Q211, Q31, B4, S4]

P4 = dict()
P4[S4] = simplify("(a - g)*(2*a - g)/((a - 3)*(a - 2))")
P4[B4] = simplify("(a - 1)*(2*a - g - 2)/((a - 3)*(a - 2))")
P4[Q31] = simplify("-2*(a - g)*(a - g - 1)/((a - 3)*(a - 2))")
P4[Q211] = simplify("-(a - g)*(5*a - g - 5)/((a - 3)*(a - 2))")
P4[K4] = simplify("2*(a - g - 1)*(2*a - g - 2)/((a - 3)*(a - 2))")

def E(n):
	return simplify(binom(n, 4)*(q4*P4[S4] + q3*P4[B4] + q2*P4[Q31] + q1*P4[Q211]))

sum1 = simplify(q4**2*P4[S4] + q3**2*P4[B4] + q2**2*P4[Q31] + q1**2*P4[Q211])

sum2 = simplify((q4*P4[S4] + q3*P4[B4] + q2*P4[Q31] + q1*P4[Q211])**2)

S = dict()
for i in range(1, 5):
	S[i] = dict()
	for j in range(1, 5):
		S[i][j] = dict()
for line in open("variance_table.txt"):
	line = line.rstrip().split(' ', 3)
	i, j, k, formula = line[0], line[1], line[2], line[3]
	S[int(i)][int(j)][int(k)] = sympify(formula)

def sum3(n):
	sum = 0
	for i in range(1, 5):
		for j in range(1, 5):
			sumk = 0
			for k in range(5, 9):
				sumk += binom(n, k)*S[i][j][k]
			sum += q[i]*q[j]*sumk
	return simplify(sum)

def Var(n):
	bin4 = binom(n, 4)
	return simplify(bin4*sum1 - bin4**2*sum2 + sum3(n))

out = open("Var_2_64.txt", "w")
for n in range(2, 65):
	print "round", n, "of 64"
	print >>out, n, Var(n)

