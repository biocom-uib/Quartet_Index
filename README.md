# Quartet Index

This is the companion software for the paper *A Balance Index for Phylogenetic Trees based on Quartets* by Tomàs Martínez, Arnau Mir, Francesc Rosselló, and Gabriel Valiente.

## Getting Started

You need to install the ETE Toolkit, a Python library to work with phylogenetic trees. See http://etetoolkit.org/ for more details.

```
sudo apt install python-numpy python-qt4 python-lxml python-six
```

```
pip install ete3
```

You also need to install SymPy, a Python library for symbolic mathematics. See http://www.sympy.org/ for more details.

```
pip install sympy
```

## Usage Instructions

The Quartet Index companion software consists of a series of Python 2 scripts.

### First script

```
python qi_1_probabilities_shapes.py
```

This script will generate files **P2.txt** through **P8.txt**. These files have one line for each tree shape of the size in the file name, 2 through 8, with the tree shape and the symbolic formula (as a function of &alpha; and &gamma;) for the variance of that tree shape under the &alpha;-&gamma; model. For instance, there are 5 tree shapes with 4 leaves, and file **P4.txt** contains the following 5 lines:

```
((*,*),(*,*)); (a - 1)*(2*a - g - 2)/((a - 3)*(a - 2))
(*,(*,(*,*))); 2*(a - g - 1)*(2*a - g - 2)/((a - 3)*(a - 2))
(*,(*,*,*)); -2*(a - g)*(a - g - 1)/((a - 3)*(a - 2))
(*,*,(*,*)); -(a - g)*(5*a - g - 5)/((a - 3)*(a - 2))
(*,*,*,*); (a - g)*(2*a - g)/((a - 3)*(a - 2))
```

Grab some coffee. It can take a few hours to generate these files. Nevertheless, the files are included in this repository, for your convenience.
<!--
user	386m14.535s
-->

### Second script

```
python qi_2_variance_table.py
```

This script will will generate file **variance_table.txt**. This file has one line for each value of i, j between 1 and 4 and for each value of k between 5 and 8, with the values of i, j, k, and the symbolic formula (as a function of &alpha; and &gamma;) for the last sum in the formula for the variance of QI<sub>n</sub> under the &alpha;-&gamma; model (Proposition 3 in the aforementioned paper). For instance, the first 4 lines of file **variance_table.txt** are:

```
1 1 5 q1**2*(-78*a**3 + 66*a**2*g + 306*a**2 + 18*a*g**2 - 324*a*g - 228*a - 6*g**3 + 18*g**2 + 228*g)/(a**3 - 9*a**2 + 26*a - 24)
1 1 6 q1**2*(2634*a**4 - 5592*a**3*g - 3876*a**3 + 3780*a**2*g**2 + 9396*a**2*g - 1914*a**2 - 888*a*g**3 - 6444*a*g**2 - 252*a*g + 3156*a + 66*g**4 + 924*g**3 + 2166*g**2 - 3156*g)/(a**4 - 14*a**3 + 71*a**2 - 154*a + 120)
1 1 7 q1**2*(3500*a**5 - 8400*a**4*g - 45840*a**4 + 6440*a**3*g**2 + 101320*a**3*g + 71100*a**3 - 1680*a**2*g**3 - 70600*a**2*g**2 - 159440*a**2*g - 8600*a**2 + 140*a*g**4 + 16280*a*g**3 + 102340*a*g**2 + 43840*a*g - 20160*a - 1160*g**4 - 14000*g**3 - 35240*g**2 + 20160*g)/(a**5 - 20*a**4 + 155*a**3 - 580*a**2 + 1044*a - 720)
1 1 8 q1**2*(1750*a**6 - 4200*a**5*g - 33250*a**5 + 3220*a**4*g**2 + 79100*a**4*g + 226770*a**4 - 840*a**3*g**3 - 59640*a**3*g**2 - 509760*a**3*g - 332710*a**3 + 70*a**2*g**4 + 14980*a**2*g**3 + 357610*a**2*g**2 + 734940*a**2*g + 85600*a**2 - 1190*a*g**4 - 80060*a*g**3 - 463830*a*g**2 - 241760*a*g + 51840*a + 5440*g**4 + 61600*g**3 + 156160*g**2 - 51840*g)/(a**6 - 27*a**5 + 295*a**4 - 1665*a**3 + 5104*a**2 - 8028*a + 5040)
```

Grab some more coffee. It can take a few more hours to generate this file. Nevertheless, the file is included in this repository, for your convenience.
<!--
user	193m58.921s
-->

### Third script

```
python qi_3_variance_alpha_gamma.py
```

This script will generate file **Var_2_64.txt**. This file has one line for each value of n between 2 and 64, with the value of n and the symbolic formula for the variance of QI<sub>n</sub> under the &alpha;-&gamma; model as a function of &alpha; and &gamma;. For instance, the first 4 lines of file **Var_2_64.txt** are:

```
2 0
3 0
4 ((a - 3)*(a - 2)*(q1**2*(a - g)*(-5*a + g + 5) + 2*q2**2*(a - g)*(-a + g + 1) - q3**2*(a - 1)*(-2*a + g + 2) + q4**2*(a - g)*(2*a - g)) - (q1*(a - g)*(-5*a + g + 5) + 2*q2*(a - g)*(-a + g + 1) - q3*(a - 1)*(-2*a + g + 2) + q4*(a - g)*(2*a - g))**2)/((a - 3)**2*(a - 2)**2)
5 (5*(a - 4)*(a - 3)*(a - 2)*(a**3 - 9*a**2 + 26*a - 24)*(q1**2*(a - g)*(-5*a + g + 5) + 2*q2**2*(a - g)*(-a + g + 1) - q3**2*(a - 1)*(-2*a + g + 2) + q4**2*(a - g)*(2*a - g)) - 25*(a - 4)*(a**3 - 9*a**2 + 26*a - 24)*(q1*(a - g)*(-5*a + g + 5) + 2*q2*(a - g)*(-a + g + 1) - q3*(a - 1)*(-2*a + g + 2) + q4*(a - g)*(2*a - g))**2 + 2*(a - 3)*(a - 2)*(-6*q3**4*(a - 2)*(a - 1)*(-3*a + 2*g + 3)*(a**3 - 9*a**2 + 26*a - 24) + (a - 4)*(a - 3)*(a - 2)*(3*q1**4*(-13*a**3 + 11*a**2*g + 51*a**2 + 3*a*g**2 - 54*a*g - 38*a - g**3 + 3*g**2 + 38*g) + q2**4*(11*a**3 - 38*a**2*g + 9*a**2 + 39*a*g**2 - 21*a*g - 20*a - 12*g**3 + 12*g**2 + 20*g)) - (a - g)*(a**3 - 9*a**2 + 26*a - 24)*(2*q1**2*q2**2*(3*(a - g)*(-7*a + 3*g + 7) + (-5*a + g + 5)*(-2*a + 3*g + 2)) - 12*q1**2*q3**2*(a - 1)*(-4*a + g + 4) + 6*q1**2*q4**2*(2*a - g)*(-9*a + g + 9) - 12*q2**2*q3**2*(a - 1)*(-3*a + 2*g + 3) + 4*q2**2*q4**2*(2*a - g)*(-2*a + 3*g + 2) + 3*q4**4*(2*a - g)*(7*a - 3*g + 3))))/((a - 4)*(a - 3)**2*(a - 2)**2*(a**3 - 9*a**2 + 26*a - 24))
```

Grab still some more coffee. It can take half an hour to generate this file. Nevertheless, the file is included in this repository, for your convenience.
<!--
user	33m59.299s
-->

### Fourth script

```
python qi_4_variance_alpha_gamma_n.py
```

This script will generate file **Var_n.txt**. This file has only one line, with the symbolic formula for the variance of QI<sub>n</sub> under the &alpha;-&gamma; model as a function of &alpha;, &gamma; and n.

It shall only take a few minutes to generate this file. Nevertheless, the file is included in this repository, for your convenience.
<!--
user	7m8.109s
-->

### Fifth script

```
python qi_5_variance.py
```

This script will take values for &alpha;, &gamma; and n, and print the value of the variance of QI<sub>n</sub> under the &alpha;-&gamma; model. For instance, for XXX, XXX and XXX, it will print:

```
xxx
```

## Authors

* **Tomàs Martínez**
* **Arnau Mir** - [Email](mailto:arnau.mir@uib.es)
* **Francesc Rosselló** - [Email](mailto:cesc.rossello@uib.es)
* **Gabriel Valiente** - [Email](mailto:valiente@cs.upc.edu)

