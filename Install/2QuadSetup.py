import os
import sympy as sym
import sys
tree = os.path.abspath(os.path.join(__file__, '../..'))
pytpath = os.path.join(tree, 'PyTransport', 'PyTransport')
sys.path.append(pytpath)
import PyTransSetup

nF = 2
nP = 2
f  = sym.symarray('f',nF)
p  = sym.symarray('p',nP)
s  = [sym.Rational(1,2) * f[i] ** 2 * p[i] ** 2 for i in range(2)]
V  = sum(s)

PyTransSetup.tol(1E-8,1E-8)

PyTransSetup.potential(V, nF, nP, simpleGeometric=True, simplePotentials=True, silent=False)
PyTransSetup.compileName('2Quad', True)