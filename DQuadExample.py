import sympy as sym
from PyTransport import PyTransSetup

nF = 2
nP = 2
f = sym.symarray('f', nF)
p = sym.symarray('p', nP)
s = [sym.Rational(1, 2) * f[i] ** 2 * p[i] ** 2 for i in range(2)]
V = sum(s)

PyTransSetup.tol(1E-8, 1E-8)

PyTransSetup.potential(V, nF, nP, simpleGeometric=True, simplePotentials=True, silent=False)
PyTransSetup.compileName3('dquad', True)