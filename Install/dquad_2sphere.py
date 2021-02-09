from path_root import *  # Sets path to import PyTransport files
import sympy as sym
from PyTransport import PyTransSetup

nF = 2
nP = 3
f = sym.symarray('f', nF)
p = sym.symarray('p', nP)
s = [sym.Rational(1, 2) * f[i] ** 2 * p[i] ** 2 for i in range(2)]
V = sum(s)
G = sym.Matrix([[p[2] ** 2, 0], [0, p[2] ** 2 * sym.sin(f[0]) ** 2]])

PyTransSetup.tol(1E-8, 1E-8)
PyTransSetup.potential(V, nF, nP, G=G, simple_fmet=True, simple_potential=True, silent=False)
PyTransSetup.compileName('dquad_2sphere', True)
