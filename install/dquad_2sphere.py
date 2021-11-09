# from path_root import *  # Sets path to import PyTransport files
import sympy as sym
from pytransport import pytrans_setup

nF = 2
nP = 3
f = sym.symarray('f', nF)
p = sym.symarray('p', nP)
s = [sym.Rational(1, 2) * f[i] ** 2 * p[i] ** 2 for i in range(2)]
V = sum(s)
G = sym.Matrix([[p[2] ** 2, 0], [0, p[2] ** 2 * sym.sin(f[0]) ** 2]])

PyTransSetup.potential(V, nF, nP, G=G, simplify_fmet=True, simplify_pot=True, simplify_covd=True, silent=False)
PyTransSetup.compile_module('PyTransTEST', True)
