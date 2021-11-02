from path_root import *  # Sets path to import PyTransport files
import sympy as sym
from PyTransport import PyTransSetup

nF = 2
nP = 2
f = sym.symarray('f', nF)
p = sym.symarray('p', nP)
s = [sym.Rational(1, 2) * f[i] ** 2 * p[i] ** 2 for i in range(2)]
V = sum(s)

PyTransSetup.potential(V, nF, nP, simplify_fmet=True, simplify_pot=True, simplify_covd=True, silent=False)
PyTransSetup.compile_module('dquad_euclidean', True, use_j=True)

import distutils