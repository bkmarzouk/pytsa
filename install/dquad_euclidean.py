from path_root import *  # Sets path to import pytransport files
import sympy as sym
from pytransport import pytrans_setup

nF = 2
nP = 2
f = sym.symarray('f', nF)
p = sym.symarray('p', nP)
s = [sym.Rational(1, 2) * f[i] ** 2 * p[i] ** 2 for i in range(2)]
V = sum(s)

pytrans_setup.potential(V, nF, nP, simplify_fmet=True, simplify_pot=True, simplify_covd=True, silent=False)
pytrans_setup.compile_module('dquad_euclidean', False)