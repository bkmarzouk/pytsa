import sympy as sym
from pytsa import pytrans_setup

# Example installation file for double quadratic field inflation with Euclidean metric

nF = 2  # Number of fields
nP = 2  # Number of params

f = sym.symarray('f', nF)  # Build symbolic array for fields
p = sym.symarray('p', nP)  # Build symbolic array for params

V = sum([sym.Rational(1, 2) * f[i] ** 2 * p[i] ** 2 for i in range(2)])  # Construct symbolic expression for potential

# Translate model into c++ source code
pytrans_setup.potential(V, nF, nP, simplify_fmet=True, simplify_pot=True, simplify_covd=True, silent=False)

# Compile module ! Should now be importable python module with prefix pyt_, e.g. import pyt_dquad_euclidean as model
pytrans_setup.compile_module('dquad_euclidean', False)
