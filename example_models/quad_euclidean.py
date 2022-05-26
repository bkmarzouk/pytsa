from path_root import *  # Sets path to import pytsa files
import sympy as sym
from pytsa import pytrans_setup

# Example installation file for single field quadratic inflation

nF = 1  # number of fields
nP = 1  # number of params

f = sym.symarray('f', nF)  # symbolic array for fields
p = sym.symarray('p', nP)  # symbolic array for params

V = sum([sym.Rational(1, 2) * f[i] ** 2 * p[i] ** 2 for i in range(nF)])  # Symbolic expression for the potential

# Run translator for building c++ source code
pytrans_setup.potential(V, nF, nP, simplify_fmet=True, simplify_pot=True, simplify_covd=True, silent=False)

# Compile model ! Should now be importable with prefix pyt_, e.g. import pyt_quad as model
pytrans_setup.compile_module('quad', False)
