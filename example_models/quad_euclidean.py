import sympy as sym
from pytsa.pytrans_setup import Translator

# Example installation file for single field quadratic inflation

nF = 1  # number of fields
nP = 1  # number of params

f = sym.symarray("f", nF)  # symbolic array for fields
p = sym.symarray("p", nP)  # symbolic array for params

V = (
    sym.Rational(1, 2) * f[0] ** 2 * p[0] ** 2
)  # Symbolic expression for the potential

# Translate model into c++ source code. Note that we pass the metric G=G
Translator(nF, nP, V)

# Compile module ! Should now be importable python module with prefix pyt
# , e.g. import pyt_dquad_2sphere as model
Translator.install("quad")
