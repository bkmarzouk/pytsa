import sympy as sym
from pytsa.pytrans_setup import Translator

# Example installation file for double quadratic field
# inflation with Euclidean metric

nF = 2  # Number of fields
nP = 2  # Number of params

f = sym.symarray("f", nF)  # Build symbolic array for fields
p = sym.symarray("p", nP)  # Build symbolic array for params

V = sum(
    [sym.Rational(1, 2) * f[i] ** 2 * p[i] ** 2 for i in range(2)]
)  # Construct symbolic expression for potential

# Translate model into c++ source code. Note that we pass the metric G=G
Translator(nF, nP, V)

# Compile module ! Should now be importable python module with prefix pyt
# , e.g. import pyt_dquad_2sphere as model
Translator.install("dquad_euclidean")
