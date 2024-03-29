import sympy as sym
from pytsa.pytrans_setup import Translator

# Example installation file for double quadratic field inflation
# with 2-sphere metric

nF = 2  # Number of fields
nP = 3  # Number of params (m^2 for each field and 3rd for the metric def.)

f = sym.symarray("f", nF)  # Build symbolic array for fields
p = sym.symarray("p", nP)  # Build symbolic array for params

V = sum(
    [sym.Rational(1, 2) * f[i] ** 2 * p[i] ** 2 for i in range(2)]
)  # Construct symbolic expression for potential

G = sym.Matrix(
    [[p[2] ** 2, 0], [0, p[2] ** 2 * sym.sin(f[0]) ** 2]]
)  # Construct symbolic expression for metric

# Translate model into c++ source code. Note that we pass the metric G=G
Translator(nF, nP, V, G=G)

# Compile module ! Should now be importable python module with prefix pyt,
# e.g. import pyt_dquad_2sphere as model
Translator.install("dquad_2sphere")
