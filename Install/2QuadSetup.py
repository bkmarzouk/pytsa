import os
import sympy as sym
import sys
tree = os.path.abspath(os.path.join(__file__, '../..'))
pytpath = os.path.join(tree, 'PyTransport', 'PyTransport')
sys.path.append(pytpath)
import PyTransSetup
nF = 2
nP = 2
f  = sym.symarray('f',nF)
p  = sym.symarray('p',nP)
s  = [sym.Rational(1,2) * f[i] ** 2 * p[i] ** 2 for i in range(2)]
V  = sum(s)

G = sym.Matrix.diag([1 for i in range(nF)])

PyTransSetup.tol(1E-8,1E-8)
PyTransSetup.potential(V,nF,nP,False, G, silent=False)
PyTransSetup.compileName('2Quad', True)


#PyTransSetup.potential(V,6,32,False,G,silent=False)
#PyTransSetup.compileName('agarwal_dmax_6pt0',True)
