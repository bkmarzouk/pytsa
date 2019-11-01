import sympy as sym

# Number of fields
nF = 3
rnf = range(nF)

# Metric
G = sym.symarray("G", (nF, nF)) # Contravariant
g = sym.symarray("g", (nF, nF)) # Covariant

# Hubble rate
H = sym.symbols("H")

# Potential derivatives
dV = sym.symarray("dV", nF)
ddV = sym.symarray("ddV", (nF, nF))

# Fields & field time derivatives, d = covariant, u = contravariant index
phi = sym.symarray("f", nF)
dphi = sym.symarray("df", nF)
ddphi = sym.symarray("ddf", nF)

# Define kinetic energy of fields
KE = sym.Rational(1, 2) * sum([g[a, b]*dphi[a]*dphi[b] for a in rnf for b in rnf])

# Curvature quantities
Riemann = sym.symarray("R", (nF, nF, nF, nF)) # Riemann tensor with first index up
Christoffel = sym.symarray("C", (nF, nF, nF)) # Christoffel symbols, 1st index ip, 2nd & 3rd down

CTerm0 = sym.zeros(nF, nF)
RTerm0 = sym.zeros(nF, nF)
KTerm1 = sym.zeros(nF, nF)
KTerm2 = sym.zeros(nF, nF)
KTerm3 = sym.zeros(nF, nF)

for i in rnf:
    for j in rnf:
        for x in rnf:
            CTerm0[i, j] += ddV[x, j]
            KTerm1[i, j] += -(3 - KE / H / H) * dphi[i] * g[j, x] * dphi[x]
            
            for k in rnf:
                CTerm0[i, j] += -Christoffel[k, x, j]*dV[k]
                RTerm0[i, j] += -Riemann[i, x, k, j] * dphi[x]*dphi[k]
                
                for a in rnf:
                    KTerm1[i, j] += -(ddphi[i] + dphi[k]*Christoffel[i, a, k]*dphi[a]) * g[j, x]*dphi[x] / H
                    
                    for b in rnf:
                        KTerm2[i, j] += dphi[i]*()

            CTerm0 *= G[x, i]