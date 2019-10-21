import sympy as sym
import numpy as np
import scipy as sp
import itertools as it

from sympy.matrices import Matrix
from RiemannComponents import get_Riemann_Dictionary


class Curvature:
    """ Compute curvature quantities for a given metric.
        `metric' should be a sympy matrix
        'coords' should be sympy symbolic args for metric coordinates
        'params' is optional, but should be symbolic like fields
        """

    def __init__(self, metric, coords, params=None, precompute=False):

        print "-- Initializing curvature object"

        if metric == "canonical": self.metric = Matrix.diag(np.ones(len(coords))) # G_{IJ} = \delta_{IJ}
        else: self.metric = metric

        self.coords = list(coords)
        if params is not None: self.params  = list(params)
        else: self.params = None

        # We will assign precomputed quantities to the following arrays
        self.Csyms = np.zeros(
            (len(coords), len(coords), len(coords)), dtype=object)
        self.Rsyms = np.zeros(
            (len(coords), len(coords), len(coords), len(coords)), dtype=object
        )

        print self.metric
        self.dim = np.shape(self.metric)
        assert self.dim[0] == self.dim[1], "Metric dimensions must be equal"

        print "-- Inverting metric"

        self.metric_inverse = self.metric.inv()

        self.computed = False
        if precompute is True:
            self.precompute()
            self.computed = True

    def Christoffel(self, a, b, c, symbolic_args=True, simplify=False):
        """

        Computes connection component for given metric.

        a, b, c can be symbolic and contained within 'coords'
        or a, b, c can be indices which correspond to vars in 'coords'

        These are Christoffel symbols of the second kind, formatted as:
        Up index,     a
        Down indices, b, c
        i.e. \Gamma^a_{bc}

        """

        # Compute index placement from symbolic array of metric parameters

        if self.metric == 0: return 0

        if symbolic_args is True:
            # If symbolic items are given, retrieve indices from coords array
            a_idx = self.coords.index(a)
            b_idx = self.coords.index(b)
            c_idx = self.coords.index(c)
        else:
            # a, b, c are coordinate indices, retrieve symb. def from coords array
            a_idx = a
            b_idx = b
            c_idx = c
            a, b, c = (self.coords[i] for i in [a_idx, b_idx, c_idx])

        if self.computed is True:
            return self.Retrieve_Christoffel(a_idx, b_idx, c_idx)

        gamma = 0

        # Iterate over ther free index in the Levi Civita connection
        for free_idx in range(0, len(self.coords)):

            d = self.coords[free_idx]

            gamma += 0.5 * self.metric_inverse[a_idx, free_idx] * ( 0 # 0.5 * g^{ad}  * (
                + sym.diff(self.metric[free_idx, c_idx], b)           # grad_b g_{dc} +
                + sym.diff(self.metric[b_idx, free_idx], c)           # grad_c g_{bd} -
                - sym.diff(self.metric[b_idx, c_idx], d)              # grad_d g_{bc}   )
            )

        if simplify is True:
            return sym.simplify(gamma)

        return gamma

    def Riemann_inertial(self, d, c, a, b, symbolic_args=True, simplify=False):
        """ Computes riemann curvature tensor, all indices down: Rdcab, IN INERTIAL COORDINATES """

        if self.metric == 0: return 0

        if symbolic_args is True:
            d_idx = self.coords.index(d)
            c_idx = self.coords.index(c)
            a_idx = self.coords.index(a)
            b_idx = self.coords.index(b)
        else:
            d_idx = d
            c_idx = c
            a_idx = a
            b_idx = b
            d, c, a, b = (self.coords[i] for i in [d_idx, c_idx, a_idx, b_idx])

        Rdcab = 0.5 * ( 0
            + sym.diff(sym.diff(self.metric[b_idx, d_idx], c), a)
            - sym.diff(sym.diff(self.metric[b_idx, c_idx], d), a)
            - sym.diff(sym.diff(self.metric[a_idx, d_idx], c), b)
            + sym.diff(sym.diff(self.metric[a_idx, c_idx], d), b)
        )

        if simplify is True:
            return sym.simplify(Rdcab)

        return Rdcab

    def Riemann(self, d, c, a, b, symbolic_args=True, simplify=False, inertial=False):

        if inertial is True:
            return self.Riemann_inertial(d, c, a, b, symbolic_args, simplify)

        if symbolic_args is True:
            d_idx = self.coords.index(d)
            c_idx = self.coords.index(c)
            a_idx = self.coords.index(a)
            d_idx = self.coords.index(b)

        else:
            d_idx = d
            c_idx = c
            a_idx = a
            b_idx = b
            d, c, a, b = (self.coords[i] for i in [d_idx, c_idx, a_idx, b_idx])

        if self.computed is True:
            if simplify is True:
                return sym.simplify(self.Retrieve_Riemann(d_idx, c_idx, a_idx, b_idx))
            else:
                return self.Retrieve_Riemann(d_idx, c_idx, a_idx, b_idx)

        term1 = 0
        term2 = 0
        term3 = 0
        term4 = 0
        for f_idx in range(len(self.coords)):
            f = self.coords[f_idx]

            term1 +=  self.metric[d_idx, f_idx] * sym.diff(
                self.Christoffel(f, b, c), a
            )

            term2 += -self.metric[d_idx, f_idx] * sym.diff(
                self.Christoffel(f, a, c), b
            )

            for e_idx in range(len(self.coords)):
                e = self.coords[e_idx]

                term3 +=  self.metric[d_idx, f_idx] * self.Christoffel(f, a, e) * self.Christoffel(e, b, c)
                term4 += -self.metric[d_idx, f_idx] * self.Christoffel(f, b, e) * self.Christoffel(e, a, c)

        if simplify is True:
            return sym.simplify(term1 + term2 + term3 + term4)

        return term1 + term2 + term3 + term4

    def precompute(self):
        max_idx = len(self.coords)

        print "-- Precomputing Christoffel symbols"

        for a in range(max_idx):
            for b in range(max_idx):
                for c in range(b + 1):
                    self.Csyms[a, b, c] = self.Christoffel(a, b, c, symbolic_args=False)

        print "-- Precomputing Riemann tensors"

        print "   Finding fundamental components"

        grd = get_Riemann_Dictionary(max_idx)

        def str_to_list(ijlk): return (int(ijlk[0]), int(ijlk[1]), int(ijlk[2]), int(ijlk[3]))

        def exch(t):  return (t[2], t[3], t[0], t[1])

        def skew1(t): return (t[0], t[1], t[3], t[2])

        def skew2(t): return (t[1], t[0], t[2], t[3])

        # Compute all independent Riemann terms
        for abcd in grd:
            operation = grd[abcd]
            a, b, c, d = str_to_list(abcd)

            if operation == "CALC":
                self.Rsyms[a, b, c, d] = self.Riemann(a, b, c, d, symbolic_args=False, simplify=False)

        print "   Assigning zero terms and imposing symmetry relations"

        # Assign all zero terms
        for abcd in grd:
            operation = grd[abcd]
            a, b, c, d = str_to_list(abcd)

            if operation == 0: self.Rsyms[a, b, c, d] = 0

        # Assign values to remaining terms via symmetries
        for abcd in grd:
            operation  = grd[abcd]
            a, b, c, d = str_to_list(abcd)
            new        = (a, b, c, d)
            minus_exp  = 0

            try: # Cyclic terms have keys like 'abcd+efgh'; hence, we should be able to convert 0th element to int
                int(operation[0])
                Rabcd1, Rabcd2 = operation.split('+')
                idx1 = str_to_list(Rabcd1); idx2 = str_to_list(Rabcd2)
                self.Rsyms[a, b, c, d] = self.Rsyms[
                                             idx1[0], idx1[1], idx1[2], idx1[3]] + self.Rsyms[
                                             idx2[0], idx2[1], idx2[2], idx2[3] ]

            except: # Error thrown implies symmetry impositions
                if operation != 0 and operation != "CALC":
                    split = operation.split('~')
                    execs = []

                    for op in split:
                        if   op == 'x'  :
                            execs.append(exch)
                        elif op == 's1' :
                            minus_exp+=1
                            execs.append(skew1)
                        elif op == 's2' :
                            minus_exp+=1
                            execs.append(skew2)
                        else            : pass

                    for e in execs: new = e(new)

                    self.Rsyms[a, b, c, d] = (-1)**minus_exp * self.Rsyms[new[0], new[1], new[2], new[3]]

        print "   Done"

    def Retrieve_Christoffel(self, a, b, c):
        if c > b:
            return self.Csyms[a, c, b]
        else:
            return self.Csyms[a, b, c]

    def Retrieve_Riemann(self, d, c, a, b):
        return self.Rsyms[d, c, a, b]

    def Riemann_up(self, d, c, a, b, up_idx=0, symbolic_args=True, simplify=False):
        n = len(self.coords)
        r = 0

        if symbolic_args is True:
            d_idx = self.coords.index(d)
            c_idx = self.coords.index(c)
            a_idx = self.coords.index(a)
            d_idx = self.coords.index(b)

        else:
            d_idx = d
            c_idx = c
            a_idx = a
            b_idx = b
            d, c, a, b = (self.coords[i] for i in [d_idx, c_idx, a_idx, b_idx])

        for free_idx in range(n):
            coords = [d_idx, c_idx, a_idx, b_idx]
            coords[up_idx] = free_idx
            r+=self.metric_inverse[d_idx, free_idx]*self.Riemann(
                coords[0], coords[1], coords[2], coords[3], symbolic_args=symbolic_args, simplify=simplify
            )

        if simplify is True: return sym.simplify(r)
        return r

#
# x, y = sym.symbols('x, y')
# metric = Matrix([
#    [y/x, x/y],
#    [x/y, 0]
# ])
#
# a=Curvature(metric, (x, y), None, precompute=True)
#
# print a.metric_inverse
#
# for i in range(2):
#    for j in range(2):
#        for k in range(2):
#            print [i, j, k], a.Christoffel(i, j, k, symbolic_args=False, simplify=True)
#            for l in range(2):
#                print [i, j, k, l], a.Riemann(i, j, k, l, symbolic_args=False, simplify=True)

# t, r, theta, phi, M = sym.symbols('t, r, \\theta, \phi, M') # define some symbolic variables
# metric = Matrix.diag(-(1-2*M/r), 1/(1-2*M/r), r**2, r**2*sym.sin(theta)**2) # define a matrix of a metric tensor components
# a = Curvature("canonical", (t, r, theta, phi), (M,), precompute=True)
#
# print sym.simplify(a.Riemann_up(0, 1, 0, 1, symbolic_args=False))
# print sym.simplify(a.Riemann_up(0, 1, 1, 0, symbolic_args=False))
# print sym.simplify(a.Riemann_up(0, 3, 0, 3, symbolic_args=False))
# print sym.simplify(a.Riemann_up(1, 3, 1, 3, symbolic_args=False))

# for i in range(4):
#     for j in range(4):
#         for k in range(4):
#             c= a.Christoffel(i, j, k, symbolic_args=False)
#             if c != 0:
#                 print c
#             for l in range(4):
#                 r = a.Riemann(i, j, k, l, symbolic_args=False)
#                 if r != 0:
#                     print r

# for i in range(4):
#     for j in range(4):
#         for k in range(4):
#             for l in range(4):
#                 p = sym.simplify(a.Riemann(i, j, l, k, symbolic_args=False))
#                 q = sym.simplify(a.Riemann_up(i, j, l, k))
#
#                 print p, "->", q


# f = sym.symbols('x t_1 t_2 p_1 p_2 psi')
# G = Matrix( [ [1., 0, 0, 0, 0, 0],
#               [0, (f[0]**2)*1./6., 0, 0, 0, 0],
#               [0, 0, (f[0]**2)*1./6., 0, 0, 0],
#               [0, 0, 0, (f[0]**2)*(1./9.)*(sym.cos(f[1])**2.)+(1./6.)*(sym.sin(f[1])**2.), (f[0]**2)*(1./9.)*sym.cos(f[1])*sym.cos(f[2]), (f[0]**2)*(1./9.)*(sym.cos(f[1]))],
#               [0, 0, 0, (f[0]**2)*(1./9.)*sym.cos(f[1])*sym.cos(f[2]), (f[0]**2)*(1./9.)*(sym.cos(f[2])**2.)+(1./6.)*(sym.sin(f[2])**2.), (f[0]**2)*(1./9.)*(sym.cos(f[2]))],
#               [0, 0, 0, (f[0]**2)*(1./9.)*(sym.cos(f[1])), (f[0]**2)*(1./9.)*(sym.cos(f[2])), (f[0]**2)*1./9.]
#             ] ) # Conifold metric in sympy notation; metric components g ~ {r, theta1, theta2, phi1, phi2, psi}
#
#
# b = Curvature(G, f, None, precompute=True)
