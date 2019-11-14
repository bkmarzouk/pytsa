import sympy as sym
import numpy as np
import scipy as sp
from sympy.matrices import Matrix

from idp import OperationMatrix


class Curvature:
    """
    Compute curvature quantities for a given metric
    `metric' should be a sympy matrix
    'coords' should be sympy symbolic args for metric coordinates that define the metric
    'params' is optional, but should be symbolic like 'coords'
    'precompute' computes all curvature quantities upon initialization
    """
    
    
    def __init__(self, metric, coords, params=None, precompute=False):
        
        print "-- Initializing curvature object: gravtools_pyt"
        
        # If canonical, assume Kronecker-delta, i.e. identity matrix representation for the metric
        if metric == "canonical":
            self.metric = Matrix.diag(np.ones(len(coords)))  # G_{IJ} = \delta_{IJ}
        else:
            self.metric = metric
        
        # Distinguish coordinates from other symbols which may be in the metric, e.g. mass parameters
        self.coords = list(coords)
        
        # Keep track of parameters, if we wish
        if params is not None:
            self.params = list(params)
        else:
            self.params = None
        
        # Initialize zeros array to store Christoffel symbols
        self.Csyms = np.zeros(
            (len(coords), len(coords), len(coords)), dtype=object
        )
        
        # Initialize zeros array to store Riemann tensors
        self.Rsyms = np.zeros(
            (len(coords), len(coords), len(coords), len(coords)), dtype=object
        )
        
        # Perform consistency check on in put metric
        self.dim = np.shape(self.metric)
        assert self.dim[0] == self.dim[1], "Metric dimensions must be NxN: {}x{}".format(self.dim[0], self.dim[1])
        
        print "-- Inverting metric"
        
        # Compute inverse metric
        self.metric_inverse = self.metric.inv()
        
        print "-- Simplifying inverse metric"
        
        ncr = range(len(self.coords))
        for a in ncr:
            for b in ncr:
                self.metric_inverse[a, b] = sym.simplify(self.metric_inverse[a, b])
        
        # Check that the metric and inverse metric are both symmetric
        assert self.metric.is_symmetric(), "Metric is not symmetric! {}".format(self.metric)
        assert self.metric_inverse.is_symmetric(), "Inverse metric is not symmetric! {}".format(self.metric_inverse)
        
        # Consistency check: Product of G * G^-1 = Identity
        identity_matrix = self.metric * self.metric_inverse
        identity_matrix = sym.simplify(identity_matrix)
        assert identity_matrix == sym.Matrix.diag([1 for c in self.coords]), "G * G^-1 != I: {}".format(identity_matrix)
        
        # Assign precomputed status as False by default
        self.computed = False
        
        # If precompute is True, compute all curvature quantities and change status
        if precompute is True:
            self.precompute()
            self.computed = True
        
        print "-- Curvature object constructed"
    
    
    def Christoffel(self, a, b, c, symbolic_args=True, simplify=False):
        """
        Computes Christoffel symbols of the second kind, \Gamma^a_{bc}
        a, b, c can correspond to index elements of the initialized coordinates or
        a, b, c can be symbolic, the same as those prescribed in initialization
        """
        
        if symbolic_args is True:
            # If symbolic items are given, retrieve indices from coords array
            a_idx = self.coords.index(a)
            b_idx = self.coords.index(b)
            c_idx = self.coords.index(c)
        else:
            # a, b, c are coordinate indices, retrieve symbolic definition from coords array
            a_idx = a
            b_idx = b
            c_idx = c
            a, b, c = (self.coords[i] for i in [a_idx, b_idx, c_idx])
        
        # If precomputed, get element from Christoffel matrix
        if self.computed is True:
            gamma = self.Csyms[a_idx, b_idx, c_idx]
        
        # Otherwise compute directly
        else:
            gamma = 0
            
            # Iterate over the repeated index, d
            for d_idx in range(0, len(self.coords)):
                d = self.coords[d_idx]
                
                gamma += sym.Rational(1, 2) * self.metric_inverse[a_idx, d_idx] * sum([
                    sym.diff(self.metric[d_idx, c_idx], b),  # partial_b g_{dc} +
                    sym.diff(self.metric[b_idx, d_idx], c),  # partial_c g_{bd} -
                    -sym.diff(self.metric[b_idx, c_idx], d)  # partial_d g_{bc}
                ])
        
        # Return Christoffel symbol with simplification as desired
        if simplify is True: gamma = sym.simplify(gamma)
        
        return gamma
    
    
    def RiemannUp(self, f, c, a, b, symbolic_args=True, simplify=False):
        """ We define the Riemann Curvature tensor, R^{f}_{cab} """
        
        if symbolic_args is False:
            (f, c, a, b) = (self.coords[i] for i in [f, c, a, b])
        
        term1 = sym.diff(self.Christoffel(f, b, c), a)
        term2 = -sym.diff(self.Christoffel(f, a, c), b)
        term3 = sum([self.Christoffel(f, a, d) * self.Christoffel(d, b, c) -
                     self.Christoffel(f, b, d) * self.Christoffel(d, a, c) for d in self.coords])
        
        Rfcab = term1 + term2 + term3
        
        if simplify is True: Rfcab = sym.simplify(Rfcab)
        
        return Rfcab
    
    
    def Riemann(self, d, c, a, b, symbolic_args=True, simplify=False):
        """ We define the Riemann Curvature tensor in full contravariant form, R_{dcab} """
        
        if symbolic_args is True:
            (d_idx, c_idx, a_idx, b_idx) = (self.coords.index(s) for s in [d, c, a, b])
        else:
            (d_idx, c_idx, a_idx, b_idx) = (d, c, a, b)
        
        if self.computed is True:
            Rdcab = self.Rsyms[d_idx, c_idx, a_idx, b_idx]
        
        else:
            Rdcab = sum([self.metric[d_idx, f_idx] * self.RiemannUp(
                f_idx, c_idx, a_idx, b_idx, symbolic_args=False) for f_idx in range(len(self.coords))])
        
        if simplify is True: Rdcab = sym.simplify(Rdcab)
        
        return Rdcab
    
    
    def precompute(self):
        """ Precomputes all components of the Riemann tensor and Christoffel symbols, implementing symmetry relations """
        
        max_idx = len(self.coords)
        
        print "-- Precomputing Christoffel symbols"
        
        for a in range(max_idx):
            for b in range(max_idx):
                for c in range(max_idx):
                    self.Csyms[a, b, c] = self.Christoffel(a, b, c, symbolic_args=False, simplify=True)
        
        print "-- Precomputing Riemann tensors"
        
        # Get operation matrix, i.e. identify fundamental components and associated symmetries
        OM = OperationMatrix(max_idx)
        
        assert np.shape(OM) == np.shape(self.Rsyms), "Operation matrix incompatible with Riemann matrix"
        
        # Get location of independent elements to compute & symmetrically related components
        calc_locs = np.where(OM == "calc")
        symm_locs = np.where(np.logical_and(OM != "calc", OM != 0))
        
        # Transpose into rows of i, j, k, l matrix coordinates / tensor indices
        calc_ijkl = np.asarray(calc_locs).T
        symm_ijkl = np.asarray(symm_locs).T
        
        # Count fundamental component computations dynamically
        c = 1
        m = len(calc_ijkl)
        
        # Compute independent components
        for ijkl in calc_ijkl:
            print "computing {}/{}".format(c, m)
            
            i, j, k, l = ijkl
            
            assert self.Rsyms[i, j, k, l] == 0, "Riemann component already has assigned value! {}".format(
                self.Rsyms[i, j, k, l]
            )
            
            Rijkl = self.Riemann(i, j, k, l, symbolic_args=False)
            
            self.Rsyms[i, j, k, l] = sym.simplify(Rijkl)
            
            c += 1
        
        # Assign symmetrically related components
        for ijkl in symm_ijkl:
            
            i, j, k, l = ijkl
            
            assert self.Rsyms[i, j, k, l] == 0, "Riemann component already has assigned value! {}".format(
                self.Rsyms[i, j, k, l]
            )
            
            # Unpack the string that points to the fundamental Matrix element(s) and multiplicative sign
            str_pointer = OM[i, j, k, l]
            
            if len(str_pointer) == 5:
                sign, p, q, r, s = str_pointer
                
                Rpqrs = self.Rsyms[int(p), int(q), int(r), int(s)]
                
                if sign == "+":
                    self.Rsyms[i, j, k, l] += Rpqrs
                elif sign == "-":
                    self.Rsyms[i, j, k, l] += -Rpqrs
                else:
                    raise ValueError, "Unrecognized sign value: {}".format(sign)
            
            elif len(str_pointer) == 10:
                sign1, p, q, r, s, sign2, t, u, v, w = str_pointer
                
                Rpqrs = self.Rsyms[int(p), int(q), int(r), int(s)]
                Rtuvw = self.Rsyms[int(t), int(u), int(v), int(w)]
                
                if sign1 == "+":
                    self.Rsyms[i, j, k, l] += Rpqrs
                elif sign1 == "-":
                    self.Rsyms[i, j, k, l] += -Rpqrs
                else:
                    raise ValueError, "Unrecognized sign value: {}".format(sign1)
                
                if sign2 == "+":
                    self.Rsyms[i, j, k, l] += Rtuvw
                elif sign2 == "-":
                    self.Rsyms[i, j, k, l] += -Rtuvw
                else:
                    raise ValueError, "Unrecognized sign value: {}".format(sign2)
            
            else:
                raise ValueError, "Unrecognized symmetry relation: {}".format(str_pointer)

# f = sym.symarray("f", 6)
# p = sym.symarray("p", 32)
#
# from sympy import Rational as rat
#
# G = sym.Matrix(
# 	[
# 		[f[0]**-2, 0, 0, 0, 0, 0],
# 		[0, rat(1, 6), 0, 0, 0, 0],
# 		[0, 0, rat(1, 6), 0, 0, 0],
# 		[0, 0, 0, rat(1, 9)*sym.cos(f[1])**2 + rat(1, 6)*sym.sin(f[1])**2, sym.cos(f[1])*sym.cos(f[2])*rat(1, 9), sym.cos(f[1])*rat(1, 9)],
# 		[0, 0, 0, sym.cos(f[1])*sym.cos(f[2])*rat(1, 9), rat(1, 9)*sym.cos(f[2])**2 + rat(1, 6)*sym.sin(f[2])**2, sym.cos(f[2])*rat(1, 9)],
# 		[0, 0, 0, sym.cos(f[1])*rat(1, 9), sym.cos(f[2])*rat(1, 9), rat(1, 9)]
# 	]
# )
# G *= p[3]**2*f[0]**2
#
# db_obj = Curvature(G, f, p, precompute=True)
#
# import pickle as pk
#
# f = open("db.curv", "wb")
#
# with f: pk.dump(db_obj, f)
