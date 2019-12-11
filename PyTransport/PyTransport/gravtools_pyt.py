import sympy as sym
import numpy as np
import scipy as sp
from sympy.matrices import Matrix

from idp import OperationMatrix


class Curvature:
    
    
    def __init__(self, metric, coords, potential, params=None):
        
        print "-- Initializing curvature object: gravtools_pyt"
        
        # If canonical, assume Kronecker-delta, i.e. identity matrix representation for the metric
        if metric == "canonical" or metric is None or metric == 0:
            self.metric = Matrix.diag(np.ones(len(coords)))  # G_{IJ} = \delta_{IJ}
            self.canonical = True
        else:
            self.metric = metric
            self.canonical = False
        
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

        # Initialize zeros array to store covariant derivative(s) of Riemann tensor
        self.covDRsyms = np.zeros(
            (len(coords), len(coords), len(coords), len(coords), len(coords)), dtype=object
        )
        
        # Initialize similar arrays to store covariant derivatives of the potential
        self.V = potential
        self.dV = np.zeros((len(coords)), dtype=object)
        self.ddV = np.zeros((len(coords), len(coords)), dtype=object)
        self.dddV = np.zeros((len(coords), len(coords), len(coords)), dtype=object)
        
        # Initialize arrays to store partial derivatives of the potential
        self.pdV = np.zeros((len(coords)), dtype=object)
        self.pdpdV = np.zeros((len(coords), len(coords)), dtype=object)
        self.pdpdpdV = np.zeros((len(coords), len(coords), len(coords)), dtype=object)
        
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
                
        print "-- Computing remaining metric index combinations"
        self.metric_updown = np.zeros((len(self.coords), len(self.coords)), dtype=object)
        for a in ncr:
            for b in ncr:
                self.metric_updown[a, b] = sym.simplify(sum([self.metric_inverse[a, c]*self.metric[c, b] for c in ncr]))
        
        # Check that the metric and inverse metric are both symmetric
        assert self.metric.is_symmetric(), "Metric is not symmetric! {}".format(self.metric)
        assert self.metric_inverse.is_symmetric(), "Inverse metric is not symmetric! {}".format(self.metric_inverse)
        
        # Consistency check: Product of G * G^-1 = Identity
        identity_matrix = self.metric * self.metric_inverse
        identity_matrix = sym.simplify(identity_matrix)
        assert identity_matrix == sym.Matrix.diag([1 for c in self.coords]), "G * G^-1 != I: {}".format(identity_matrix)
        
        print "-- Precomputing symbolic expressions"
        
        self.precompute()
        
        print "-- Curvature object constructed"
    
    
    def Christoffel(self, a, b, c, symbolic_args=True):
        """
        Computes Christoffel symbols of the second kind, \Gamma^a_{bc}
        a, b, c can correspond to index elements of the initialized coordinates or
        a, b, c can be symbolic, the same as those prescribed in initialization
        """
        
        if self.canonical is True: return 0
        
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
        

        gamma = 0
        
        # Iterate over the repeated index, d
        for d_idx in range(0, len(self.coords)):
            d = self.coords[d_idx]
            
            gamma += sym.Rational(1, 2) * self.metric_inverse[a_idx, d_idx] * sum([
                sym.diff(self.metric[d_idx, c_idx], b),  # partial_b g_{dc} +
                sym.diff(self.metric[b_idx, d_idx], c),  # partial_c g_{bd} -
                -sym.diff(self.metric[b_idx, c_idx], d)  # partial_d g_{bc}
            ])
        
        # Return Christoffel symbol with simplification
        gamma = sym.simplify(gamma)
        
        return gamma
        
    
    def RiemannUp(self, f, c, a, b, symbolic_args=True):
        """ We define the Riemann Curvature tensor, R^{f}_{cab} """

        if self.canonical is True: return 0
        
        if symbolic_args is False:
            f_idx, c_idx, a_idx, b_idx = f, c, a, b
            (f, c, a, b) = (self.coords[i] for i in [f, c, a, b])
        
        term1 =  sym.diff(self.Csyms[f_idx, b_idx, c_idx], a)
        term2 = -sym.diff(self.Csyms[f_idx, a_idx, c_idx], b)
        term3 = sum([self.Csyms[f_idx, a_idx, d_idx] * self.Csyms[d_idx, b_idx, c_idx] -
                     self.Csyms[f_idx, b_idx, d_idx] * self.Csyms[d_idx, a_idx, c_idx]
                     for d_idx in range(len(self.coords))])
        
        Rfcab = term1 + term2 + term3
        
        Rfcab = sym.simplify(Rfcab)
        
        return Rfcab
    
    
    def Riemann(self, d, c, a, b, symbolic_args=True):
        """ We define the Riemann Curvature tensor in full contravariant form, R_{dcab} """
        
        if self.canonical is True:
            return 0
        
        if symbolic_args is True:
            (d_idx, c_idx, a_idx, b_idx) = (self.coords.index(s) for s in [d, c, a, b])
        else:
            (d_idx, c_idx, a_idx, b_idx) = (d, c, a, b)

        Rdcab = sum([self.metric[d_idx, f_idx] * self.RiemannUp(
            f_idx, c_idx, a_idx, b_idx, symbolic_args=False) for f_idx in range(len(self.coords))])
        
        Rdcab = sym.simplify(Rdcab)
        
        return Rdcab


    def CovDRiemann(self, e, d, c, a, b, symbolic_args=True):
    
        if self.canonical is True: return 0
    
        if symbolic_args is True:
            (e_idx, d_idx, c_idx, a_idx, b_idx) = (self.coords.index(s) for s in [e, d, c, a, b])
        else:
            e_idx, d_idx, c_idx, a_idx, b_idx = e, d, c, a, b
            (e, d, c, a, b) = (self.coords[idx] for idx in [e_idx, d_idx, c_idx, a_idx, b_idx])
    
        partialD_term = sym.diff(self.Rsyms[d_idx, c_idx, a_idx, b_idx], e)
    
        connexion_terms = -sum([
            self.Csyms[x_idx, d_idx, e_idx] * self.Rsyms[x_idx, c_idx, a_idx, b_idx] +
            self.Csyms[x_idx, c_idx, e_idx] * self.Rsyms[d_idx, x_idx, a_idx, b_idx] +
            self.Csyms[x_idx, a_idx, e_idx] * self.Rsyms[d_idx, c_idx, x_idx, b_idx] +
            self.Csyms[x_idx, b_idx, e_idx] * self.Rsyms[d_idx, c_idx, a_idx, x_idx]
            for x_idx in range(len(self.coords))
        ])
    
        covDeRdcab = partialD_term + connexion_terms
    
        covDeRdcab = sym.simplify(covDeRdcab)
    
        return covDeRdcab


    def pD_matrices(self, simplify=False):
    
        # Get list of coordinates for partial derivatives
        coords = [f for f in self.coords]
    
        # Convert potential into list of sub-expressions (i.e. seperated by +/-)
        V_list = self.V.args
    
        # Simplify the potential term by term and then sum
        print "computing V"
        if simplify: V_list = [sym.simplify(item) for item in V_list]
        self.V = sum(V_list)
    
        # Repeat process for 1st, 2nd and 3rd order partial derivatives
        for fi in coords:
            print "computing dV", fi
            i = coords.index(fi)
            di = [sym.diff(item, fi) for item in V_list]
            if simplify: di = [sym.simplify(item) for item in di]
            self.pdV[i] = sum(di)
        
            for fii in coords:
                print "computing ddV", fi, fii
                ii = coords.index(fii)
                dii = [sym.diff(item, fii) for item in di]
                if simplify: dii = [sym.simplify(item) for item in dii]
                self.pdpdV[i, ii] = sum(dii)
            
                for fiii in coords:
                    print "computing dddV", fi, fii, fiii
                    iii = coords.index(fiii)
                    diii = [sym.diff(item, fiii) for item in dii]
                    if simplify: diii = [sym.simplify(item) for item in diii]
                    self.pdpdpdV[i, ii, iii] = sum(diii)
                    

    def Pot(self, simplify=False):
        if simplify: self.V = sym.simplify(self.V)
        return self.V


    def dPot(self, a, symbolic_args=True, simplify=False):
        
        if symbolic_args is False:
            a_idx = a
            a = self.coords[a_idx]
        else:
            a_idx = self.coords.index(a)
        
        # partial_aV = sym.diff(self.V, a)
        # if simplify: partial_aV = sym.simplify(partial_aV)
        partial_aV = self.pdV[a_idx]
        
        return partial_aV
   
    
    def ddPot(self, a, b, symbolic_args=True, simplify=False):
    
        if symbolic_args is False:
            a_idx = a
            b_idx = b
            a = self.coords[a_idx]
            b = self.coords[b_idx]
        else:
            a_idx = self.coords.index(a)
            b_idx = self.coords.index(b)
    
        # partial_abV = sym.diff(self.dV[a_idx], b)
        partial_abV = self.pdpdV[a_idx, b_idx]
        
        if self.canonical is False:
            connexion_term = -sum([self.Csyms[k_idx, a_idx, b_idx]*self.dV[k_idx] for k_idx in range(len(self.coords))])
            partial_abV += connexion_term
        
        if simplify: partial_abV = sym.simplify(partial_abV)
        
        return partial_abV
    
    
    def dddPot(self, a, b, c, symbolic_args=True, simplify=False):
        
        if symbolic_args is False:
            a_idx = a
            b_idx = b
            c_idx = c
            a = self.coords[a_idx]
            b = self.coords[b_idx]
            c = self.coords[c_idx]
        else:
            a_idx = self.coords.index(a)
            b_idx = self.coords.index(b)
            c_idx = self.coords.index(c)
            
        # partial_term = sym.diff(self.ddV[b_idx, c_idx], a)
        partial_term = self.pdpdpdV[a_idx, b_idx, c_idx]
        
        if self.canonical is False:
            
            partial_term += -sum([sym.diff(self.Csyms[k_idx, b_idx, c_idx] * self.pdV[k_idx], a)
                                  for k_idx in range(len(self.coords))])
            
            connexion_term = -sum([self.Csyms[d_idx, a_idx, b_idx]*self.ddV[d_idx, c_idx] +
                                   self.Csyms[d_idx, a_idx, c_idx]*self.ddV[d_idx, b_idx]
                                   for d_idx in range(len(self.coords))])
            
            total = partial_term + connexion_term
            
        else: total = partial_term
        
        if simplify: total = sym.simplify(total)
        
        return total

    
    def precompute(self):
        
        max_idx = len(self.coords)
        
        print "-- Precomputing Christoffel symbols"
        
        for a in range(max_idx):
            for b in range(max_idx):
                for c in range(max_idx):
                    self.Csyms[a, b, c] = self.Christoffel(a, b, c, symbolic_args=False)
        
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
            
            self.Rsyms[i, j, k, l] = Rijkl
            
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
        
        rnf = range(max_idx)
        
        print "-- Computing covariant derivative of Riemann tensor"
        
        for i in rnf:
            for j in rnf:
                for k in rnf:
                    for l in rnf:
                        for m in rnf:
                            self.covDRsyms[i, j, k, l, m] = self.CovDRiemann(i, j, k, l, m, symbolic_args=False)
                            
        print "-- Computing partial derivative matrices"
        self.pD_matrices()
        
        print "-- Computing 1st potential derivatives"
        c = 1
        for i in rnf:
            print "computing: {} / {}".format(c, len(self.coords))
            self.dV[i] = self.dPot(i, symbolic_args=False, simplify=False)
            c+=1
            
        print "-- Computing 2nd potential derivatives"
        c=1
        for i in rnf:
            for j in rnf:
                print "computing: {} / {}".format(c, len(self.coords)*len(self.coords))
                c+=1
                self.ddV[i, j] = self.ddPot(i, j, symbolic_args=False, simplify=False)
        
        print "-- Computing 3rd potential derivatives"
        c=1
        for i in rnf:
            for j in rnf:
                for k in rnf:
                    print "computing: {} / {}".format(c, len(self.coords)*len(self.coords)*len(self.coords))
                    c+=1
                    self.dddV[i, j, k] = self.dddPot(i, j, k, symbolic_args=False, simplify=False)

    
    
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
# s = [
# 	2*p[0]*p[1]**4*(-27*p[0]*p[1]**4/(32*sym.pi**2*p[3]**4*f[0]**4) + 1),
# 	p[3]**2*f[0]**2*(2*p[0]*p[1]**4 + p[4])/3,
# 	9*p[5]*p[2]*p[3]**2*f[0]**2*(2*p[0]*p[1]**4 + p[4])*sym.cos(f[2])/(4*sym.pi**(3/2)),
# 	3*sym.sqrt(15)*p[6]*p[2]*p[3]**2*f[0]**(-2 + 2*sym.sqrt(10))*(2*p[0]*p[1]**4 + p[4])*(3*sym.cos(f[2])**2 - 1)/(8*sym.pi**(3/2)),
# 	9*p[7]*p[2]*p[3]**2*f[0]**2*(2*p[0]*p[1]**4 + p[4])*sym.cos(f[1])/(4*sym.pi**(3/2)),
# 	9*sym.sqrt(3)*p[8]*p[2]*p[3]**2*f[0]**(-2 + 2*sym.sqrt(7))*(2*p[0]*p[1]**4 + p[4])*sym.cos(f[1])*sym.cos(f[2])/(4*sym.pi**(3/2)),
# 	9*sym.sqrt(5)*p[9]*p[2]*p[3]**2*f[0]**(-2 + 2*sym.sqrt(13))*(2*p[0]*p[1]**4 + p[4])*(3*sym.cos(f[2])**2 - 1)*sym.cos(f[1])/(8*sym.pi**(3/2)),
# 	3*sym.sqrt(15)*p[10]*p[2]*p[3]**2*f[0]**(-2 + 2*sym.sqrt(10))*(2*p[0]*p[1]**4 + p[4])*(3*sym.cos(f[1])**2 - 1)/(8*sym.pi**(3/2)),
# 	9*sym.sqrt(5)*p[11]*p[2]*p[3]**2*f[0]**(-2 + 2*sym.sqrt(13))*(2*p[0]*p[1]**4 + p[4])*(3*sym.cos(f[1])**2 - 1)*sym.cos(f[2])/(8*sym.pi**(3/2)),
# 	9*sym.sqrt(3)*p[2]*p[3]**2*f[0]**3*(-p[22]*sym.sin(f[3] - f[4] + f[5]) + p[12]*sym.cos(f[3] - f[4] + f[5]))*(2*p[0]*p[1]**4 + p[4])*sym.sin(f[1])*sym.sin(f[2])*sym.tan(f[2]/2)/(8*sym.pi**(3/2)*sym.tan(f[1]/2)),
# 	9*sym.sqrt(5)*p[2]*p[3]**2*f[0]**5*(-p[23]*sym.sin(f[3] - f[4] + f[5]) + p[13]*sym.cos(f[3] - f[4] + f[5]))*(2*p[0]*p[1]**4 + p[4])*(sym.sin(f[2]) + sym.sin(2*f[2]))*sym.sin(f[1])*sym.tan(f[2]/2)/(8*sym.pi**(3/2)*sym.tan(f[1]/2)),
# 	9*sym.sqrt(5)*p[2]*p[3]**2*f[0]**5*(p[24]*sym.sin(f[3] - f[4] + f[5]) - p[14]*sym.cos(f[3] - f[4] + f[5]))*(2*p[0]*p[1]**4 + p[4])*(sym.sin(f[1]) - sym.sin(2*f[1]))*sym.sin(f[2])*sym.tan(f[2]/2)/(8*sym.pi**(3/2)*sym.tan(f[1]/2)),
# 	3*sym.sqrt(3)*p[2]*p[3]**2*f[0]**(3/2)*(-p[25]*sym.sin(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2) + p[15]*sym.cos(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2))*(2*p[0]*p[1]**4 + p[4])*sym.sqrt(sym.sin(f[1]))*sym.sqrt(sym.sin(f[2]))*sym.sqrt(sym.tan(f[2]/2))/(2*sym.pi**(3/2)*sym.sqrt(sym.tan(f[1]/2))),
# 	3*sym.sqrt(6)*p[2]*p[3]**2*f[0]**(7/2)*(-p[26]*sym.sin(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2) + p[16]*sym.cos(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2))*(2*p[0]*p[1]**4 + p[4])*(3*sym.cos(f[2]) + 1)*sym.sqrt(sym.sin(f[1]))*sym.sqrt(sym.sin(f[2]))*sym.sqrt(sym.tan(f[2]/2))/(4*sym.pi**(3/2)*sym.sqrt(sym.tan(f[1]/2))),
# 	9*p[2]*p[3]**2*f[0]**(-2 + sym.sqrt(241)/2)*(-p[27]*sym.sin(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2) + p[17]*sym.cos(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2))*(2*p[0]*p[1]**4 + p[4])*(-5*sym.sin(f[2])**2 + 2*sym.cos(f[2]) + 4)*sym.sqrt(sym.sin(f[1]))*sym.sqrt(sym.sin(f[2]))*sym.sqrt(sym.tan(f[2]/2))/(4*sym.pi**(3/2)*sym.sqrt(sym.tan(f[1]/2))),
# 	3*sym.sqrt(6)*p[2]*p[3]**2*f[0]**(7/2)*(-p[28]*sym.sin(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2) + p[18]*sym.cos(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2))*(2*p[0]*p[1]**4 + p[4])*(3*sym.cos(f[1]) - 1)*sym.sqrt(sym.sin(f[1]))*sym.sqrt(sym.sin(f[2]))*sym.sqrt(sym.tan(f[2]/2))/(4*sym.pi**(3/2)*sym.sqrt(sym.tan(f[1]/2))),
# 	3*sym.sqrt(3)*p[2]*p[3]**2*f[0]**(-2 + sym.sqrt(193)/2)*(-p[29]*sym.sin(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2) + p[19]*sym.cos(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2))*(2*p[0]*p[1]**4 + p[4])*(3*sym.cos(f[1]) - 1)*(3*sym.cos(f[2]) + 1)*sym.sqrt(sym.sin(f[1]))*sym.sqrt(sym.sin(f[2]))*sym.sqrt(sym.tan(f[2]/2))/(4*sym.pi**(3/2)*sym.sqrt(sym.tan(f[1]/2))),
# 	3*sym.sqrt(3)*p[2]*p[3]**2*f[0]**(9/2)*(-p[30]*sym.sin(sym.Rational(3, 2)*f[3] - sym.Rational(3, 2)*f[4] + 3*f[5]/2) + p[20]*sym.cos(sym.Rational(3, 2)*f[3] - sym.Rational(3, 2)*f[4] + 3*f[5]/2))*(2*p[0]*p[1]**4 + p[4])*sym.sin(f[1])**(3/2)*sym.sin(f[2])**(3/2)*sym.tan(f[2]/2)**(3/2)/(4*sym.pi**(3/2)*sym.tan(f[1]/2)**(3/2)),
# 	9*p[2]*p[3]**2*f[0]**(-2 + sym.sqrt(241)/2)*(p[31]*sym.sin(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2) - p[21]*sym.cos(f[3]/2 - sym.Rational(1, 2)*f[4] + f[5]/2))*(2*p[0]*p[1]**4 + p[4])*(5*sym.sin(f[1])**2 + 2*sym.cos(f[1]) - 4)*sym.sqrt(sym.sin(f[1]))*sym.sqrt(sym.sin(f[2]))*sym.sqrt(sym.tan(f[2]/2))/(4*sym.pi**(3/2)*sym.sqrt(sym.tan(f[1]/2))),
# 	0
# ]
# V = sum(s)
#
# V_list = V.args
#
# coords = [fi for fi in f]
#
# # Repeat process for 1st, 2nd and 3rd order partial derivatives
# for fi in coords:
#     print "computing dV", fi
#     i = coords.index(fi)
#     print "-- Differentiating"
#     di_dif = [sym.diff(sym.expand_trig(item), fi) for item in V_list]
#     print "-- Simplifying"
#     di_sim = [sym.simplify(item) for item in di_dif]
#     print sum(di_sim)
#
#     for fii in coords:
#         print "computing ddV", fi, fii
#         ii = coords.index(fii)
#         print "-- Differentiating"
#         dii_dif = [sym.diff(sym.expand_trig(item), fii) for item in di_dif]
#         print "-- Simplifying"
#         dii_sim = [sym.simplify(item) for item in dii_dif]
#         print sum(dii_sim)
#
#         for fiii in coords:
#             print "computing dddV", fi, fii, fiii
#             iii = coords.index(fiii)
#             print "-- Differentiating"
#             diii_dif = [sym.diff(sym.expand_trig(item), fiii) for item in dii_dif]
#             print "-- Simplifying"
#             diii_sim = [sym.simplify(item) for item in diii_dif]
#             print sum(diii_sim)

