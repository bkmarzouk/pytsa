import sympy as sym
import numpy as np
import scipy as sp
from sympy.matrices import Matrix

from PyTransport.idp import OperationMatrix


class curvatureObject:

    def __init__(self, metric, coords, potential, params=None, simpleGeometric=True, simplePotentials=False):

        print("-- Initializing curvature object")

        # Distinguish coordinates from other symbols which may be in the metric, e.g. mass parameters
        self.coords = list(coords)

        # Get number of coords and range of coord indices
        self.nC = len(self.coords)
        self.nCr = range(self.nC)

        # Keep track of parameters, if we wish
        if params is not None:
            self.params = list(params)
        else:
            self.params = None

        # If canonical, assume Kronecker-delta, i.e. identity matrix representation for the metric
        if metric == "canonical" or metric is None or metric == 0:
            self.metric = sym.eye(self.nC)  # G_{IJ} = \delta_{IJ}
            self.canonical = True
        else:

            # Check that metric representation is square matrix
            assert metric.is_symmetric(), "Metric must be represented by a square NxN matrix."
            self.metric = metric
            self.canonical = False

        print("-- Inverting metric")

        # Compute inverse metric
        self.metric_inverse = self.metric.inv()

        print("-- Simplifying inverse metric")


        for a in self.nCr:
            for b in self.nCr:
                if b >= a:
                    self.metric_inverse[a, b] = sym.simplify(self.metric_inverse[a, b])
                else:
                    self.metric_inverse[a, b] = self.metric_inverse[b, a]

        print("-- Computing remaining metric index combinations")

        self.metric_Ud = Matrix.zeros(len(self.coords))  # Up down
        self.metric_dU = Matrix.zeros(len(self.coords))  # down Up
        for a in self.nCr:
            for b in self.nCr:
                self.metric_Ud[a, b] = sym.simplify(
                    sum([self.metric_inverse[a, c] * self.metric[c, b] for c in self.nCr]))
                self.metric_dU[a, b] = sym.simplify(
                    sum([self.metric[a, c] * self.metric_inverse[c, b] for c in self.nCr]))

        assert self.metric_Ud.is_symmetric(), 0
        assert self.metric_dU.is_symmetric(), 0

        # Consistency check: Product of G * G^-1 = Identity
        identity_matrix = self.metric * self.metric_inverse
        identity_matrix = sym.simplify(identity_matrix)
        assert identity_matrix == sym.eye(self.nC), "G * G^-1 != I: {}".format(identity_matrix)

        # We now generate the matrix format that is ready by _PyTransport backend

        # Construct half of combined matrix
        hg1 = np.empty((self.nC, 2 * self.nC), dtype=object)
        hg2 = np.empty((self.nC, 2 * self.nC), dtype=object)
        for i in self.nCr:
            for j in self.nCr:
                hg1[i, j] = self.metric_inverse[i, j]
                hg1[i, self.nC + j] = self.metric_Ud[i, j]
                hg2[i, j] = self.metric_dU[i, j]
                hg2[i, self.nC + j] = self.metric[i, j]

        combined_matrix_flattened = np.vstack((hg1, hg2)).flatten()
        sympy_matrix_flattened = sym.symarray("G", 2 * self.nC * 2 * self.nC)

        for ii in range(2 * self.nC * 2 * self.nC):
            sympy_matrix_flattened[ii] = sym.simplify(combined_matrix_flattened[ii])

        self.G_array = sympy_matrix_flattened

        # Indicate whether to simplify further symbolic expressions (though we always demand the metric identity)
        self.simpleG = simpleGeometric
        self.simpleV = simplePotentials

        # Initialize zeros array to store Christoffel symbols
        self._christoffelSymbols = np.zeros(
            (len(coords), len(coords), len(coords)), dtype=object
        )

        # Initialize zeros array to store Riemann tensors
        self._riemannTensor = np.zeros(
            (len(coords), len(coords), len(coords), len(coords)), dtype=object
        )

        # Initialize zeros array to store covariant derivative(s) of Riemann tensor
        self._riemannTensorCovD = np.zeros(
            (len(coords), len(coords), len(coords), len(coords), len(coords)), dtype=object
        )

        # Initialize arrays to store partial derivatives of the potential
        self._pdV = np.zeros((len(coords)), dtype=object)
        self._pdpdV = np.zeros((len(coords), len(coords)), dtype=object)
        self._pdpdpdV = np.zeros((len(coords), len(coords), len(coords)), dtype=object)

        # Initialize arrays to store covariant derivatives of the potential
        print("-- Simplifying initial potential")
        self.V = potential if simplePotentials is False else sym.simplify(potential)
        self.dV = np.zeros((len(coords)), dtype=object)
        self.ddV = np.zeros((len(coords), len(coords)), dtype=object)
        self.dddV = np.zeros((len(coords), len(coords), len(coords)), dtype=object)

        print("-- Precomputing symbolic expressions")


        self.computed = False
        self.precompute()

        # update flag: calls to functions will now load precompute results
        self.computed = True

        print("-- Curvature object constructed")


    def _handleCoordIn(self, c):
        """ Maps coordinate input arg -> coordinate symbol, coordinate index """

        # Check whether c is an appropriate integer
        valid_int = type(c) in [int, np.int64, np.int32, np.int16] and c < self.nC

        # Return coordinate symbol, index, depending on c type
        if valid_int:
            c = int(c)
            return self.coords[c], c
        elif c in self.coords:
            return c, self.coords.index(c)
        else:
            raise AttributeError("Cannot interpret coordinate arg: {} <type>: {}".format(c, type(c)))

    def _computePartialDerivatives(self):
        """ Compute partial derivatives of the potential up until the 3rd derivative """

        c = 1
        NN = self.nC ** 3 + self.nC ** 2 + self.nC

        for ii in self.nCr:
            c += 1

            # Compute first partial derivative of V
            fii = self.coords[ii]

            _pdV = sym.diff(self.V, fii)

            if self.simpleV: _pdV = sym.simplify(_pdV)

            self._pdV[ii] = _pdV

            print("partial derivative {} / {}".format(c, NN))

            for jj in self.nCr:

                if jj >= ii:

                    # Compute second partial derivative of V
                    fjj = self.coords[jj]

                    _pdpdV = sym.diff(self._pdV[ii], fjj)

                    if self.simpleV: _pdpdV = sym.simplify(_pdpdV)

                    self._pdpdV[ii, jj] = _pdpdV

                else:
                    self._pdpdV[ii, jj] = self._pdpdV[jj, ii]

                c += 1
                print("partial derivative {} / {}".format(c, NN))

                for kk in self.nCr:

                    if kk >= jj:

                        fkk = self.coords[kk]

                        _pdpdpdV = sym.diff(self._pdpdV[ii, jj], fkk)

                        if self.simpleV: _pdpdpdV = sym.simplify(_pdpdpdV)

                        self._pdpdpdV[ii, jj, kk] = _pdpdpdV

                    else:
                        self._pdpdpdV[ii, jj, kk] = self._pdpdpdV[ii, kk, jj]

                    c += 1
                    print("partial derivative {} / {}".format(c, NN))

    def getChristoffel(self, a, b, c):
        """ Compute Christoffel symbol of the second kind, \Gamma^a_{bc} """

        if self.canonical is True: return 0

        a, a_idx = self._handleCoordIn(a)
        b, b_idx = self._handleCoordIn(b)
        c, c_idx = self._handleCoordIn(c)

        if self.computed is True: return self._christoffelSymbols[a_idx, b_idx, c_idx]

        gamma = 0

        # Iterate over the repeated index, d
        for d_idx in self.nCr:
            d = self.coords[d_idx]

            gamma += sym.Rational(1, 2) * self.metric_inverse[a_idx, d_idx] * sum([
                sym.diff(self.metric[d_idx, c_idx], b),  # partial_b g_{dc} +
                sym.diff(self.metric[b_idx, d_idx], c),  # partial_c g_{bd} -
                -sym.diff(self.metric[b_idx, c_idx], d)  # partial_d g_{bc}
            ])

        # Return Christoffel symbol with simplification
        if self.simpleG: gamma = sym.simplify(gamma)

        return gamma

    def getRiemannUp(self, f, c, a, b):
        """ Compute Riemann tensor with first index up R^{f}_{cab} """

        if self.canonical is True: return 0

        f, f_idx = self._handleCoordIn(f)
        c, c_idx = self._handleCoordIn(c)
        a, a_idx = self._handleCoordIn(a)
        b, b_idx = self._handleCoordIn(b)

        term1 = sym.diff(self._christoffelSymbols[f_idx, b_idx, c_idx], a)
        term2 = -sym.diff(self._christoffelSymbols[f_idx, a_idx, c_idx], b)
        term3 = sum([self._christoffelSymbols[f_idx, a_idx, d_idx] * self._christoffelSymbols[d_idx, b_idx, c_idx] -
                     self._christoffelSymbols[f_idx, b_idx, d_idx] * self._christoffelSymbols[d_idx, a_idx, c_idx]
                     for d_idx in self.nCr])

        Rfcab = term1 + term2 + term3

        if self.simpleG: Rfcab = sym.simplify(Rfcab)

        return Rfcab

    def getRiemann(self, d, c, a, b):
        """ Compute Riemann tensor with all indices down R_{dcab} """

        if self.canonical is True: return 0

        d, d_idx = self._handleCoordIn(d)
        c, c_idx = self._handleCoordIn(c)
        a, a_idx = self._handleCoordIn(a)
        b, b_idx = self._handleCoordIn(b)

        if self.computed is True: return self._riemannTensor[d_idx, c_idx, a_idx, b_idx]

        Rdcab = sum([self.metric[d_idx, f_idx] * self.getRiemannUp(
            f_idx, c_idx, a_idx, b_idx) for f_idx in self.nCr])

        if self.simpleG: Rdcab = sym.simplify(Rdcab)

        return Rdcab

    def getCovDRiemann(self, e, d, c, a, b):
        """ Compute the covariant derivative of the Riemann tensor \nabla_e R_{dcab}"""

        if self.canonical is True: return 0

        e, e_idx = self._handleCoordIn(e)
        d, d_idx = self._handleCoordIn(d)
        c, c_idx = self._handleCoordIn(c)
        a, a_idx = self._handleCoordIn(a)
        b, b_idx = self._handleCoordIn(b)

        if self.computed is True: return self._riemannTensorCovD[e_idx, d_idx, c_idx, a_idx, b_idx]

        partialD_term = sym.diff(self._riemannTensor[d_idx, c_idx, a_idx, b_idx], e)

        connexion_terms = -sum([
            self._christoffelSymbols[x_idx, d_idx, e_idx] * self._riemannTensor[x_idx, c_idx, a_idx, b_idx] +
            self._christoffelSymbols[x_idx, c_idx, e_idx] * self._riemannTensor[d_idx, x_idx, a_idx, b_idx] +
            self._christoffelSymbols[x_idx, a_idx, e_idx] * self._riemannTensor[d_idx, c_idx, x_idx, b_idx] +
            self._christoffelSymbols[x_idx, b_idx, e_idx] * self._riemannTensor[d_idx, c_idx, a_idx, x_idx]
            for x_idx in self.nCr
        ])

        covDeRdcab = partialD_term + connexion_terms

        if self.simpleG: covDeRdcab = sym.simplify(covDeRdcab)

        return covDeRdcab

    def getCovDV(self, a):
        """ Compute 1st covariant derivative of the potential: \nabla_a V = \partial_a V """

        a, a_idx = self._handleCoordIn(a)

        if self.computed is True: return self.dV[a_idx]

        partial_aV = self._pdV[a_idx]

        if self.simpleV: partial_aV = sym.simplify(partial_aV)

        return partial_aV

    def getCovDDV(self, a, b):
        """ Compute 2nd covariant derivative of the potential, \nabla_a \nabla_b V """

        a, a_idx = self._handleCoordIn(a)
        b, b_idx = self._handleCoordIn(b)

        if self.computed is True: return self.ddV[a_idx, b_idx]

        partial_abV = self._pdpdV[a_idx, b_idx]

        if self.canonical is False:
            connexion_term = -sum(
                [self._christoffelSymbols[k_idx, a_idx, b_idx] * self.dV[k_idx] for k_idx in self.nCr])
            partial_abV += connexion_term

        if self.simpleV: partial_abV = sym.simplify(partial_abV)

        return partial_abV

    def getCovDDDV(self, a, b, c):
        """ Compute 3rd covariant derivative of the potential, \nabla_a \nabla_b \nabla_c V """

        a, a_idx = self._handleCoordIn(a)
        b, b_idx = self._handleCoordIn(b)
        c, c_idx = self._handleCoordIn(c)

        if self.computed is True: return self.dddV[a_idx, b_idx, c_idx]

        partial_abcV = self._pdpdpdV[a_idx, b_idx, c_idx]

        if self.canonical is False:

            t1 = -sum([sym.diff(self._christoffelSymbols[k_idx, b_idx, c_idx] * self._pdV[k_idx], a)
                       for k_idx in self.nCr])

            t2 = -sum([self._christoffelSymbols[d_idx, a_idx, b_idx] * self.ddV[d_idx, c_idx] +
                       self._christoffelSymbols[d_idx, a_idx, c_idx] * self.ddV[d_idx, b_idx]
                       for d_idx in self.nCr])

            total = partial_abcV + t1 + t2

        else:
            total = partial_abcV

        if self.simpleV: total = sym.simplify(total)

        return total

    def precompute(self):
        """ Precomputes Christoffel Symbols, Riemann Tensor and derivatives of the potential"""

        print("-- Precomputing Christoffel symbols")

        for a in self.nCr:
            for b in self.nCr:
                for c in self.nCr:
                    if c >= b:
                        self._christoffelSymbols[a, b, c] = self.getChristoffel(a, b, c)
                    else:
                        self._christoffelSymbols[a, b, c] = self._christoffelSymbols[a, c, b]


        print("-- Precomputing Riemann tensor components")

        # Get operation matrix, i.e. identify fundamental components and associated symmetries
        OM = OperationMatrix(self.nC)

        assert np.shape(OM) == np.shape(self._riemannTensor), "Operation matrix incompatible with Riemann matrix"

        # Get location of independent elements to compute & symmetrically related components
        calc_locs = np.where(OM == "calc")
        symm_locs = np.where(np.logical_and(OM != "calc", OM != 0))

        # Transpose into rows of i, j, k, l matrix coordinates / tensor indices
        calc_ijkl = np.asarray(calc_locs).T
        symm_ijkl = np.asarray(symm_locs).T

        # Count fundamental component computations dynamically
        c = 0
        m = len(calc_ijkl)

        # Compute independent components
        for ijkl in calc_ijkl:
            i, j, k, l = ijkl

            assert self._riemannTensor[i, j, k, l] == 0, "Riemann component already has assigned value! {}".format(
                self._riemannTensor[i, j, k, l]
            )

            Rijkl = self.getRiemann(i, j, k, l)

            self._riemannTensor[i, j, k, l] = Rijkl

            c += 1


        # Assign symmetrically related components
        for ijkl in symm_ijkl:

            i, j, k, l = ijkl

            assert self._riemannTensor[i, j, k, l] == 0, "Riemann component already has assigned value! {}".format(
                self._riemannTensor[i, j, k, l]
            )

            # Unpack the string that points to the fundamental Matrix element(s) and multiplicative sign
            str_pointer = OM[i, j, k, l]

            if len(str_pointer) == 5:
                sign, p, q, r, s = str_pointer

                Rpqrs = self._riemannTensor[int(p), int(q), int(r), int(s)]

                if sign == "+":
                    self._riemannTensor[i, j, k, l] += Rpqrs
                elif sign == "-":
                    self._riemannTensor[i, j, k, l] += -Rpqrs
                else:
                    raise ValueError("Unrecognized sign value: {}".format(sign))

            elif len(str_pointer) == 10:
                sign1, p, q, r, s, sign2, t, u, v, w = str_pointer

                Rpqrs = self._riemannTensor[int(p), int(q), int(r), int(s)]
                Rtuvw = self._riemannTensor[int(t), int(u), int(v), int(w)]

                if sign1 == "+":
                    self._riemannTensor[i, j, k, l] += Rpqrs
                elif sign1 == "-":
                    self._riemannTensor[i, j, k, l] += -Rpqrs
                else:
                    raise ValueError("Unrecognized sign value: {}".format(sign1))

                if sign2 == "+":
                    self._riemannTensor[i, j, k, l] += Rtuvw
                elif sign2 == "-":
                    self._riemannTensor[i, j, k, l] += -Rtuvw
                else:
                    raise ValueError("Unrecognized sign value: {}".format(sign2))

            else:
                raise ValueError("Unrecognized symmetry relation: {}".format(str_pointer))

        print("-- Computing covariant derivative of Riemann tensor")

        c = 1

        for i in self.nCr:
            for j in self.nCr:
                for k in self.nCr:
                    for l in self.nCr:
                        for m in self.nCr:
                            print("{} / {}".format(c, self.nC ** 5))
                            self._riemannTensorCovD[i, j, k, l, m] = self.getCovDRiemann(i, j, k, l, m)
                            c += 1

        print("-- Computing partial derivatives")
        self._computePartialDerivatives()

        print("-- Computing potential 1st derivatives")
        c = 1
        for i in self.nCr:
            print("computing: {} / {}".format(c, len(self.coords)))
            self.dV[i] = self.getCovDV(i)
            c += 1

        print("-- Computing potential 2nd derivatives")
        c = 1
        for i in self.nCr:
            for j in self.nCr:
                print("potential derivative {} / {}".format(c, len(self.coords) * len(self.coords)))
                c += 1
                if j >= i:
                    self.ddV[i, j] = self.getCovDDV(i, j)
                else:
                    self.ddV[i, j] = self.ddV[j, i]

        print("-- Computing potential 3rd derivatives")
        c = 1
        for i in self.nCr:
            for j in self.nCr:
                for k in self.nCr:
                    print("computing: {} / {}".format(c, len(self.coords) * len(self.coords) * len(self.coords)))
                    c += 1
                    if k >= j:
                        self.dddV[i, j, k] = self.getCovDDDV(i, j, k)
                    else:
                        self.dddV[i, j, k] = self.dddV[i, k, j]

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
# a = curvatureObject(G, f, V, p)
