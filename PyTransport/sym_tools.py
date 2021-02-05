import numpy as np
import numpy as np
import sympy as sym
import os
import cache_tools


def OperationMatrix(N):
    """
    Computes matrix of operations for independent components of the Riemann tensor

    :param N: Number of dimensions (fields)
    :return: Array of operatios defined by string kwargs
    """

    # Initialize N^4 array
    OM = np.zeros([N, N, N, N], dtype=object)

    # Initialize counters for non-zero independent components
    idp_ijij = 0  # R_ijij
    idp_ijki = 0  # R_ijki
    idp_ijkl = 0  # R_ijkl

    # Track components which are reducible via the Bianchi cyclic identity
    bterms = []

    # Iterate over N coords
    for i in range(N):

        jrange = list(range(N))
        jrange.remove(i)

        # Iterate over (N-1) coords
        for j in jrange:

            # Construct independent components to have i<j
            if i < j:
                OM[i, j, i, j] = "calc"
                idp_ijij += 1

                # Prescribe dependent components in terms of fundamental component * +/-
                subs = [i, j, i, j]
                pijij = "+{}{}{}{}".format(*subs)
                mijij = "-{}{}{}{}".format(*subs)

                # Skewing of first or last pair incurs a minus sign
                OM[j, i, i, j] = mijij
                OM[i, j, j, i] = mijij

                # Skewing both pairs of terms leaves the sign unchanged
                OM[j, i, j, i] = pijij

                # Since symmetric over index pairs, exchange symmetry need not be considered

            krange = list(range(N))
            krange.remove(i)
            krange.remove(j)

            # Iterate over (N-2) coords
            for k in krange:

                # Choose independent solutions to have k < j
                if k < j:
                    OM[i, j, k, i] = "calc"
                    idp_ijki += 1

                    # Prescribe dependent components in terms of fundamental component & +/-
                    subs = [i, j, k, i]
                    pijki = "+{}{}{}{}".format(*subs)
                    mijki = "-{}{}{}{}".format(*subs)

                    # Skewing first and last pair incurs a minus sign
                    OM[j, i, k, i] = mijki
                    OM[i, j, i, k] = mijki
                    OM[j, i, i, k] = pijki

                    # Since first and last pairs carry different index pairs, consider exchange symmetries
                    OM[k, i, i, j] = pijki
                    OM[k, i, j, i] = mijki
                    OM[i, k, i, j] = mijki
                    OM[i, k, j, i] = pijki

                lrange = list(range(N))
                lrange.remove(i)
                lrange.remove(j)
                lrange.remove(k)

                # iterate over (N-3) coords
                for l in lrange:

                    # Choose irreducible representation to be those with i < j & k < l in first & 2nd index pairs.
                    # To avoid double counting under exchange symmetry, we items with lowest index in first index pair
                    if i < j and k < l and min(i, j) <= min(k, l):

                        OM[i, j, k, l] = "calc"
                        idp_ijkl += 1

                        # Prescribe dependent components in terms of fundamental component & +/-
                        subs = [i, j, k, l]
                        pijkl = "+{}{}{}{}".format(*subs)
                        mijkl = "-{}{}{}{}".format(*subs)

                        # Skewing first and last pair incurs a minus sign
                        OM[j, i, k, l] = mijkl
                        OM[i, j, l, k] = mijkl
                        OM[j, i, l, k] = pijkl

                        # Since first and last pairs carry different index paris, consider exchange symmetry
                        OM[k, l, i, j] = pijkl
                        OM[k, l, j, i] = mijkl
                        OM[l, k, i, j] = mijkl
                        OM[l, k, j, i] = pijkl

                        bterms.append([i, j, k, l])

                    else:
                        pass

    c = 0

    # Copy terms which could be related via the Bianchi identity
    bterms_copy = bterms

    def exch(i, j, k, l):
        return [k, l, i, j]

    for item in bterms_copy:
        i, j, k, l = item

        # Cyclic permute item in bterms
        b1 = [i, k, l, j]
        b2 = [i, l, j, k]

        # Cyclic permute item in bterms and skew indices
        b1s = [i, k, j, l]
        b2s = [i, l, k, j]

        # We have then the following possible representations for the Riemann tensor
        con1 = (b1 in bterms or exch(*b1) in bterms) and (b2 in bterms or exch(*b2) in bterms)
        con2 = (b1s in bterms or exch(*b1s) in bterms) and (b2s in bterms or exch(*b2s) in bterms)
        con3 = (b1 in bterms or exch(*b1) in bterms) and (b2s in bterms or exch(*b2s) in bterms)
        con4 = (b1s in bterms or exch(*b1s) in bterms) and (b2 in bterms or exch(*b2) in bterms)

        if con1:
            bterms.remove(item);
            idp_ijkl += -1
            OM[i, j, k, l] = "-{}{}{}{}".format(*b1) + "-{}{}{}{}".format(*b2)

        elif con2:
            bterms.remove(item);
            idp_ijkl += -1
            OM[i, j, k, l] = "+{}{}{}{}".format(*b1s) + "+{}{}{}{}".format(*b2s)

        elif con3:
            bterms.remove(item);
            idp_ijkl += -1
            OM[i, j, k, l] = "-{}{}{}{}".format(*b1) + "+{}{}{}{}".format(*b2s)

        elif con4:
            bterms.remove(item);
            idp_ijkl += -1
            OM[i, j, k, l] = "+{}{}{}{}".format(*b1s) + "-{}{}{}{}".format(*b2)
        else:
            pass

    # Enforce analytic consistency check on the number of independent items
    assert 2 * idp_ijij == N * (N - 1), str(idp_ijij) + "/" + str(N * (N - 1) / 2)
    assert 2 * idp_ijki == N * (N - 1) * (N - 2), str(idp_ijki) + "/" + str(N * (N - 1) * (N - 2) / 2)
    assert 12 * idp_ijkl == N * (N - 1) * (N - 2) * (N - 3), str(idp_ijkl) + "/" + str(
        N * (N - 1) * (N - 2) * (N - 3) / 12)

    print("-- Independent components: {}".format(idp_ijij + idp_ijki + idp_ijkl))

    # Return matrix
    return OM


def FieldSpaceSym(n_fields: int, n_params: int, metric: sym.Matrix, recache=False, simplify=True):
    try:
        fs_derived = cache_tools.load_fmet(metric, delete=recache)
        print("-- Derived field space calculations loaded from cache")

    except OSError:
        print("-- Performing field space calculations ...")
        fs_derived = _FieldSpaceSym(n_fields, n_params, metric, simplify)

    return fs_derived

def PotentialSym(n_fields: int, n_params: int, potential: sym.Expr, recache=False, simplify=True):
    try:
        pot_derived = cache_tools.load_pot(potential, delete=recache)
        print("-- Derived potential calculations loaded from cache")

    except OSError:
        print("-- Performing potential calculations ...")
        pot_derived = _PotentialSym(n_fields, n_params, potential, simplify)

    return pot_derived


class _FieldSpaceSym(object):
    """ Handles symbolic calculations relevant to the field space metric"""

    def __init__(self, n_fields: int, n_params: int, metric: sym.Matrix, simplify=True):

        self.nf = n_fields
        self.nf_range = list(range(n_fields))
        self.f_syms = sym.symarray("f", n_fields)  # field symbols
        self.p_syms = sym.symarray("f", n_params) if n_params > 0 else None  # param symbols

        self.simplify = simplify

        self._configure_metric(metric)
        self._configure_christoffel_symbols()
        self._configure_riemann_tensor()

        cache_tools.cache_fmet(metric, self)

    def _coord2symb(self, coord):

        if coord < 0 or coord > self.nf - 1:
            raise IndexError("coord number {} not defined".format(coord))

        return self.f_syms[coord]

    def _symb2coord(self, symb):
        for coord, s in enumerate(self.f_syms):
            if symb == s:
                return coord

        raise IndexError("field symbol {} not defined".format(symb))

    def get_symb_and_coord(self, symb_or_coord):

        if isinstance(symb_or_coord, np.int64):
            return self._coord2symb(symb_or_coord), symb_or_coord
        elif isinstance(symb_or_coord, np.int32):
            return self._coord2symb(symb_or_coord), symb_or_coord
        elif isinstance(symb_or_coord, np.int16):
            return self._coord2symb(symb_or_coord), symb_or_coord
        elif isinstance(symb_or_coord, int):
            return self._coord2symb(symb_or_coord), symb_or_coord
        else:
            return symb_or_coord, self._symb2coord(symb_or_coord)

    def _configure_metric(self, metric):

        kronecker = sym.eye(self.nf)

        if metric is None:
            self.metric = kronecker
            self.canonical = True
        else:
            assert metric.is_symmetric()
            self.metric = metric
            self.canonical = False

        assert self.metric.shape[0] == self.nf, "Metric dimensions do not match number of fields: {}".format(self.nf)

        self.metric_inverse = metric.inv()

        # Simplify inverse metric, can in general be complicated
        for a in self.nf_range:
            for b in self.nf_range:
                if b >= a:
                    self.metric_inverse[a, b] = sym.simplify(self.metric_inverse[a, b])
                else:
                    self.metric_inverse[b, a] = self.metric_inverse[a, b]  # symmetry

        identity = sym.simplify(self.metric * self.metric_inverse)

        assert identity == kronecker, "Failed matrix inversion: g g^-1 != I"

        # Pack matrix into format to be read into PyTransport

        self.metric_ud = sym.Matrix.zeros(self.nf)  # Up down indices
        self.metric_du = sym.Matrix.zeros(self.nf)  # down Up indices

        # Construct half of combined matrix
        hg1 = np.empty((self.nf, 2 * self.nf), dtype=object)
        hg2 = np.empty((self.nf, 2 * self.nf), dtype=object)
        for i in self.nf_range:
            for j in self.nf_range:
                hg1[i, j] = self.metric_inverse[i, j]
                hg1[i, self.nf + j] = self.metric_ud[i, j]
                hg2[i, j] = self.metric_du[i, j]
                hg2[i, self.nf + j] = self.metric[i, j]

        combined_matrix_flattened = np.vstack((hg1, hg2)).flatten()
        sympy_matrix_flattened = sym.symarray("G", 2 * self.nf * 2 * self.nf)

        for ii in range(2 * self.nf * 2 * self.nf):
            sympy_matrix_flattened[ii] = sym.simplify(combined_matrix_flattened[ii])

        assert self.metric_ud.is_symmetric(), 0
        assert self.metric_du.is_symmetric(), 0

        self.G_array = sympy_matrix_flattened

    def _compute_christoffel(self, a, b, c):

        if self.canonical:
            christoffel = 0

        else:
            a, a_idx = self.get_symb_and_coord(a)
            b, b_idx = self.get_symb_and_coord(b)
            c, c_idx = self.get_symb_and_coord(c)

            christoffel = 0

            # Iterate over the repeated index, d
            for d_idx, d in enumerate(self.f_syms):
                christoffel += sym.Rational(1, 2) * self.metric_inverse[a_idx, d_idx] * sum([
                    sym.diff(self.metric[d_idx, c_idx], b),  # partial_b g_{dc} +
                    sym.diff(self.metric[b_idx, d_idx], c),  # partial_c g_{bd} -
                    -sym.diff(self.metric[b_idx, c_idx], d)  # partial_d g_{bc}
                ])

        return sym.simplify(christoffel) if self.simplify else christoffel

    def _configure_christoffel_symbols(self):

        print("-- Computing Christoffel symbols")

        if self.canonical:
            self.christoffel_symbols = 0

        else:
            self.christoffel_symbols = np.zeros((self.nf, self.nf, self.nf), dtype=object)

            for a in self.nf_range:
                for b in self.nf_range:
                    for c in self.nf_range:
                        if c >= b:  # symmetry on lower indices
                            self.christoffel_symbols[a, b, c] = self._compute_christoffel(a, b, c)
                        else:
                            self.christoffel_symbols[a, b, c] = self.christoffel_symbols[a, c, b]

    def _compute_mixed_riemann_tensor(self, f, c, a, b):
        """
        Computes element of Riemann tensor with first index raised and rest lower

        :return:R^f_{cab}
        """

        if self.canonical:
            Rfcab = 0
        else:
            f, f_idx = self.get_symb_and_coord(f)
            c, c_idx = self.get_symb_and_coord(c)
            a, a_idx = self.get_symb_and_coord(a)
            b, b_idx = self.get_symb_and_coord(b)

            term1 = sym.diff(self.christoffel_symbols[f_idx, b_idx, c_idx], a)
            term2 = -sym.diff(self.christoffel_symbols[f_idx, a_idx, c_idx], b)
            term3 = sum([self.christoffel_symbols[f_idx, a_idx, d_idx] * self.christoffel_symbols[d_idx, b_idx, c_idx] -
                         self.christoffel_symbols[f_idx, b_idx, d_idx] * self.christoffel_symbols[d_idx, a_idx, c_idx]
                         for d_idx in self.nf_range])

            Rfcab = term1 + term2 + term3

        return sym.simplify(Rfcab) if self.simplify else Rfcab

    def _compute_lowered_riemann_tensor(self, d, c, a, b):
        """
        Computes element of Riemann tensor with all indices lowered

        :return:R_{dcab}
        """

        if self.canonical:
            Rdcab = 0
        else:

            print(d, type(d))
            d, d_idx = self.get_symb_and_coord(d)
            c, c_idx = self.get_symb_and_coord(c)
            a, a_idx = self.get_symb_and_coord(a)
            b, b_idx = self.get_symb_and_coord(b)

            Rdcab = sum([self.metric[d_idx, f_idx] * self._compute_mixed_riemann_tensor(
                f_idx, c_idx, a_idx, b_idx) for f_idx in self.nf_range])

        return sym.simplify(Rdcab) if self.simplify else Rdcab

    def _configure_riemann_tensor(self):

        print("-- Computing Riemann tensor")

        self.riemann_tensor = np.zeros((self.nf, self.nf, self.nf, self.nf), dtype=object)

        om = OperationMatrix(self.nf)

        assert self.riemann_tensor.shape == om.shape, "Operation matrix incompatible shape with Riemann tensor"

        calc_locs = np.where(om == "calc")
        symm_locs = np.where(np.logical_and(om != "calc", om != 0))

        # Transpose into rows of i, j, k, l matrix coordinates / tensor indices
        calc_ijkl = np.asarray(calc_locs).T
        symm_ijkl = np.asarray(symm_locs).T

        # Count fundamental component computations dynamically
        c = 0
        m = len(calc_ijkl)

        # Compute independent components
        for ijkl in calc_ijkl:
            i, j, k, l = ijkl

            assert self.riemann_tensor[i, j, k, l] == 0, "Riemann component already has assigned value! {}".format(
                self.riemann_tensor[i, j, k, l]
            )

            Rijkl = self._compute_lowered_riemann_tensor(i, j, k, l)

            self.riemann_tensor[i, j, k, l] = Rijkl

            c += 1

            print("computed {}/{}".format(c, m))

        # Assign symmetrically related components
        for ijkl in symm_ijkl:

            i, j, k, l = ijkl

            assert self.riemann_tensor[i, j, k, l] == 0, "Riemann component already has assigned value! {}".format(
                self.riemann_tensor[i, j, k, l]
            )

            # Unpack the string that points to the fundamental Matrix element(s) and multiplicative sign
            str_pointer = om[i, j, k, l]

            if len(str_pointer) == 5:
                sign, p, q, r, s = str_pointer

                Rpqrs = self.riemann_tensor[int(p), int(q), int(r), int(s)]

                if sign == "+":
                    self.riemann_tensor[i, j, k, l] += Rpqrs
                elif sign == "-":
                    self.riemann_tensor[i, j, k, l] += -Rpqrs
                else:
                    raise ValueError("Unrecognized sign value: {}".format(sign))

            elif len(str_pointer) == 10:
                sign1, p, q, r, s, sign2, t, u, v, w = str_pointer

                Rpqrs = self.riemann_tensor[int(p), int(q), int(r), int(s)]
                Rtuvw = self.riemann_tensor[int(t), int(u), int(v), int(w)]

                if sign1 == "+":
                    self.riemann_tensor[i, j, k, l] += Rpqrs
                elif sign1 == "-":
                    self.riemann_tensor[i, j, k, l] += -Rpqrs
                else:
                    raise ValueError("Unrecognized sign value: {}".format(sign1))

                if sign2 == "+":
                    self.riemann_tensor[i, j, k, l] += Rtuvw
                elif sign2 == "-":
                    self.riemann_tensor[i, j, k, l] += -Rtuvw
                else:
                    raise ValueError("Unrecognized sign value: {}".format(sign2))

            else:
                raise ValueError("Unrecognized symmetry relation: {}".format(str_pointer))


class _PotentialSym(object):

    def __init__(self, n_fields: int, n_params: int, potential: sym.Expr, simplify=True):

        self.nf = n_fields
        self.nf_range = list(range(n_fields))
        self.f_syms = sym.symarray("f", n_fields)  # field symbols
        self.p_syms = sym.symarray("f", n_params) if n_params > 0 else None  # param symbols

        self.simplify = simplify
        self._configure_potential(potential)
        cache_tools.cache_pot(potential, self)

    def _configure_potential(self, potential):

        print("-- Configuring potential")

        self.potential = potential if not self.simplify else sym.simplify(potential)
        self.pd_potential = np.zeros((self.nf), dtype=object)
        self.pdpd_potential = np.zeros((self.nf, self.nf), dtype=object)
        self.pdpdpd_potential = np.zeros((self.nf, self.nf, self.nf), dtype=object)

        print("-- Computing partial derivatives")

        for ii, f_ii in enumerate(self.f_syms):

            pdv = sym.diff(potential, f_ii)

            if self.simplify:
                pdv = sym.simplify(pdv)

            self.pd_potential[ii] = pdv

            for jj, f_jj in enumerate(self.f_syms):

                if jj >= ii:
                    pdpdv = sym.diff(self.pd_potential[ii], f_jj)

                    if self.simplify:
                        pdpdv = sym.simplify(pdpdv)

                    self.pdpd_potential[ii, jj] = pdpdv

                else:
                    self.pdpd_potential[ii, jj] = self.pdpd_potential[jj, ii]

                for kk, f_kk in enumerate(self.f_syms):

                    if kk >= jj:

                        pdpdpdv = sym.diff(self.pdpd_potential[ii, jj], f_kk)

                        if self.simplify:
                            pdpdpdv = sym.simplify(pdpdpdv)

                        self.pdpdpd_potential[ii, jj, kk] = pdpdpdv

                    else:
                        self.pdpdpd_potential[ii, jj, kk] = self.pdpdpd_potential[ii, kk, jj]


nF = 2  # number of fields needed to define the double quadratic potential
nP = 3  # number of parameters needed to define the double quadtartic potential
f = sym.symarray('f', nF)  # an array representing the nF fields present for this model
p = sym.symarray('p',
                 nP)  # an array representing the nP parameters needed to define this model (that we might wish to change) if we don't
# wish to change them they could be typed explicitly in the potential below

V = 1. / 2. * p[0] ** 2.0 * f[0] ** 2.0 + 1. / 2. * p[1] ** 2.0 * f[
    1] ** 2.0  # this is the potential written in sympy notation
G = sym.Matrix([[p[2] ** 2.0, 0], [0, p[2] ** 2.0 * sym.sin(f[0]) ** 2.0]])

fs = FieldSpaceSym(2, 2, G, recache=False)
ps = PotentialSym(2, 2, V)

# n = 6
# OM = OperationMatrix(n)
#
# print OM

# for i in range(n):
#     for j in range(n):
#         for k in range(n):
#             for l in range(n):
#                 if len(list(set([i, j, k, l])))==4:
#                     print [i, j, k, l], OM[i, j, k, l]

# print OM[0, 3, 1, 3], OM[3, 1, 0, 3]
# print OM[0, 3, 2, 3], OM[3, 2, 0, 3]
# print OM[0, 4, 1, 4], OM[4, 1, 0, 4]
# print OM[0, 4, 2, 4], OM[4, 2, 0, 4]
# print OM[0, 5, 1, 5], OM[5, 1, 0, 5]
# print OM[0, 5, 2, 5], OM[5, 2, 0, 5]
