import numpy as np
from scipy.special import binom

def OperationMatrix(N):

    """ Computes and N X N X N X N matrix corresponding to the components of the Riemann tensor in N dimensions """

    # Initialize N^4 array
    OM = np.zeros([N, N, N, N], dtype=object)

    # Initialize counters for non-zero independent components
    idp_ijij = 0 # R_ijij
    idp_ijki = 0 # R_ijki
    idp_ijkl = 0 # R_ijkl

    # Track components which are reducible via the Bianchi cyclic identity
    bterms = []

    # Iterate over N coords
    for i in range(N):

        jrange = range(N); jrange.remove(i)

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

            krange = range(N); krange.remove(i); krange.remove(j)

            # Iterate over (N-2) coords
            for k in krange:

                # Choose independent solutions to have k < j
                if k < j:

                    OM[i, j, k, i] = "calc"
                    idp_ijki += 1

                    # Prescribe dependent components in terms of fundamental component & +/-
                    subs  = [i, j, k, i]
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

                lrange = range(N); lrange.remove(i); lrange.remove(j); lrange.remove(k)

                # iterate over (N-3) coords
                for l in lrange:

                    # Choose irreducible representation to be those with i < j & k < l in first & 2nd index pairs.
                    # To avoid double counting under exchange symmetry, we items with lowest index in first index pair
                    if i < j and k < l and min(i, j) <= min(k, l):

                        OM[i, j, k, l] = "calc"
                        idp_ijkl+=1

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

                    else: pass

    c = 0

    # Copy terms which could be related via the Bianchi identity
    bterms_copy = bterms

    def exch(i, j, k, l): return [k, l, i, j]

    for item in bterms_copy:
        i, j, k, l = item

        # Cyclic permute item in bterms
        b1 = [i, k, l, j]
        b2 = [i, l, j, k]

        # Cyclic permute item in bterms and skew indices
        b1s = [i, k, j, l]
        b2s = [i, l, k, j]

        # We have then the following possible representations for the Riemann tensor
        con1 = (b1  in bterms or exch(*b1)  in bterms) and (b2  in bterms or exch(*b2)  in bterms)
        con2 = (b1s in bterms or exch(*b1s) in bterms) and (b2s in bterms or exch(*b2s) in bterms)
        con3 = (b1  in bterms or exch(*b1)  in bterms) and (b2s in bterms or exch(*b2s) in bterms)
        con4 = (b1s in bterms or exch(*b1s) in bterms) and (b2  in bterms or exch(*b2)  in bterms)

        if   con1:
            bterms.remove(item); idp_ijkl+=-1
            OM[i, j, k, l] = "-{}{}{}{}".format(*b1) + "-{}{}{}{}".format(*b2)

        elif con2:
            bterms.remove(item); idp_ijkl+=-1
            OM[i, j, k, l] = "+{}{}{}{}".format(*b1s) + "+{}{}{}{}".format(*b2s)

        elif con3:
            bterms.remove(item); idp_ijkl+=-1
            OM[i, j, k, l] = "-{}{}{}{}".format(*b1) + "+{}{}{}{}".format(*b2s)

        elif con4:
            bterms.remove(item); idp_ijkl+=-1
            OM[i, j, k, l] = "+{}{}{}{}".format(*b1s) + "-{}{}{}{}".format(*b2)
        else: pass

    # Enforce analytic consistency check on the number of independent items
    assert 2  *idp_ijij == N *(N - 1), str(idp_ijij) + "/" + str(N * (N - 1)/2)
    assert 2  *idp_ijki == N *(N - 1) * (N - 2), str(idp_ijki) + "/" + str(N * (N - 1) * (N - 2)/2)
    assert 12 *idp_ijkl == N *(N - 1) * (N - 2) * (N - 3), str(idp_ijkl) + "/" + str(N * (N - 1) * (N - 2) * (N - 3)/12)

    print "-- Independent components: {}".format(idp_ijij+idp_ijki+idp_ijkl)

    # Return matrix
    return OM

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

