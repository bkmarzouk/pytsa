import numpy as np

nF = 4

G = np.empty((2*nF*2*nF), dtype="object")

def t(i, j):
    
    if i > 0:
        idx1 = "a"
    else:
        idx1 = "A"
        
    if j > 0:
        idx2 = "b"
    else:
        idx2 = "B"
        
    return "G{}{}[{},{}]".format(idx1, idx2, i, j)
        
# for i in range(2 * nF):
#     for j in range(2 * nF):
#         if i < nF:
#             ii = -i - 1
#         else:
#             ii = i - (nF - 1)
#         if j < nF:
#             jj = -j - 1
#         else:
#             jj = j - (nF - 1)
#
#         G[(2 * nF) * i + j] = t(ii, jj)
#
# print G.reshape((2*nF, 2*nF))

# stupid = np.empty((nF**3), dtype=object)
#
# for i in range(nF):
#     for j in range(nF):
#         for k in range(nF):
#             # print i, j*nF, k*nF*nF, "d{}d{}d{}V".format(i, j, k)
#             stupid[i+j*nF+k*nF*nF] = "d{}d{}d{}V".format(i, j, k)
#
# s_3d = stupid.reshape((nF, nF, nF))
#
# sense = np.empty((nF, nF, nF), dtype=object)
#
# for i in range(nF):
#     for j in range(nF):
#         for k in range(nF):
#             sense[i][j][k] = "d{}d{}d{}V".format(i, j, k)
#             print sense[i, j, k] == sense[i][j][k]
#
# print stupid == sense.T.flatten()

Gamma_array_pyt = np.empty((2 * nF * 2 * nF * 2 * nF), dtype=object)

# populate connexion matrix
for i in range(2 * nF):
    for j in range(2 * nF):
        for k in range(2 * nF):
            if i < nF:
                ii = -i - 1
            else:
                ii = i - (nF - 1)
            if j < nF:
                jj = -j - 1
            else:
                jj = j - (nF - 1)
            if k < nF:
                kk = -k - 1
            else:
                kk = k - (nF - 1)
            
            if kk < 0 or jj < 0 or ii > 0:
                Gamma_array_pyt[(2 * nF) * (2 * nF) * i + (2 * nF) * j + k] = 0
            else:
                Gamma_array_pyt[(2 * nF) * (2 * nF) * i + (2 * nF) * j + k] = "C[{},{},{}]".format(
                    ii,jj,kk, i, j, k)
                
            print [i, j, k], Gamma_array_pyt[(2 * nF) * (2 * nF) * i + (2 * nF) * j + k]
                    
# print Gamma_array_pyt.reshape((2*nF, 2*nF, 2*nF))