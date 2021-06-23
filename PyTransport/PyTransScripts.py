# This file is part of _PyTransport.

# _PyTransport is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# _PyTransport is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with _PyTransport.  If not, see <http://www.gnu.org/licenses/>.


# python code contains some useful scripts to use with the compiled PyTrans module.

import gc
import os
import pickle as pk
import sys
import timeit

import numpy as np
import sympy as sym
from PyTransport.gravtools_pyt import curvatureObject
from scipy import interpolate
from scipy.interpolate import UnivariateSpline
from scipy.optimize import minimize_scalar


# this script finds initial conditions at least NBMassless e-folds before the massless point
# back must finely sample the backgroudn evolution for the initial conditions to be close to exactly NBMassless before

def unPackAlp(threePtOut, MTE):
    nF = MTE.nF()
    if np.size(threePtOut[0, :]) != 1 + 4 + 2 * nF + 6 * 2 * nF * 2 * nF + 2 * nF * 2 * nF * 2 * nF:
        print("\n\n\n\n warning array you asked to unpack is not of correct dimension \n\n\n\n")
        return np.nan, np.nan, np.nan, np.nan, np.nan
    zetaMs = threePtOut[:, 1:5]
    sig1R = threePtOut[:, 1 + 4 + 2 * nF:1 + 4 + 2 * nF + 2 * nF * 2 * nF]
    sig2R = threePtOut[:, 1 + 4 + 2 * nF + 2 * nF * 2 * nF:1 + 4 + 2 * nF + 2 * 2 * nF * 2 * nF]
    sig3R = threePtOut[:, 1 + 4 + 2 * nF + 2 * 2 * nF * 2 * nF:1 + 4 + 2 * nF + 3 * 2 * nF * 2 * nF]
    alpha = threePtOut[:, 1 + 4 + 2 * nF + 6 * 2 * nF * 2 * nF:]

    return zetaMs, np.reshape(sig1R, (np.size(threePtOut[:, 0]), 2 * nF, 2 * nF)), np.reshape(sig2R, (
        np.size(threePtOut[:, 0]), 2 * nF, 2 * nF)), np.reshape(sig3R,
                                                                (
                                                                    np.size(threePtOut[:, 0]), 2 * nF,
                                                                    2 * nF)), np.reshape(
        alpha, (np.size(threePtOut[:, 0]), 2 * nF, 2 * nF, 2 * nF))


def unPackSig(twoPtOut, MTE):
    nF = MTE.nF()
    if np.size(twoPtOut[0, :]) != 1 + 1 + 2 * nF + 2 * nF * 2 * nF:
        print("\n\n\n\n warning array you asked to unpack is not of correct dimension \n\n\n\n")
        return np.nan, np.nan
    zeta = twoPtOut[:, 1]
    sig = twoPtOut[:, 1 + 1 + 2 * nF:1 + 1 + 2 * nF + 2 * nF * 2 * nF]

    return zeta, np.reshape(sig, (np.size(twoPtOut[:, 0]), 2 * nF, 2 * nF))


def arr_eval_spl_N(N, arr, k=3):
    """
    Estimate time ordered array of values at exact time N, via interpolating through
    array data on the interval [N-0.5, N+0.5]

    :param N: evaluation time (efolds, N)
    :param arr: time ordered array of evolving quantities. Zeroth column should contain 'N' values.
    Zeroth should contain quantities at the earliest time
    :param k: order (knots) used to fit spline
    :return: array quantities at the time step N
    """
    if N in arr.T[0]:
        return arr[np.where(N == arr.T[0])]

    assert N > 1, N
    assert N < arr[-1][0] - 1, N

    out = np.zeros(arr[0].shape)

    for step in arr:
        N_step = step[0]

        if N_step > N + 0.5:
            break

        if N_step > N - 0.5:
            out = np.vstack([out, step])

    out = out[1:]

    Nspl = out.T[0]

    fspl = [UnivariateSpline(Nspl, f, k=k) for f in out.T[1:]]

    return np.array([N] + [spl(N) for spl in fspl])


def arr_eval_spl_val(val, arr, col_idx, k=3, spl_steps=5):
    """
    Estimate time ordered array of values at exact time N, via interpolating through
    array data on the interval [N-0.5, N+0.5]

    :param N: evaluation time (efolds, N)
    :param arr: time ordered array of evolving quantities. Zeroth column should contain 'N' values.
    Zeroth should contain quantities at the earliest time
    :param k: order (knots) used to fit spline
    :return: array quantities at the time step N
    """
    if val in arr.T[col_idx]:
        return arr[np.where(val == arr.T[col_idx])]

    idx = None

    for idx, step in enumerate(arr):

        if step[col_idx] > val:
            break

    if idx is None:
        raise ValueError("Unable to locate appropriate value range within input array")

    out = arr[np.max([0, idx - spl_steps // 2]): np.min([len(arr), idx + spl_steps // 2])]

    if len(out) < 2:
        raise ValueError("Unable to locate appropriate value range within input array")

    if k >= len(out):
        k = len(out) - 1

    val_spl = out.T[col_idx]

    spls = [UnivariateSpline(val_spl, col, k=k) for col in out.T]

    return np.array([s(val) for s in spls])


def arr_eval_nearest_N(N, arr):
    """
    Estimate time ordered array of values at approximate time N, via locating the nearest
    available time step

    :param N: evaluation time (efolds, N)
    :param arr: time ordered array of evolving quantities. Zeroth column should contain 'N' values.
    Zeroth should contain quantities at the earliest time
    :return: array quantities at the time step N
    """
    if N in arr.T[0]:
        return arr[np.where(N == arr.T[0])]

    assert N > 1, N
    assert N < arr[-1][0] - 1, N

    for ii, step in enumerate(arr):

        if step[0] > N:
            jj = ii - 1
            dNii = step[0] - N
            dNjj = N - arr[jj][0]
            return arr[jj if dNjj < dNii else ii]


def rescaleBack(bg, Nr=70.0, reposition=True):
    """
    
    Repositions background evolution to begin Nr efolds before the end of inflation.
    Helps avoid momenta growing exponentially large.
    
    If Nr exceeds the length of inflation, the original background is returned.
    
    Finally, we also ensure that the first 10 efolds of inflation are finely sampled enough s.t.
    accurate estimates of ICsBM can be made
    
    """

    # Get final efold of background evolution
    Nend = bg[-1][0]
    Nstart = bg[0][0]

    if reposition:

        # If background evolution is sufficiently long
        if (Nend - Nstart) > Nr:

            # Copy final row of evolution
            out = np.copy(bg[-1])

            # Iterate over back steps from last to first
            for row in np.flipud(bg[:-1]):

                # Compute difference between current efold and end of inflation
                dN = Nend - row[0]

                # Stack ontop of output array
                out = np.vstack((row, out))

                # If dN has reached the target reposition
                if dN > Nr:
                    # Stop background repositioning, and redefine N efolution to start at zero from this point
                    out[:, 0] -= out[0][0]

                    break
        else:
            out = bg

    else:
        out = bg

    # Cross check final field values haven't changed
    assert np.all(out[-1][1:] == bg[-1][1:]), [out[-1][1:], bg[-1][1:]]

    # Finally, we check that the background hasn't been too finely sampled, s.t back[ii] == back[ii+1]
    killRows = np.array([
        ii for ii in range(len(out) - 1) if np.all(out[ii] == out[ii + 1])
    ])

    killRows = killRows[::-1]

    for ii in killRows:
        out = np.delete(out, ii, axis=0)

    return out


def ICsBM(NBMassless, k, back, params, MTE, fit="spl"):

    assert fit in ["spl", "nearest"], fit

    masses_array = evolveMasses(back, params, MTE)

    m_array = np.zeros((len(back), 3))

    for idx, (m_eigs, step) in enumerate(zip(masses_array, back)):
        m_array[idx][0] = m_eigs[0]  # N
        m_array[idx][1] = np.max(m_eigs[1:])  # M^2 / H^2
        m_array[idx][2] = k ** 2 / MTE.H(step[1:], params) ** 2 / np.exp(2 * step[0])  # k^2 / (aH)^2

    sign_init = np.sign(m_array[0][1] - m_array[0][2])

    idx = None
    N = None

    for idx, nmk in enumerate(m_array):
        N, M, K = nmk

        if np.sign(M - K) != sign_init:
            break

    if idx is None:
        raise ValueError("Unable to locate massless transition")

    N_start_subhorizon = N - NBMassless

    if N_start_subhorizon < 0:
        raise ValueError("Unable to locate massless initial condition")

    if fit == "nearest":
        out = arr_eval_nearest_N(N_start_subhorizon, back)

    out = arr_eval_spl_N(N_start_subhorizon, back)

    return out[0], out[1:]


def ICsBE(NBExit, k, back, params, MTE):
    nF = np.size(back[0, 1:]) // 2
    kvaH = -1.;
    jj = 0
    while (kvaH < 0.0 and jj < np.size(back[:, 0]) - 1):
        H = MTE.H(back[jj, 1:1 + 2 * nF], params)
        kvaH = -k + np.exp(back[jj, 0]) * H
        jj = jj + 1
    if jj == np.size(back[:, 0]):
        print("\n\n\n\n warning exit condition not found \n\n\n\n")
        return np.nan, np.nan
    NExit = back[jj - 2, 0]
    ll = 0
    Ncond = -1
    while (Ncond < 0 and ll < np.size(back[:, 0]) - 1):
        Ncond = back[ll, 0] - (NExit - NBExit)
        ll = ll + 1

    if ll == np.size(back[:, 0]) or (NExit - back[0, 0]) < NBExit:
        print("\n\n\n\n warning initial condition not found \n\n\n\n")
        return np.nan, np.nan

    NexitMinus = back[ll - 2, 0]
    backExitMinus = back[ll - 2, 1:]

    return NexitMinus, backExitMinus


# find the earliest condition between the massless one and the horizon exit one
def ICs(NB, k, back, params, MTE):
    NBEs, fieldBE = ICsBE(NB, k, back, params, MTE)
    NBMs, fieldBM = ICsBM(NB, k, back, params, MTE)

    if (NBEs < NBMs):
        return NBEs, fieldBE
    return NBMs, fieldBM


# calculates the power spectrum at each element in kA at the end of the background evolution (back)
def pSpectra(kA, back, params, NB, tols, MTE):
    zzOut = np.array([])
    times = np.array([])
    num = np.size(kA)

    for ii in range(0, num):
        print("\n \n \n performing " + str(ii + 1) + " of " + str(num) + "\n \n \n")
        k = kA[ii]
        Nstart, backExitMinus = ICs(NB, k, back, params, MTE)
        start_time = timeit.default_timer()

        if Nstart == np.nan:
            twoPt = numpy.empty((2, 2))
            twoPt[:] = np.nan
        else:
            t = np.linspace(Nstart, back[-1, 0], 10)
            # run solver for this triangle
            twoPt = MTE.sigEvolve(t, k, backExitMinus, params, tols,
                                  True, tmax_2pf, flagReturn)  # all data from three point run goes into threePt array
        zzOut = np.append(zzOut, twoPt[-1, 1])
        times = np.append(times, timeit.default_timer() - start_time)

    return zzOut, times


# calculates the power spectrum at each element in kA at the end of the background evolution (back) in a manner
# suitable to be called over many processes
def pSpecMpi(kA, back, params, NB, tols, MTE):
    from mpi4py import MPI
    comm = MPI.COMM_WORLD

    rank = comm.Get_rank()
    size = comm.Get_size()
    points = np.size(kA)
    num = points / size;

    if float(points) / size != float(points // size):
        if rank == 0:
            print("\n \n \n warning! number of points is divisable by number of processes, exiting \n \n \n ")
        return (np.empty, np.empty)

    kOutL = kA[rank * num:rank * num + num]

    zzL, timesL = pSpectra(kOutL, back, params, NB, tols, MTE)

    if rank != 0:
        comm.Send(zzL, dest=0)
        comm.Send(timesL, dest=0)

    if rank == 0:
        zzzOut = np.array([])
        zzOut = np.array([])
        timesOut = np.array([])
        zzOut = np.append(zzOut, zzL)
        timesOut = np.append(timesOut, timesL)

        for jj in range(1, size):
            comm.Recv(zzL, source=jj)
            comm.Recv(timesL, source=jj)
            zzOut = np.append(zzOut, zzL)
            timesOut = np.append(timesOut, timesL)

        return (zzOut, timesOut)
    else:
        return (np.empty, np.empty)


# calculates the power spectrum and bisecpturm in equilateral configuration at each element in kA at the end of the
# background evolution (back)
def eqSpectra(kA, back, params, NB, tols, MTE):
    zzzOut = np.array([])
    zzOut = np.array([])
    times = np.array([])
    num = np.size(kA)

    for ii in range(0, num):
        print("\n \n \n performing " + str(ii + 1) + " of " + str(num) + "\n \n \n")
        k = kA[ii]
        Nstart, backExitMinus = ICs(NB, k, back, params, MTE)
        t = np.linspace(Nstart, back[-1, 0], 10)
        k1 = k;
        k2 = k;
        k3 = k;
        # run solver for this triangle
        start_time = timeit.default_timer()

        if Nstart == np.nan:
            nF = MTE.nF();
            threePt = numpy.empty((2, 5))
            threePt[:] = np.nan
        else:
            threePt = MTE.alphaEvolve(t, k1, k2, k3, backExitMinus, params, tols,
                                      True)  # all data from three point run goes into threePt array
        zzOut = np.append(zzOut, threePt[-1, 1])
        zzzOut = np.append(zzzOut, threePt[-1, 4])
        times = np.append(times, timeit.default_timer() - start_time)
    return zzOut, zzzOut, times


# calculates the power spectrum and bisecpturm in equilateral configuration at each element in kA at the end of the
# background evolution (back) in a manner suitable to be run accoss many processes
def eqSpecMpi(kA, back, params, NB, tols, MTE):
    from mpi4py import MPI
    comm = MPI.COMM_WORLD

    rank = comm.Get_rank()
    size = comm.Get_size()
    points = np.size(kA)
    num = points / size

    if float(points) / size != float(points // size):
        if rank == 0:
            print("\n \n \n warning! number of points is divisable by number of processes, exiting \n \n \n ")
        return (np.empty, np.empty, np.empty)

    kOutL = kA[rank * num:rank * num + num]

    zzL, zzzL, timesL = eqSpectra(kOutL, back, params, NB, tols, MTE)

    if rank != 0:
        comm.Send(zzzL, dest=0)
        comm.Send(zzL, dest=0)
        comm.Send(timesL, dest=0)

    if rank == 0:
        zzzOut = np.array([])
        zzOut = np.array([])
        timesOut = np.array([])
        zzzOut = np.append(zzzOut, zzzL)
        zzOut = np.append(zzOut, zzL)
        timesOut = np.append(timesOut, timesL)

        for jj in range(1, size):
            comm.Recv(zzzL, source=jj)
            comm.Recv(zzL, source=jj)
            comm.Recv(timesL, source=jj)

            zzzOut = np.append(zzzOut, zzzL)
            zzOut = np.append(zzOut, zzL)
            timesOut = np.append(timesOut, timesL)

        return (zzOut, zzzOut, timesOut)
    else:
        return (np.empty, np.empty, np.empty)


# calcualtes the bispectrum in the alpha beta notation for a given kt at every value of the alphaIn and betIn arrays.
# The bispectrum is given an nsnaps times always incuding the final time of the evolution (back)
def alpBetSpectra(kt, alphaIn, betaIn, back, params, NB, nsnaps, tols, MTE):
    Hin = np.zeros(np.size(back[:, 0]))
    for jj in range(0, np.size(back[:, 0])):
        Hin[jj] = MTE.H(back[jj, 1:], params)
    aH = np.exp(back[:, 0]) * Hin
    positions = np.argsort(aH)
    Nexit = interpolate.splev(kt / 3., interpolate.splrep(aH[positions], back[positions, 0], s=1e-15), der=0)
    Nend = back[-1, 0]
    snaps = np.linspace(Nexit - (NB - .1), Nend, nsnaps)
    if (nsnaps == 1 or nsnaps == 0):
        snaps = np.array([Nend])
        nsnaps = 1

    biAOut = np.zeros([np.size(alphaIn), np.size(betaIn), np.size(snaps)])
    zz1 = np.zeros([np.size(alphaIn), np.size(betaIn), np.size(snaps)])
    zz2 = np.zeros([np.size(alphaIn), np.size(betaIn), np.size(snaps)])
    zz3 = np.zeros([np.size(alphaIn), np.size(betaIn), np.size(snaps)])
    times = np.zeros([np.size(alphaIn), np.size(betaIn)])

    for l in range(0, np.size(alphaIn)):
        alpha = alphaIn[l]
        for j in range(0, np.size(betaIn)):
            print("\n \n \n performing " + str(l + 1) + " " + str(j + 1) + " of " + str(np.size(alphaIn)) + " " + str(
                np.size(betaIn)) + "\n \n \n")
            timebefore = timeit.default_timer()
            beta = betaIn[j]
            k1 = kt / 2. - beta * kt / 2.
            k2 = kt / 4. * (1. + alpha + beta)
            k3 = kt / 4. * (1. - alpha + beta)
            if alpha > -(1 - beta) and alpha < 1 - beta:
                kM = min(k1, k2, k3)
                Nstart, backExitMinus = ICs(NB, kM, back, params, MTE)
                # run solver for this triangle
                t = np.concatenate((np.array([Nstart]), snaps))
                if Nstart == np.nan:
                    threePt = numpy.empty((2, 5))
                    threePt[:] = np.nan
                else:
                    threePt = MTE.alphaEvolve(t, k1, k2, k3, backExitMinus, params, tols, True)
                zzz = threePt[:, :5]

                for ii in range(1, nsnaps + 1):
                    biAOut[l, j, ii - 1] = zzz[ii, 4]
                    zz1[l, j, ii - 1] = zzz[ii, 1]
                    zz2[l, j, ii - 1] = zzz[ii, 2]
                    zz3[l, j, ii - 1] = zzz[ii, 3]
            else:
                for ii in range(0, nsnaps + 1):
                    biAOut[l, j, ii - 1] = np.nan
                    zz1[l, j, ii - 1] = np.nan
                    zz2[l, j, ii - 1] = np.nan
                    zz3[l, j, ii - 1] = np.nan

            times[l, j] = timeit.default_timer() - timebefore
    return (biAOut, zz1, zz2, zz3, times, snaps)


# performs the same task as alpBetSpectra but in a manner suitable to be spread across many processes
def alpBetSpecMpi(kt, alpha, beta, back, params, NB, nsnaps, tols, MTE):
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    side = np.size(alpha)
    Nbefore = NB
    rank = comm.Get_rank()
    size = comm.Get_size()
    num = side / size;

    if float(side) / size != float(side // size):
        if rank == 0:
            print("\n \n \n warning! number size of alpha must be divisable by number of processes, exiting \n \n \n ")
        return (np.empty, np.empty, np.empty, np.empty, np.empty)
    else:

        alphaL = alpha[rank * num:rank * num + num]

        BzL, Pz1L, Pz2L, Pz3L, timesL, snaps = alpBetSpectra(kt, alphaL, beta, back, params, Nbefore, nsnaps, tols, MTE)

        if rank != 0:
            comm.send(Pz1L, dest=0)
            comm.send(Pz2L, dest=0)
            comm.send(Pz3L, dest=0)
            comm.send(BzL, dest=0)
            comm.send(timesL, dest=0)

        if rank == 0:
            Bztot = np.zeros([np.size(alpha), np.size(beta), np.size(BzL[0, 0, :])])
            Pz1tot = np.zeros([np.size(alpha), np.size(beta), np.size(BzL[0, 0, :])])
            Pz2tot = np.zeros([np.size(alpha), np.size(beta), np.size(BzL[0, 0, :])])
            Pz3tot = np.zeros([np.size(alpha), np.size(beta), np.size(BzL[0, 0, :])])
            timestot = np.zeros([np.size(alpha), np.size(beta)])
            Bztot[0:num, :, :] = BzL
            Pz1tot[0:num, :, :] = Pz1L
            Pz2tot[0:num, :, :] = Pz2L
            Pz3tot[0:num, :, :] = Pz3L
            timestot[0:num, :] = timesL

            for jj in range(1, size):
                Pz1L = comm.recv(source=jj)
                Pz2L = comm.recv(source=jj)
                Pz3L = comm.recv(source=jj)
                BzL = comm.recv(source=jj)
                timesL = comm.recv(source=jj)

                Bztot[jj * num:jj * num + num, :, :] = BzL
                Pz1tot[jj * num:jj * num + num, :, :] = Pz1L
                Pz2tot[jj * num:jj * num + num, :, :] = Pz2L
                Pz3tot[jj * num:jj * num + num, :, :] = Pz3L
                timestot[jj * num:jj * num + num, :] = timesL
            return (Bztot, Pz1tot, Pz2tot, Pz3tot, timestot, snaps)
        else:
            return (np.empty, np.empty, np.empty, np.empty, np.empty, snaps)


def kexitN(Nexit, back, params, MTE):
    backExitArr = np.vstack([np.zeros(4) for ii in range(2 * MTE.nF())])
    Narr = np.zeros(4)

    count = 0

    for ii in range(len(back) - 1):

        row = back[ii]

        N, fdf = row[0], row[1:]

        Narr = np.roll(Narr, -1)
        Narr[-1] = N

        backExitArr = np.roll(backExitArr, -1, axis=1)
        backExitArr[:, -1] = fdf

        count += 1

        if N > Nexit:

            row = back[ii + 1]

            N, fdf = row[0], row[1:]

            Narr = np.roll(Narr, -1)
            Narr[-1] = N

            backExitArr = np.roll(backExitArr, -1, axis=1)
            backExitArr[:, -1] = fdf

            count += 1

            if count < 4:

                for jj in range(4 - count):
                    row = back[ii + 1 + jj + 1]

                    N, fdf = row[0], row[1:]

                    Narr = np.roll(Narr, -1)
                    Narr[-1] = N

                    backExitArr = np.roll(backExitArr, -1, axis=1)
                    backExitArr[:, -1] = fdf

                    count += 1

            break

    if ii == len(back):
        print("** Warning, kExit not found in background **")
        return np.nan

    fdfSpl = np.array([UnivariateSpline(Narr, row) for row in backExitArr])

    # fdfOut = np.array([spl(0) for spl in fdfSpl])
    fdfOut = np.array([spl(Nexit) for spl in fdfSpl])

    HOut = MTE.H(fdfOut, params)

    # if HOut is None:

    #    HEvo = []
    #    Nevo = []

    #    for row in back:
    #        N, fdf = row[0], row[1:]

    #        H = MTE.H(fdf, params)

    k = np.exp(Nexit + np.log(HOut))

    # print MTE.V(fdfOut[:6], params), HOut, k

    return k


def kexitPhi(PhiExit, n, back, params, MTE):
    nF = np.size(back[0, 1:]) // 2
    backExit = np.zeros(2 * nF)
    positions = np.argsort(back[:, n])
    Nexit = interpolate.splev(PhiExit, interpolate.splrep(back[positions, n], back[positions, 0], s=1e-15), der=0)

    for i in range(1, 2 * nF + 1):
        backExit[i - 1] = interpolate.splev(Nexit, interpolate.splrep(back[:, 0], back[:, i], s=1e-15), der=0)
    k = np.exp(Nexit) * MTE.H(backExit, params);
    return k


def matchKExitN(back, params, MTE, k=0.002):
    const = 55.75 - np.log(k / 0.05)

    kExitArr = np.zeros(4, dtype=float)
    HArr = np.zeros(4, dtype=float)

    count = 0

    # Iterate over background step
    for ii in range(len(back) - 4):

        row = back[ii]
        N, fdf = row[0], row[1:]

        kExit = kexitN(N, back, params, MTE)

        kDelta = kExit - k
        H = MTE.H(fdf, params)

        kExitArr = np.roll(kExitArr, -1)
        HArr = np.roll(HArr, -1)

        kExitArr[-1] = kExit
        HArr[-1] = H

        count += 1

        if kDelta > 0:

            row = back[ii + 1]
            N, fdf = row[0], row[1:]

            kExit = kexitN(N, back, params, MTE)

            H = MTE.H(fdf, params)

            kExitArr = np.roll(kExitArr, -1)
            HArr = np.roll(HArr, -1)

            kExitArr[-1] = kExit
            HArr[-1] = H

            count += 1

            if count < 4:
                for jj in range(4 - count):
                    row = back[ii + 1 + jj + 1]
                    N, fdf = row[0], row[1:]

                    kExit = kexitN(N, back, params, MTE)

                    H = MTE.H(fdf, params)

                    kExitArr = np.roll(kExitArr, -1)
                    HArr = np.roll(HArr, -1)

                    kExitArr[-1] = kExit
                    HArr[-1] = H

                    count += 1

            break

    if ii == len(back) - 4:
        print("\n\n\n\n warning unable to compute pivot time Nk with background \n\n\n\n")

        return np.nan

    try:

        print(kExitArr)
        print(HArr)

        HkSpl = UnivariateSpline(kExitArr, HArr)

        Hk = HkSpl(k)

        Hend = MTE.H(back[-1][1:], params)

        Nk = const + np.log(243.5 * 3 ** 0.25 * np.sqrt(Hk)) + np.log(np.sqrt(Hk / Hend))

    except ValueError:

        print("\n\n\n\n warning error fitting spline for pivot exit scale: {}\n\n\n\n".format(k))

        return np.nan

    return Nk


def evolveMasses(back, params, MTE, scale_eigs=False, hess_approx=False, covariant=False):
    """ Computes the mass matrix along a given background evolution, M^{I}_{J} """

    # Build empty array to populate with efold number + eigenvalues of mass-matrix at time step
    eigs = np.empty((len(back), 1 + MTE.nF()))

    for idx, step in enumerate(back):

        # Unpack efolding and fields as np arrays
        N, fieldsdotfields = step[0], step[1:1 + 2 * MTE.nF()]

        # Compute mass matrix
        Mij = MTE.massMatrix(fieldsdotfields, params, hess_approx, covariant)

        # Compute and sort eigenvalues
        masses = np.linalg.eigvals(Mij)
        masses = np.sort(masses)

        # If scale eigs m^2 -> m
        if scale_eigs: masses = np.sign(masses) * np.sqrt(np.abs(masses))

        # Assign row values for mass-matrix
        eigs[idx] = np.concatenate((np.array([N]), masses))

    return eigs


def spectralIndex(back, pvals, Nexit, tols, subevo, MTE, returnRunning=True,
                  errorReturn=False, tmax=None, useMatchEq=True, return_logarr=False):
    """

    Simple based method to compute the spectral index of a mode that exits the horizon at a time Nend - Nexit

    If kPivot is None, pivot scale is treated as k(Nexit), else kPivot

    If returnRunning is True, returns the first running of the spectral index

    """

    # Get end of background evolution
    Nend = back.T[0][-1]

    if useMatchEq:
        Nexit = matchKExitN(back, pvals, MTE)

    if np.isnan(Nexit):

        print("\n\n\n\n warning error computin Nexit from matchin \n\n\n\n")

        if errorReturn:
            return ValueError, "ics"

        raise ValueError(Nexit)

    # We will compute the 2pf based on momenta that horizon exit 1.2-efolds above & below k pivot,
    # defined as the mode that exits the horizon at a time Nexit
    Npivot = Nend - Nexit
    Nscatter = 0.6
    Nsteps = 2
    DeltaN = [Npivot - Nscatter + ii * Nscatter / float(Nsteps) for ii in range(Nsteps * 2 + 1)]

    # Build numpy array of momenta
    kVals = np.array([])

    # Iterate over efolds about exit
    for NN in DeltaN:

        # Get horizon exit mode
        k = kexitN(NN, back, pvals, MTE)

        # if momenta has become infinite, or subject to numerical error
        if np.isinf(k) or np.isnan(k):

            if errorReturn:
                return ValueError, "k"

            raise ValueError(k)

        # Append to momenta values for spline calculation
        kVals = np.append(kVals, k)

    # Check momenta are strictly monotonically increasing
    if not all(kVals[i] < kVals[i + 1] for i in range(len(kVals) - 1)):

        if errorReturn:
            return ValueError, "k"

        raise ValueError(kVals)

    kExit = kVals[Nsteps]

    # Get pivot scale, midpoint in spline
    kPivot = kVals[Nsteps]

    if np.isnan(kPivot) or np.isinf(kPivot):
        print("\n\n\n\n warning kPivot parameter must be float: {} \n\n\n\n".format(kPivot))

        if errorReturn:
            return ValueError, "k"

        raise ValueError(kVals)

    # Build power spectra values
    pZetaVals = np.array([])

    # Iterate over k-scales
    for k in kVals:

        # Build initial conditions for mode based on massless condition
        Nstart, ICs = ICsBM(subevo, k, back, pvals, MTE)

        # Check ICs are valid
        if type(ICs) != np.ndarray and np.isnan(ICs):

            if errorReturn:
                return ValueError, "ics"

            raise ValueError(ICs)

        # If no time out is given, assume no kill time
        if tmax is None:
            tmax = -1

        twoPf = MTE.sigEvolve(np.array([Nstart, Nend]), k, ICs, pvals, tols, False, tmax, True)

        # Handle failed calculation
        if type(twoPf) is tuple:

            if errorReturn:
                return ValueError, "2pf", twoPf

            raise ValueError(twoPf)

        # Append power spedctrum at end of background
        pZetaVals = np.append(pZetaVals, twoPf.T[1][-1])

    #  Build log arrays in k and Pzeta
    arrLogK = np.log(kVals / kPivot)
    arrLogPz = np.log(pZetaVals)

    if return_logarr:
        return arrLogK, arrLogPz

    # Build 4-pt spline from log values
    twoPtSpline = UnivariateSpline(arrLogK, arrLogPz, k=4)

    # Differentiate for ns & running
    nsSpline = twoPtSpline.derivative()
    alphaSpline = nsSpline.derivative()

    # Define functions to map exit scale to observables subject to kPicot
    ns_ = lambda k: nsSpline(np.log(k / kPivot)) + 4.0
    alpha_ = lambda k: alphaSpline(np.log(k / kPivot))

    ns = ns_(kExit)
    alpha = alpha_(kExit)

    if returnRunning:
        return np.array([ns, alpha])

    return ns


def fNL(back, pvals, Nexit, tols, subevo, MTE, alpha=None, beta=None, stdConfig=None,
        errorReturn=False, tmax=None, useMatchEq=True):
    """
    
    Simple wrapper for the reduced bispectrum fNL at the end of inflation, subject to a configuration
    defined by alpha, beta (Fergusson & Shellard), centred about the horizon exit time Nend - Nexit
    
    """

    if alpha is None and beta is None:
        assert stdConfig is not None, "Without supplying alpha and beta parameters, we must pass a standard config. " \
                                      "name: 'squeezed', 'folded', 'equilateral' [or abrv. 'sq', 'fo', 'eq']"

    # If tmax is none, assign no maximum integration time
    if tmax is None:
        tmax = -1

    # Get end of background evolution
    Nend = back.T[0][-1]

    if useMatchEq:
        Nexit = matchKExitN(back, pvals, MTE)

    if np.isnan(Nexit):

        if errorReturn:
            return ValueError, "ics"

        raise ValueError(Nexit)

    Npivot = Nend - Nexit

    kExit = kexitN(Npivot, back, pvals, MTE)

    if kExit == np.inf or kExit == np.nan:

        if errorReturn:
            return ValueError, "k"

        raise ValueError(kExit)

    # Build Fourier triangle via Fergusson Shellard convention;
    _k1 = lambda k, alpha, beta: k * (1. - beta) / 2.
    _k2 = lambda k, alpha, beta: k * (1. + alpha + beta) / 4.
    _k3 = lambda k, alpha, beta: k * (1. - alpha + beta) / 4.

    if stdConfig is not None:

        if stdConfig in ["eq", "equal", "equilateral"]:
            alpha, beta = 0, 1. / 3.

        elif stdConfig in ["fo", "fold", "folded"]:
            alpha, beta = -0.5, 0.5

        elif stdConfig in ["sq", "squeeze", "squeezed"]:
            alpha, beta = 0.9, 0.01

        else:
            assert alpha is not None, "Must supply alpha: {}".format(alpha)
            assert beta is not None, "Must supply beta: {}".format(beta)

    kt = 3 * kExit

    # Compute momentum triangle
    k1 = _k1(kt, alpha, beta)
    k2 = _k2(kt, alpha, beta)
    k3 = _k3(kt, alpha, beta)

    # Check momenta TODO: return error for bad k here
    for k in [k1, k2, k3]:

        # if momenta has become infinite, or subject to numerical error
        if np.isinf(k) or np.isnan(k):

            if errorReturn:
                return ValueError, "k"

            raise ValueError(k)

    # Find larges scale and compute ICs
    kmin = np.min([k1, k2, k3])

    Nstart, ICs = ICsBM(subevo, kmin, back, pvals, MTE)

    if type(ICs) != np.ndarray and np.isnan(ICs):

        if errorReturn:
            return ValueError, "ics"

        raise ValueError(ICs)

    # Compute three-point function up until end of background
    threePt = MTE.alphaEvolve(np.array([Nstart, Nend]), k1, k2, k3, ICs, pvals, tols, True, tmax, True)

    # If flag has been returned when computing 3pf, return flag data
    if type(threePt) is tuple:

        if errorReturn:
            return ValueError, "3pf", threePt

        raise ValueError(threePt)

    # Compute amplitude of fNL at end of inglation
    Pz1, Pz2, Pz3, Bz = [threePt.T[i][-1] for i in range(1, 5)]
    fNL = (5. / 6.) * Bz / (Pz1 * Pz2 + Pz2 * Pz3 + Pz1 * Pz3)

    return fNL


def GetCurvatureObject(MTE):
    """
    Returns the class 'curvature' object associate with the MTE module
    """

    # Find curvature records directory and create file instance to lead curvature class obj.
    dir = os.path.dirname(__file__)
    curv_dir = os.path.join(dir, 'PyTrans', 'CurvatureRecords')
    curv_name = str(MTE).split("'")[1][7:] + ".curvature"
    curv_path = os.path.join(curv_dir, curv_name)
    curv_file = open(curv_path, 'rb')

    with curv_file:
        curv_obj = pk.load(curv_file)

    return curv_obj
