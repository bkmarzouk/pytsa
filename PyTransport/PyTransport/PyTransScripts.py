# This file is part of PyTransport.

# PyTransport is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# PyTransport is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with PyTransport.  If not, see <http://www.gnu.org/licenses/>.


# python code contains some useful scripts to use with the compiled PyTrans module.

import gc
import os
import pickle as pk
import sys
import timeit

import numpy as np
import sympy as sym
from gravtools_pyt import curvatureObject
from scipy import interpolate
from scipy.interpolate import UnivariateSpline
from scipy.optimize import minimize_scalar


# this script finds initial conditions at least NBMassless e-folds before the massless point
# back must finely sample the backgroudn evolution for the initial conditions to be close to exactly NBMassless before

def unPackAlp(threePtOut, MTE):
    nF = MTE.nF()
    if np.size(threePtOut[0, :]) != 1 + 4 + 2 * nF + 6 * 2 * nF * 2 * nF + 2 * nF * 2 * nF * 2 * nF:
        print ("\n\n\n\n warning array you asked to unpack is not of correct dimension \n\n\n\n")
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
        print ("\n\n\n\n warning array you asked to unpack is not of correct dimension \n\n\n\n")
        return np.nan, np.nan
    zeta = twoPtOut[:, 1]
    sig = twoPtOut[:, 1 + 1 + 2 * nF:1 + 1 + 2 * nF + 2 * nF * 2 * nF]
    
    return zeta, np.reshape(sig, (np.size(twoPtOut[:, 0]), 2 * nF, 2 * nF))


def rescaleBack(bg, Nr=70.0):
    """
    
    Repositions background evolution to begin Nr efolds before the end of inflation.
    Helps avoid momenta growing exponentially large.
    
    If Nr exceeds the length of inflation, the original background is returned.
    
    """
    
    # Get final efold of background evolution
    Nend = bg[-1][0]
    Nstart = bg[0][0]
    
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


def ICsBM(NBMassless, k, back, params, MTE, optimize=True):
    """
    
    Computes initial conditions subject to NBMassless efolds before massless condition is realised.
    
    Massless condition: m^2 = (k/a)^2 : m is larges eigenvalue of the mass-squared matrix
    
    optimize runs a spline based estimate of the massless condition
    
    """
    
    # If running optimization, we need to make sure we don't have N_{ii} = N_{ii+1}, s.t.
    # we build splines from strictly monotonic sequences
    evoN = back.T[0]
    
    if optimize:
        killRows = np.array([ii for ii in range(len(back) - 1) if back.T[0][ii] == back.T[0][ii + 1]])
        killRows = killRows[::-1]
        for ii in killRows:
            back = np.delete(back, ii, axis=0)
    
    # Get eigenvalue evolution
    evoEigs = evolveMasses(back, params, MTE)  # mi^2 / H^2
    
    # get fields-dotfields evolution
    evoFields = back[:, 1:]
    Nmassless = None
    
    # Initialize empty arrays for splin-point calculation
    _Nr = np.array([0, 0, 0], dtype=float)
    _mEffr = np.array([0, 0, 0], dtype=float)
    
    for ii in range(len(back) - 1):
        
        # Get N, eigs and field vals. at background step
        N, eigs, fieldsdotfields = evoN[ii], evoEigs[ii], evoFields[ii]
        
        # Compute effective mass
        mEff = np.max(eigs) - (k * np.exp(-N) / MTE.H(fieldsdotfields, params)) ** 2
        
        # Maintin 3x3 array of spline data, permuting, then replacing the last (column) element with the current value
        _mEffr = np.roll(_mEffr, -1)
        _Nr = np.roll(_Nr, -1)
        _mEffr[-1] = mEff
        _Nr[-1] = N
        
        # When mEff turns over to a non-negative value
        if not mEff < 0:
            # Define this to be the massless efold
            Nmassless = N
            
            # Append one further time step of data to the arrays used for spline calculation
            N, eigs, fieldsdotfields = evoN[ii + 1], evoEigs[ii + 1], evoFields[ii + 1]
            mEff = np.max(eigs) - (k * np.exp(-N) / MTE.H(fieldsdotfields, params)) ** 2
            _mEffr = np.append(_mEffr, mEff)
            _Nr = np.append(_Nr, N)
            
            break
    
    if optimize:
        # Build spline of effective mass as a function of N
        m_N_spline = UnivariateSpline(_Nr, _mEffr, k=3)
        
        # We want to find the N that gives an effective mass of zero. This is achived by minimizing
        # the abs val. of the spline function
        minfunc = lambda N: abs(m_N_spline(N))
        Nmassless = minimize_scalar(minfunc).x  # .x evaluates result
    
    Nstart = evoN[0]
    
    if Nmassless is None:
        print ("\n\n\n\n massless condition not found \n\n\n\n")
        return np.nan, np.nan
    
    if Nmassless - NBMassless < Nstart:
        print ("\n\n\n\n Insufficient background to compute NB efolds before massless \n\n\n\n")
        return np.nan, np.nan
    
    # We now find the background data NBMassless efolds before the massless condition.
    
    icsN = None
    icsFields = None
    
    # Initialize arrays for spline based calculation
    _icsFieldsr = np.vstack((np.zeros(3, dtype=float) for ii in range(2 * MTE.nF())))
    _icsNr = np.zeros(3, dtype=float)
    _Neffr = np.zeros(3, dtype=float)
    
    # Iterate over background steps
    for ii in range(len(evoN) - 1):
        
        # Unpack N and field-dot-fields
        N, fieldsdotfields = evoN[ii], evoFields[ii]
        
        # When zero, we have our ICs
        Neff = N - (Nmassless - NBMassless)
        
        # Permute field array and update right most col with fielddotfield vals
        _icsFieldsr = np.roll(_icsFieldsr, -1, axis=1)
        _icsFieldsr[:, -1] = np.array([fdf for fdf in fieldsdotfields])
        
        # Permute ics N array and update with current efold
        _icsNr = np.roll(_icsNr, -1)
        _icsNr[-1] = N
        
        # Permute effective N array and update with current Neff
        _Neffr = np.roll(_Neffr, -1)
        _Neffr[-1] = Neff
        
        # when N - Nmassless - Nstart >= NB 
        if not Neff < 0:
            icsN = N
            icsFields = fieldsdotfields
            
            _icsFieldsr = np.hstack((_icsFieldsr, np.array([[fdf] for fdf in evoFields[ii + 1]])))
            _icsNr = np.append(_icsNr, evoN[ii + 1])
            _Neffr = np.append(_Neffr, (evoN[ii + 1] - Nstart) - (Nmassless - NBMassless))
            
            break
    
    if optimize:
        
        Neff_Nic_spline = UnivariateSpline(_Neffr, _icsNr, k=3)
        Neff_Fic_splines = [UnivariateSpline(_Neffr, row) for row in _icsFieldsr]
        
        icsN = Neff_Nic_spline(0)
        icsFields = np.array([spl(0) for spl in Neff_Fic_splines])
    
    if icsN is None:
        print ("\n\n\n\n massless condition not found (BUG!) \n\n\n\n")
        return np.nan, np.nan
    
    return icsN, icsFields


def ICsBE(NBExit, k, back, params, MTE):
    nF = np.size(back[0, 1:]) // 2
    kvaH = -1.;
    jj = 0
    while (kvaH < 0.0 and jj < np.size(back[:, 0]) - 1):
        H = MTE.H(back[jj, 1:1 + 2 * nF], params)
        kvaH = -k + np.exp(back[jj, 0]) * H
        jj = jj + 1
    if jj == np.size(back[:, 0]):
        print ("\n\n\n\n warning exit condition not found \n\n\n\n")
        return np.nan, np.nan
    NExit = back[jj - 2, 0]
    ll = 0
    Ncond = -1
    while (Ncond < 0 and ll < np.size(back[:, 0]) - 1):
        Ncond = back[ll, 0] - (NExit - NBExit)
        ll = ll + 1
    
    if ll == np.size(back[:, 0]) or (NExit - back[0, 0]) < NBExit:
        print  ("\n\n\n\n warning initial condition not found \n\n\n\n")
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
        print ("\n \n \n performing " + str(ii + 1) + " of " + str(num) + "\n \n \n")
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
            print ("\n \n \n warning! number of points is divisable by number of processes, exiting \n \n \n ")
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
        print ("\n \n \n performing " + str(ii + 1) + " of " + str(num) + "\n \n \n")
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
            print ("\n \n \n warning! number of points is divisable by number of processes, exiting \n \n \n ")
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
            print ("\n \n \n performing " + str(l + 1) + " " + str(j + 1) + " of " + str(np.size(alphaIn)) + " " + str(
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
            print ("\n \n \n warning! number size of alpha must be divisable by number of processes, exiting \n \n \n ")
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
    nF = np.size(back[0, 1:]) // 2
    backExit = np.zeros(2 * nF)
    
    for i in range(1, 2 * nF + 1):
        backExit[i - 1] = interpolate.splev(Nexit, interpolate.splrep(back[:, 0], back[:, i], s=1e-15), der=0)
    k = np.exp(Nexit) * MTE.H(backExit, params);
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


def spectralIndex(back, pvals, Nexit, tols, subevo, MTE, kPivot=None, returnRunning=True, errorReturn=False, tmax=None):
    """

    Simple spline based method to compute the spectral index of a mode that exits the horizon at a time Nend - Nexit

    If kPivot is None, pivot scale is treated as k(Nexit), else kPivot

    If returnRunning is True, returns the first running of the spectral index

    """
    
    # Get end of background evolution
    Nend = back.T[0][-1]
    
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
                return ValueError, k
            
            raise ValueError, k
        
        # Append to momenta values for spline calculation
        kVals = np.append(kVals, k)
    
    # Check momenta are strictly monotonically increasing
    if not all(kVals[i] < kVals[i + 1] for i in range(len(kVals) - 1)):
        
        if errorReturn:
            return ValueError, kVals
        
        raise ValueError, kVals
    
    kExit = kVals[Nsteps]
    
    # Get pivot scale, midpoint in spline
    if kPivot is None:
        kPivot = kVals[Nsteps]
        
    else:
        assert type(kPivot) == float, "kPivot parameter must be float: {}".format(kPivot)
    
    # Build power spectra values
    pZetaVals = np.array([])
    
    # Iterate over k-scales
    for k in kVals:
        
        # Build initial conditions for mode based on massless condition
        Nstart, ICs = ICsBM(subevo, k, back, pvals, MTE)
        
        # Check ICs are valid
        if type(ICs) != np.ndarray and np.isnan(ICs):
            
            if errorReturn:
                return ValueError, ICs
            
            raise ValueError, ICs
        
        # If no time out is given, assume no kill time
        if tmax is None:
            tmax = -1
        
        twoPf = MTE.sigEvolve(np.array([Nstart, Nend]), k, ICs, pvals, tols, False, tmax, True)
        
        # Handle failed calculation
        if type(twoPf) is tuple:
            
            if errorReturn:
                return ValueError, twoPf
            
            raise ValueError, twoPf
        
        # Append power spedctrum at end of background
        pZetaVals = np.append(pZetaVals, twoPf.T[1][-1])
    
    #  Build log arrays in k and Pzeta
    arrLogK = np.log(kVals / kPivot)
    arrLogPz = np.log(pZetaVals)
    
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


def fNL(back, pvals, Nexit, tols, subevo, alpha, beta, MTE, errorReturn=False, tmax=None):
    """
    
    Simple wrapper for the reduced bispectrum fNL at the end of inflation, subject to a configuration
    defined by alpha, beta (Fergusson & Shellard), centred about the horizon exit time Nend - Nexit
    
    """
    
    # If tmax is none, assign no maximum integration time
    if tmax is None: tmax_3pf = -1
    
    Nend = back.T[0][-1]
    Npivot = Nend - Nexit
    
    kExit = kexitN(Npivot, back, pvals, MTE)
    
    # Build Fourier triangle via Fergusson Shellard convention
    k1 = kExit / 2. - beta * kExit / 2.
    k2 = kExit * (1. + alpha + beta) / 4.
    k3 = kExit * (1. - alpha + beta) / 4.
    
    # Find largest scale; exits horizon first hence defines ICs
    kmin = np.min([k1, k2, k3])
    Nstart, ICs = ICsBM(subevo, kmin, back, pvals, MTE)
    
    if type(ICs) != np.ndarray and np.isnan(ICs):
        
        if errorReturn:
            return ValueError, ICs
        
        raise ValueError, ICs
    
    # Compute three-point function up until end of background
    threePt = MTE.alphaEvolve(np.array([Nstart, Nend]), k1, k2, k3, ICs, pvals, tols, True, tmax_3pf, True)
    
    # If flag has been returned when computing 3pf, return flag data
    if type(threePt) is tuple:
        
        if errorReturn:
            return ValueError, threePt
        
        raise ValueError, threePt
    
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