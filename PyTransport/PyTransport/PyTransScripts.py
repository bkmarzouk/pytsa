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

import numpy as np
from scipy import interpolate
import timeit
import sys
from gravtools_pyt import curvatureObject
import os
import pickle as pk
import gc

from scipy.interpolate import UnivariateSpline as US
import sympy as sym


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
    alp = threePtOut[:, 1 + 4 + 2 * nF + 6 * 2 * nF * 2 * nF:]
    
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


def ICsBM(NBMassless, k, back, params, MTE):
    nF = np.size(back[0, 1:]) // 2
    massEff = -1;
    
    # calculate the element of back for which -k^2/a^2 + M^2 for M the largest eigenvalue of mass matrix
    jj = 0
    while (massEff < 0 and jj < np.size(back[:, 0]) - 1):
        w, v = np.linalg.eig(MTE.ddV(back[jj, 1:1 + nF], params))
        eigen = np.max(w)
        massEff = -k ** 2 * np.exp(-2.0 * back[jj, 0]) + eigen
        jj = jj + 1
    if jj == np.size(back[:, 0]):
        print ("\n\n\n\n warning massless condition not found \n\n\n\n")
        return np.nan, np.nan
    NMassless = back[jj - 2, 0]
    backExitMinus = np.zeros(2 * nF)
    ll = 0
    Ncond = -1.
    while (Ncond < 0.0 and ll < np.size(back[:, 0]) - 1):
        Ncond = back[ll, 0] - (NMassless - NBMassless)
        ll = ll + 1
    
    if ll == np.size(back[:, 0]) or (NMassless - back[0, 0]) < NBMassless:
        print ("\n\n\n\n warning initial condition not found \n\n\n\n")
        return np.nan, np.nan
    
    NexitMinus = back[ll - 2, 0]
    backExitMinus = back[ll - 2, 1:]
    
    return NexitMinus, backExitMinus


# this script finds initial conditions at least NBExit e-folds before horizon exit of k
# back must finely sample the backgroudn evolution for the initial conditions to be close to exactly NBMassless before

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
                                  True)  # all data from three point run goes into threePt array
        zzOut = np.append(zzOut, twoPt[-1, 1])
        times = np.append(times, timeit.default_timer() - start_time)
    
    return zzOut, times


# calculates the power spectrum at each element in kA at the end of the background evolution (back) in a manner suitable to be called over many processes
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


# calculates the power spectrum and bisecpturm in equilateral configuration at each element in kA at the end of the background evolution (back)
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


# calculates the power spectrum and bisecpturm in equilateral configuration at each element in kA at the end of the background evolution (back) in a manner suitable to be run accoss many processes
def eqSpecMpi(kA, back, params, NB, tols, MTE):
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    
    rank = comm.Get_rank()
    size = comm.Get_size()
    points = np.size(kA)
    num = points / size;
    
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


# calcualtes the bispectrum in the alpha beta notation for a given kt at every value of the alphaIn and betIn arrays. The bispectrum is given an nsnaps times always incuding the final time of the evolution (back)
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


def ExtendedBackEvolve(initial, params, MTE, Nstart=0, Next=1, adpt_step=4e-3,
                       tols=np.array([1e-12, 1e-12]), tmax_bg=-1):
    """ Simple iterative procedure to extend canonical background evolution past epsilon;
        differs from exit=True routine in backEvolve by attempting to re-integrate after integrator limit """
    
    # Define number of iterations to attempt
    n_iter = Next / adpt_step

    """
    
    Integrator flag defs.
    
    NOTE: Failure upon integrating background quantities return tuple: (flag, Efold)

    int SHORT = -50;          // Inflation too short
    int KEXIT = -49;          // Unable to find Fourier mode
    int FEOI = -48;           // Integration failure in feoi
    int BACK = -47;           // Integration failure in background
    int VIOLATED = -46;       // Field space position violates model
    int ETERNAL = -45;        // Unable to find end of inflation
    int TIMEOUT = -44;        // Integration time exceeded
    
    """
    
    # Get fiducial end of inflation, i.e. when epsilon=1
    Nepsilon = MTE.findEndOfInflation(initial, params, tols, Nstart, 12000, tmax_bg)
    
    # Return flag
    if type(Nepsilon) is tuple: return Nepsilon[0]
    
    # Define initial efolding range and compute background
    Nspace_init = np.linspace(Nstart, Nepsilon, 10000)
    BG_epsilon = MTE.backEvolve(Nspace_init, initial, params, tols, True, tmax_bg)

    # If extended integration immediately fails, return background up until integrator limit
    if type(BG_epsilon) is tuple and BG_epsilon[0] in [-47, -44]: return BG_epsilon
    
    # We will store extensions to the background evolution in the following list
    extensions = [BG_epsilon]
    
    # Iteration counter
    c = 0
    
    # Adapt to higher precision
    tols = np.array([1e-12, 1e-12])
    
    while c < n_iter:
        
        # Define new ICs with last background step
        N_init, bg_init = extensions[-1][-1][0], extensions[-1][-1][1:]
        
        # Extend target efold range
        N_space = np.linspace(N_init, N_init + adpt_step, 100)
        
        # Try and compute extension to background evolution
        bg_ext = MTE.backEvolve(N_space, bg_init, params, tols, False, tmax_bg)
        
        # If a integration flag is returned, break out of loop
        if type(bg_ext) is tuple and bg_ext[0] in [-47, -44]: break
        
        # Otherwise append extension to list
        extensions.append(bg_ext[1:]) # [1:] skips first step in extension (this is already in the list)
        
        # Increment counter
        c += 1
    
    # If no extensions are found, simply
    if len(extensions) == 1:
        return BG_epsilon, Nepsilon
    
    else:
        # Reconstruct extensions to 2d numpy array (as backEvolve)
        stacked = np.vstack((bg for bg in extensions))
        
        # return the background, as well as the e-folding where slow-roll was violated
        return stacked, Nepsilon


def evolveMasses(back, params, MTE, scale_eigs=False, hess_approx=False, covariant=False):
    """ Computes the mass matrix along a given background evolution, M^{I}_{J} """
    
    # Build empty array to populate with efold number + eigenvalues of mass-matrix at time step
    eigs = np.empty((len(back), 1+MTE.nF()))
    
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


def GetCurvatureObject(MTE):
    """ Returns the class 'curvature' object associate with the MTE module """
    
    # Find curvature records directory and create file instance to lead curvature class obj.
    dir = os.path.dirname(__file__)
    curv_dir = os.path.join(dir, 'PyTrans', 'CurvatureRecords')
    curv_name = str(MTE).split("'")[1][7:] + ".curvature"
    curv_path = os.path.join(curv_dir, curv_name)
    curv_file = open(curv_path, 'rb')
    
    with curv_file:
        curv_obj = pk.load(curv_file)
    
    return curv_obj