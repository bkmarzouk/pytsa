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
from gravtools_pyt import Curvature
import os
import pickle as pk

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


def derivative(x, y, k=3, s=0, M=100):
    y_spl = interpolate.UnivariateSpline(x, y, k=k, s=s)
    
    x_fine = np.linspace(x[0], x[-1], M * len(x));
    s = x_fine[1] - x_fine[0]
    y_fine = y_spl(x_fine)
    
    dy = [(y_fine[i + 1] - y_fine[i]) / s for i in range(len(y_fine) - 1)]
    
    dy_spl = interpolate.UnivariateSpline(x_fine[:-1], dy, k=k, s=s)
    return y_spl, y_spl.derivative()


def dotdotfields(back, MTE, params):
    """ Computes second time derivatives of the fields. Uses analytic form of background evolution equations. """
    
    # backT = back.T
    # nF = MTE.nF()
    #
    # # Get Hubble values to transform N derivatives to t derivatives
    # Hseries = np.asarray([MTE.H(np.array([step[1:]][0]), params) for step in back])
    #
    # dV = [MTE.dV(fvals[1:1 + nF], params) for fvals in back]
    #
    # phidotdot = np.vstack(
    #     (-3. * Hseries * backT[1 + nF + i] - np.asarray(dV).T[i] for i in range(nF))
    # ).T
    #
    # return phidotdot
    
    """ Returns second time derivatives of the fields using analytic form of the background equations of motion
        as a contravariant object"""
    
    # Find curvature records directory and create file instance to lead curvature class obj.
    dir = os.path.dirname(__file__)
    curv_dir = os.path.join(dir, 'PyTrans', 'CurvatureRecords')
    curv_name = str(MTE).split("'")[1][7:] + ".curvature"
    curv_path = os.path.join(curv_dir, curv_name)
    curv_file = open(curv_path, 'r')
    with curv_file:
        curv_obj = pk.load(curv_file)
    
    nF = MTE.nF()

    # Get Hubble values to transform N derivatives to t derivatives
    Hubble = np.array([MTE.H(np.array([step[1:]][0]), params) for step in back])
    
    # Get field time derivatives from background data
    phi = np.array([step[1:1+nF] for step in back])
    phidot = np.array([step[1 + nF:1+2*nF] for step in back])
    
    subs = []
    for s in range(len(back)):
        # = {}
        
        for i in range(nF):
            subs[curv_obj.coords[i]] = phi[i]

    # Get derivatives of the potential from background data
    dV_ = np.array([[MTE.dV(fvals[1:1 + nF], params) for fvals in back]])
    
    # Raise an index from dV
    Ginv = curv_obj.metric_inverse
    dV = np.array([Ginv[a, b]*dV_[b] for a in range(nF) for b in range(nF)])
    
    # Compute second time derivative of the fields
    phidotdot = np.vstack(
        (-3. * Hubble[i] * phidot[i] - dV[i] for i in range(nF))
    )

    return phidotdot


def EvolveKE(back, params, MTE):
    """ Computes the Kinetic Energy of the fields along the background trajectory """
    
    # Find curvature records directory and create file instance to lead curvature class obj.
    dir = os.path.dirname(__file__)
    curv_dir = os.path.join(dir, 'PyTrans', 'CurvatureRecords')
    curv_name = str(MTE).split("'")[1][7:] + ".curvature"
    curv_path = os.path.join(curv_dir, curv_name)
    curv_file = open(curv_path, 'r')
    with curv_file:
        curv_obj = pk.load(curv_file)
    
    # Inverse metric
    metric = curv_obj.metric
    
    # Get number of fields and params
    nF = MTE.nF()
    nP = MTE.nP()
    
    print "-- Building dict"
    
    # Build lambda functions for metric
    pdict = dict()
    for i in range(nP):
        pdict[curv_obj.params[i]] = params[i]
    
    print "-- Building lambda functions"
    
    metric_lambdas = np.empty((nF, nF), dtype=object)
    
    for i in range(nF):
        for j in range(nF):
            metric_lambdas[i, j] = sym.lambdify([f for f in curv_obj.coords], metric[i, j].subs(pdict), "numpy")
    
    print "-- Computing K.E. evolution"
    
    
    # Build dictionary for symbolic substitutions
    def g_ab(a, b, step):
        fields = list(step[1:1 + nF])
        m_eval = metric_lambdas[a, b](*fields)
        return m_eval
    
    
    # Kinetic Energy = 1/2 * g_{ab} * dt phi^a dt phi^b
    backKE = np.array(
        [0.5 * sum([g_ab(a, b, step) * step[1 + nF + a] * step[1 + nF +  b] for a in range(nF) for b in range(nF)]) for
         step in back])
    
    # return 2 x N numpy array, with the the first col containing efold number and the second the fields K.E.
    return np.vstack([back.T[0], backKE]).T


def KineticEnergy(fields, dotfields, params, MTE, curv_obj):
    """ Computes the kinetic energy of the fields. Assumes canonical lagrangian kinetic term """
    
    if curv_obj is None:
        # Find curvature records directory and create file instance to lead curvature class obj.
        dir = os.path.dirname(__file__)
        curv_dir = os.path.join(dir, 'PyTrans', 'CurvatureRecords')
        curv_name = str(MTE).split("'")[1][7:] + ".curvature"
        curv_path = os.path.join(curv_dir, curv_name)
        curv_file = open(curv_path, 'r')
        with curv_file:
            curv_obj = pk.load(curv_file)
    
    # Inverse metric
    metric = curv_obj.metric
    
    nF = MTE.nF()
    nP = MTE.nP()
    KE = 0
    
    # Build dictionary for symbolic substitutions
    sub_dict = dict()
    for i in range(nF):
        sub_dict[curv_obj.coords[i]] = fields[i]
    if curv_obj.params is not None:
        for i in range(nP):
            sub_dict[curv_obj.params[i]] = params[i]
    
    # Iterate over field indices
    for I in range(nF):
        for J in range(nF):
            # Contract field velocities with field metric to get kinetic energy
            KE += 0.5 * metric[I, J].subs(sub_dict) * dotfields[I] * dotfields[J]
    
    return KE


def ExtendedBackEvolve(initial, params, MTE, Nstart=0, Next=1, adpt_step=1e-4, tols=np.array([1e-10, 1e-10])):
    # Define number of iterations to attempt
    n_iter = Next / adpt_step
    
    # Get fiducial end of inflation, i.e. when epsilon=1
    Nepsilon = MTE.findEndOfInflation(initial, params, tols, Nstart, 12000)
    
    # If integration failure, return None
    if type(Nepsilon) == type('S32'):
        return None
    
    # Define initial efolding range and compute background
    Nspace_init = np.linspace(Nstart, Nepsilon, 10000)
    BG_epsilon = MTE.backEvolve(Nspace_init, initial, params, tols, True)
    
    idx_epsilon = len(BG_epsilon)
    
    # We will store extensions to the background evolution in the following list
    extensions = [BG_epsilon]
    
    c = 0
    tols = np.array([1e-12, 1e-12])  # Adapt to higher precision
    while c < n_iter:
        
        # Define new initial conditions to be the field data at the previous background segment
        N_init = extensions[-1][-1][0]
        bg_init = extensions[-1][-1][1:]
        
        # define new efolding range and compute extension to background
        N_space = np.linspace(N_init, N_init + adpt_step, 100)
        bg_ext = MTE.backEvolve(N_space, bg_init, params, tols, False)
        
        # If we fail to compute the backgroud, break out of iterative cycle
        if type(bg_ext) != np.ndarray:
            break
        
        # Otherwise add segment to the list of background data
        extensions.append(bg_ext[1:])
        
        c += 1
    
    if len(extensions) == 1:
        
        print "No Extension", BG_epsilon[-1][0]
        
        return BG_epsilon, Nepsilon
    
    else:
        # Stack together into numpy ndarray, structured consistently with typical backEvolve
        stacked = np.vstack((bg for bg in extensions))
        
        print "Extended background by N = {}".format(stacked.T[0][-1] - Nepsilon)
        
        # return the background, as well as the efolding where slow-roll was violated
        return stacked, Nepsilon


def MijEvolve(back, params, MTE, DropKineticTerms=False, scale_eigs=False):
    """ Computes the mass matrix along a given background evolution, M^{I}_{J} """
    
    # Find curvature records directory and create file instance to lead curvature class obj.
    dir = os.path.dirname(__file__)
    curv_dir = os.path.join(dir, 'PyTrans', 'CurvatureRecords')
    curv_name = str(MTE).split("'")[1][7:] + ".curvature"
    curv_path = os.path.join(curv_dir, curv_name)
    curv_file = open(curv_path, 'rb')
    
    with curv_file:
        curv_obj = pk.load(curv_file)
    
    # Number of fields and parameters
    nF = MTE.nF()
    nP = MTE.nP()
    
    """ Get evolution of background quantities """
    
    # Potential & derivatives efold evolution
    V_evo = [MTE.V(step[1:1 + nF], params) for step in back] # Time series of potential
    dV_evo = [MTE.dV(step[1:1 + nF], params) for step in back]  # Time series of 1st deriv. pot.
    ddV_evo = [MTE.ddV(step[1:1 + nF], params) for step in back]  # Time series of 2nd deriv. pot.
    
    # Hubble rate
    hubble = [MTE.H(np.concatenate([step[1:1 + nF], step[1 + nF:1+2*nF]]), params) for step in back]
    
    # Transpose background evolution into horizontal components
    backT = back.T

    # Get symbolic definitions for background
    f_syms = curv_obj.coords        # fields
    v_syms = sym.symarray("v", nF)  # velocities
    p_syms = curv_obj.params        # parameters
    fp_syms = [f for f in f_syms] + [v for v in v_syms] + [p for p in p_syms]
    
    # Get Christoffel symbols and Riemann tensor
    csyms = curv_obj.Csyms
    rsyms = curv_obj.Rsyms
    
    # Get field space metric and inverse
    G = curv_obj.metric
    Ginv = curv_obj.metric_inverse
    
    # Build array and populate with field space lambda functions
    Glambdas = np.empty(np.shape(G), dtype=object)
    Ginvlambdas = np.empty(np.shape(Ginv), dtype=object)
    clambdas = np.empty(np.shape(csyms), dtype=object)
    rlambdas = np.empty(np.shape(rsyms), dtype=object)

    # Define coordinate index range to iterate over
    rnf = range(nF)

    print "-- lambdifying"
    
    # Lambdify symbolic expressions
    for a in rnf:
        for b in rnf:
            Glambdas[a, b] = sym.lambdify(fp_syms, G[a, b], "numpy")
            Ginvlambdas[a, b] = sym.lambdify(fp_syms, Ginv[a, b], "numpy")
            for c in rnf:
                clambdas[a, b, c] = sym.lambdify(fp_syms, csyms[a, b, c], "numpy")
                for d in rnf:
                    rlambdas[a, b, c, d] = sym.lambdify(fp_syms, rsyms[a, b, c, d], "numpy")
    print "-- done"
    
    eig_evo_raw = [] # Record raw eigenvalues of mass-matrix
    eig_evo_hub = [] # Record Hubble normalized eigenvalues of mass-matrix
    
    s = 0
    lb = len(back)
    
    def lp(*terms):
        # Logarithmic product
        
        signs = np.array([np.sign(item) for item in terms])
        logs  = np.array([np.log(abs(item)) for item in terms])
        
        return np.prod(signs)*np.exp(np.sum(logs))
    
    for s in range(lb):
        
        print s+1, "/", lb
        
        # Background step
        step = back[s]
        
        # Fields and field time derivatives
        phi = step[1:1+nF]
        dtphi = step[1+nF:1+2*nF]
        
        # Values for symbolic evaluations
        sym_subs = [x for x in phi] + [x for x in dtphi] + [x for x in params]
        
        # Hubble rate at current step
        H = hubble[s]
        
        # Get values for potential and derivatives
        V = V_evo[s]
        dV = dV_evo[s]
        ddV = ddV_evo[s]
        
        # Initialize covariant Hessian, Riemann and Kinetic term
        covhess = np.zeros((nF, nF), dtype=float)
        riemann = np.zeros((nF, nF), dtype=float)
        kinetic = np.zeros((nF, nF), dtype=float)
        
        def covdt_dtphi_up(idx):
            return -3.*H*dtphi[idx] - sum([Ginvlambdas[k, idx](*sym_subs)*dV[k] for k in rnf])
        
        def covdt_dtphi_down(idx):
            return sum([-3.*H*Glambdas[idx, k](*sym_subs)*dtphi[k] for k in rnf]) - dV[idx]
        
        def dtphi_down(idx):
            return sum([Glambdas[idx, k](*sym_subs)*dtphi[k] for k in rnf])
        
        Hdot = V - 3.*H**2
        
        for a in rnf:
            for b in rnf:
                
                covhess[a, b] +=  sum([Ginvlambdas[k, a](*sym_subs)*ddV[k, b]
                                      for k in rnf])
                
                covhess[a, b] += -sum([Ginvlambdas[a, l](*sym_subs)*clambdas[k, l, b](*sym_subs)*dV[k]
                                       for k in rnf for l in rnf])
                
                riemann[a, b] += -sum([Ginvlambdas[a, m](*sym_subs)*rlambdas[m, k, l, b](*sym_subs)*dtphi[k]*dtphi[l]
                                        for k in rnf for l in rnf for m in rnf])
                
                kinetic[a, b] += -(3.-Hdot / H / H)*dtphi[a]*dtphi_down(b) - (dtphi[a]*covdt_dtphi_down(b) +
                                                                            covdt_dtphi_up(a)*dtphi_down(b))/H
        
        # Combine terms to produce mass-squared-matrix
        MIj = covhess + riemann + kinetic

        # Compute eigenvalues
        eig, eigv = np.linalg.eig(MIj)
        
        # Sort eigenvalues into
        eig_sorted = np.sort(eig)
        eig_sorted_sqrt = np.sqrt(np.abs(eig_sorted)) * np.sign(eig_sorted)

        eig_sorted_hubble = eig_sorted_sqrt / H
        
        eig_evo_raw.append(eig_sorted_sqrt)
        eig_evo_hub.append(eig_sorted_hubble)

    eigstack_raw = np.vstack([np.asarray(eig) for eig in eig_evo_raw])
    eigstack_raw = np.hstack([np.asarray(backT[0]).reshape(len(back), 1), eigstack_raw])

    eigstack_hub = np.vstack([np.asarray(eig) for eig in eig_evo_hub])
    eigstack_hub = np.hstack([np.asarray(backT[0]).reshape(len(back), 1), eigstack_hub])

    # Return matrix in typical PyTransport format
    return eigstack_raw, eigstack_hub