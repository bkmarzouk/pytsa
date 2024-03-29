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


# python code contains some useful scripts to use with the compiled pyt module.

import timeit

import numpy as np
from scipy import interpolate
from scipy.interpolate import UnivariateSpline
from . import spline_tools as splt


def print_warning(msg):
    print(f"[WARNING] {msg}")


def unPackAlp(threePtOut, MTE):
    nF = MTE.nF()

    _size = 1 + 4 + 2 * nF + 6 * 2 * nF * 2 * nF + 2 * nF * 2 * nF * 2 * nF

    if np.size(threePtOut[0, :]) != _size:
        print_warning(
            "warning array you asked to unpack is not of " "correct dimension"
        )
        return np.nan, np.nan, np.nan, np.nan, np.nan

    zetaMs = threePtOut[:, 1:5]
    sig1R = threePtOut[:, 1 + 4 + 2 * nF : 1 + 4 + 2 * nF + 2 * nF * 2 * nF]
    sig2R = threePtOut[
        :,
        1
        + 4
        + 2 * nF
        + 2 * nF * 2 * nF : 1
        + 4
        + 2 * nF
        + 2 * 2 * nF * 2 * nF,
    ]
    sig3R = threePtOut[
        :,
        1
        + 4
        + 2 * nF
        + 2 * 2 * nF * 2 * nF : 1
        + 4
        + 2 * nF
        + 3 * 2 * nF * 2 * nF,
    ]
    alpha = threePtOut[:, 1 + 4 + 2 * nF + 6 * 2 * nF * 2 * nF :]

    return (
        zetaMs,
        np.reshape(sig1R, (np.size(threePtOut[:, 0]), 2 * nF, 2 * nF)),
        np.reshape(sig2R, (np.size(threePtOut[:, 0]), 2 * nF, 2 * nF)),
        np.reshape(sig3R, (np.size(threePtOut[:, 0]), 2 * nF, 2 * nF)),
        np.reshape(alpha, (np.size(threePtOut[:, 0]), 2 * nF, 2 * nF, 2 * nF)),
    )


def unPackSig(twoPtOut, MTE):
    nF = MTE.nF()
    if np.size(twoPtOut[0, :]) != 1 + 1 + 2 * nF + 2 * nF * 2 * nF:
        print_warning(
            "warning array you asked to unpack is not of " "correct dimension"
        )

        return np.nan, np.nan
    zeta = twoPtOut[:, 1]
    sig = twoPtOut[:, 1 + 1 + 2 * nF : 1 + 1 + 2 * nF + 2 * nF * 2 * nF]

    return zeta, np.reshape(sig, (np.size(twoPtOut[:, 0]), 2 * nF, 2 * nF))


def arr_eval_spl_N(N, arr, k=3):
    """
    Estimate time ordered array of values at exact time N, via interpolating
    through array data on the interval [N-0.5, N+0.5]

    :param N: evaluation time (efolds, N)
    :param arr: time ordered array of evolving quantities.
                Zeroth column should contain 'N' values.
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
    Estimate time ordered array of values at exact time N,
    via interpolating through array data on the interval [N-0.5, N+0.5]

    :param N: evaluation time (efolds, N)
    :param arr: time ordered array of evolving quantities.
                Zeroth column should contain 'N' values.
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
        raise ValueError(
            "Unable to locate appropriate value range within input array"
        )

    out = arr[
        np.max([0, idx - spl_steps // 2]) : np.min(
            [len(arr), idx + spl_steps // 2]
        )
    ]

    if len(out) < 2:
        raise ValueError(
            "Unable to locate appropriate value range within input array"
        )

    if k >= len(out):
        k = len(out) - 1

    val_spl = out.T[col_idx]

    spls = [UnivariateSpline(val_spl, col, k=k) for col in out.T]

    return np.array([s(val) for s in spls])


def arr_eval_nearest_N(N, arr):
    """
    Estimate time ordered array of values at approximate time N,
    via locating the nearest available time step

    :param N: evaluation time (efolds, N)
    :param arr: time ordered array of evolving quantities.
                Zeroth column should contain 'N' values.
                Zeroth should contain quantities at the earliest time
    :return: array quantities at the time step N
    """

    Nseries = arr.T[0]

    idx_out = None

    for idx, Nstep in enumerate(Nseries):
        if Nstep > N:
            idx_out = idx

    return arr[idx_out]


def adjust_back(bg, Nr=80.0, reposition=True, remove_duplicate_steps=True):
    """

    Repositions background evolution to begin Nr efolds before the end of
    inflation.

    (Helps avoid momenta growing exponentially large.)

    If Nr exceeds the length of inflation, the original background is returned.

    Finally, we also ensure that the first 10 efolds of inflation are finely
    sampled enough s.t. accurate estimates of ICsBM can be made

    """

    # For very high density stepper, sometimes there is no change in
    # time-series values from N[ii]->N[ii+1]
    # Remove these array values to avoid numerical stiffness
    # and non-monotonicity
    if remove_duplicate_steps:
        for ii in range(1, len(bg))[::-1]:
            if np.all(bg[ii] == bg[ii - 1]):
                bg = np.delete(bg, ii, axis=0)

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
                    # Stop background repositioning,
                    # and redefine N efolution to start at zero from this point
                    out[:, 0] -= out[0][0]

                    break
        else:
            out = bg

    else:
        out = bg

    # Cross-check final field values haven't changed
    assert np.all(out[-1][1:] == bg[-1][1:]), [out[-1][1:], bg[-1][1:]]

    # Finally, we check that the background hasn't been too finely sampled,
    # s.t back[ii] == back[ii+1]
    killRows = np.array(
        [ii for ii in range(len(out) - 1) if np.all(out[ii] == out[ii + 1])]
    )

    killRows = killRows[::-1]

    for ii in killRows:
        out = np.delete(out, ii, axis=0)

    return out


def find_massless_condition(MTE, N_before_massless, k, back, params):
    kexit_evolution = np.ones(len(back), dtype=float) * k**2
    _mass_evolution = evolve_mass_spectrum(back, params, MTE)

    mass_evolution = np.zeros(kexit_evolution.shape, dtype=float)

    for idx, step in enumerate(_mass_evolution):
        mass_evolution[idx] = step[1:].max()

    for idx, step in enumerate(back):
        N, fdf = step[0], step[1:]
        kexit_evolution[idx] /= MTE.H(fdf, params) ** 2 * np.exp(2 * N)

    diff_evolution = np.abs(mass_evolution - kexit_evolution)

    idx_massless = np.where(diff_evolution == diff_evolution.min())[0][0]

    Nstart = back.T[0][idx_massless] - N_before_massless

    if Nstart < 0:
        return np.nan, np.nan

    new_ics = splt.approx_row_closest(Nstart, back, 0)

    return new_ics[0], new_ics[1:]


def ICsBM(NBMassless, k, back, params, MTE):
    nF = np.size(back[0, 1:]) // 2
    massEff = -1

    # calculate the element of back for which -k^2/a^2 + M^2 for M the largest
    # eigenvalue of mass matrix
    jj = 0
    while massEff < 0 and jj < np.size(back[:, 0]) - 1:
        w, v = np.linalg.eig(MTE.ddV(back[jj, 1 : 1 + nF], params))
        eigen = np.max(w)
        massEff = -(k**2) * np.exp(-2.0 * back[jj, 0]) + eigen
        jj = jj + 1
    if jj == np.size(back[:, 0]):
        print("\n\n\n\n warning massless condition not found \n\n\n\n")
        return np.nan, np.nan
    NMassless = back[jj - 2, 0]
    ll = 0
    Ncond = -1.0
    while Ncond < 0.0 and ll < np.size(back[:, 0]) - 1:
        Ncond = back[ll, 0] - (NMassless - NBMassless)
        ll = ll + 1

    if ll == np.size(back[:, 0]) or (NMassless - back[0, 0]) < NBMassless:
        print("\n\n\n\n warning initial condition not found \n\n\n\n")
        return np.nan, np.nan

    NexitMinus = back[ll - 2, 0]
    backExitMinus = back[ll - 2, 1:]

    return NexitMinus, backExitMinus


def ICsBE(NBExit, k, back, params, MTE):
    nF = np.size(back[0, 1:]) // 2
    kvaH = -1.0
    jj = 0
    while kvaH < 0.0 and jj < np.size(back[:, 0]) - 1:
        H = MTE.H(back[jj, 1 : 1 + 2 * nF], params)
        kvaH = -k + np.exp(back[jj, 0]) * H
        jj = jj + 1
    if jj == np.size(back[:, 0]):
        print("\n\n\n\n warning exit condition not found \n\n\n\n")
        return np.nan, np.nan
    NExit = back[jj - 2, 0]
    ll = 0
    Ncond = -1
    while Ncond < 0 and ll < np.size(back[:, 0]) - 1:
        Ncond = back[ll, 0] - (NExit - NBExit)
        ll = ll + 1

    if ll == np.size(back[:, 0]) or (NExit - back[0, 0]) < NBExit:
        print_warning("Initial condition not found")
        return np.nan, np.nan
    NexitMinus = back[ll - 2, 0]
    backExitMinus = back[ll - 2, 1:]

    return NexitMinus, backExitMinus


# find the earliest condition between the massless one and the horizon exit one
def ICs(NB, k, back, params, MTE):
    NBEs, fieldBE = ICsBE(NB, k, back, params, MTE)
    NBMs, fieldBM = ICsBM(NB, k, back, params, MTE)

    if NBEs < NBMs:
        return NBEs, fieldBE
    return NBMs, fieldBM


# calculates the power spectrum at each element in kA at the end of the
# background evolution (back)
def pSpectra(kA, back, params, NB, tols, MTE, tmax_2pf=-1, flagReturn=False):
    zzOut = np.array([])
    times = np.array([])
    num = np.size(kA)

    for ii in range(0, num):
        print(
            "\n \n \n performing "
            + str(ii + 1)
            + " of "
            + str(num)
            + "\n \n \n"
        )
        k = kA[ii]
        Nstart, backExitMinus = ICs(NB, k, back, params, MTE)
        start_time = timeit.default_timer()

        if Nstart == np.nan:
            twoPt = np.empty((2, 2))
            twoPt[:] = np.nan
        else:
            t = np.linspace(Nstart, back[-1, 0], 10)
            # run solver for this triangle
            twoPt = MTE.sigEvolve(
                t, k, backExitMinus, params, tols, True, tmax_2pf, flagReturn
            )  # all data from three point run goes into threePt array
        zzOut = np.append(zzOut, twoPt[-1, 1])
        times = np.append(times, timeit.default_timer() - start_time)

    return zzOut, times


# calculates the power spectrum at each element in kA at the end of the
# background evolution (back) in a manner suitable over many processes
def pSpecMpi(kA, back, params, NB, tols, MTE):
    from mpi4py import MPI

    comm = MPI.COMM_WORLD

    rank = comm.Get_rank()
    size = comm.Get_size()
    points = np.size(kA)
    num = points / size

    if float(points) / size != float(points // size):
        if rank == 0:
            print_warning(
                "warning! number of points is divisable by number"
                "of processes, exiting"
            )
        return np.nan, np.nan

    kOutL = kA[rank * num : rank * num + num]

    zzL, timesL = pSpectra(kOutL, back, params, NB, tols, MTE)

    if rank != 0:
        comm.Send(zzL, dest=0)
        comm.Send(timesL, dest=0)

    if rank == 0:
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


# calculates the power spectrum and bisecpturm in equilateral configuration
# at each element in kA at the end of the background evolution (back)
def eqSpectra(kA, back, params, NB, tols, MTE):
    zzzOut = np.array([])
    zzOut = np.array([])
    times = np.array([])
    num = np.size(kA)

    for ii in range(0, num):
        print(
            "\n \n \n performing "
            + str(ii + 1)
            + " of "
            + str(num)
            + "\n \n \n"
        )
        k = kA[ii]
        Nstart, backExitMinus = ICs(NB, k, back, params, MTE)
        t = np.linspace(Nstart, back[-1, 0], 10)
        k1 = k
        k2 = k
        k3 = k
        # run solver for this triangle
        start_time = timeit.default_timer()

        if Nstart == np.nan:
            threePt = np.empty((2, 5))
            threePt[:] = np.nan
        else:
            threePt = MTE.alphaEvolve(
                t, k1, k2, k3, backExitMinus, params, tols, True
            )  # all data from three point run goes into threePt array
        zzOut = np.append(zzOut, threePt[-1, 1])
        zzzOut = np.append(zzzOut, threePt[-1, 4])
        times = np.append(times, timeit.default_timer() - start_time)
    return zzOut, zzzOut, times


# calculates the power spectrum and bisecpturm in equilateral configuration at
# each element in kA at the end of the background evolution (back) in a manner
# suitable to be run accoss many processes
def eqSpecMpi(kA, back, params, NB, tols, MTE):
    from mpi4py import MPI

    comm = MPI.COMM_WORLD

    rank = comm.Get_rank()
    size = comm.Get_size()
    points = np.size(kA)
    num = points / size

    if float(points) / size != float(points // size):
        if rank == 0:
            print_warning(
                "warning! number of points is divisable by "
                "number of processes, exiting"
            )
        return (np.empty, np.empty, np.empty)

    kOutL = kA[rank * num : rank * num + num]

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


# calcualtes the bispectrum in the alpha beta notation for a given kt at every
# value of the alphaIn and betIn arrays. The bispectrum is given an nsnaps
# times always incuding the final time of the evolution (back)
def alpBetSpectra(kt, alphaIn, betaIn, back, params, NB, nsnaps, tols, MTE):
    Hin = np.zeros(np.size(back[:, 0]))
    for jj in range(0, np.size(back[:, 0])):
        Hin[jj] = MTE.H(back[jj, 1:], params)
    aH = np.exp(back[:, 0]) * Hin
    positions = np.argsort(aH)
    Nexit = interpolate.splev(
        kt / 3.0,
        interpolate.splrep(aH[positions], back[positions, 0], s=1e-15),
        der=0,
    )
    Nend = back[-1, 0]
    snaps = np.linspace(Nexit - (NB - 0.1), Nend, nsnaps)
    if nsnaps == 1 or nsnaps == 0:
        snaps = np.array([Nend])
        nsnaps = 1

    biAOut = np.zeros([np.size(alphaIn), np.size(betaIn), np.size(snaps)])
    zz1 = np.zeros([np.size(alphaIn), np.size(betaIn), np.size(snaps)])
    zz2 = np.zeros([np.size(alphaIn), np.size(betaIn), np.size(snaps)])
    zz3 = np.zeros([np.size(alphaIn), np.size(betaIn), np.size(snaps)])
    times = np.zeros([np.size(alphaIn), np.size(betaIn)])

    for alpha_in_idx in range(0, np.size(alphaIn)):
        alpha = alphaIn[alpha_in_idx]
        for j in range(0, np.size(betaIn)):
            print(
                "\n \n \n performing "
                + str(alpha_in_idx + 1)
                + " "
                + str(j + 1)
                + " of "
                + str(np.size(alphaIn))
                + " "
                + str(np.size(betaIn))
                + "\n \n \n"
            )
            timebefore = timeit.default_timer()
            beta = betaIn[j]
            k1 = kt / 2.0 - beta * kt / 2.0
            k2 = kt / 4.0 * (1.0 + alpha + beta)
            k3 = kt / 4.0 * (1.0 - alpha + beta)
            if alpha > -(1 - beta) and alpha < 1 - beta:
                kM = min(k1, k2, k3)
                Nstart, backExitMinus = ICs(NB, kM, back, params, MTE)
                # run solver for this triangle
                t = np.concatenate((np.array([Nstart]), snaps))
                if Nstart == np.nan:
                    threePt = np.empty((2, 5))
                    threePt[:] = np.nan
                else:
                    threePt = MTE.alphaEvolve(
                        t, k1, k2, k3, backExitMinus, params, tols, True
                    )
                zzz = threePt[:, :5]

                for ii in range(1, nsnaps + 1):
                    biAOut[alpha_in_idx, j, ii - 1] = zzz[ii, 4]
                    zz1[alpha_in_idx, j, ii - 1] = zzz[ii, 1]
                    zz2[alpha_in_idx, j, ii - 1] = zzz[ii, 2]
                    zz3[alpha_in_idx, j, ii - 1] = zzz[ii, 3]
            else:
                for ii in range(0, nsnaps + 1):
                    biAOut[alpha_in_idx, j, ii - 1] = np.nan
                    zz1[alpha_in_idx, j, ii - 1] = np.nan
                    zz2[alpha_in_idx, j, ii - 1] = np.nan
                    zz3[alpha_in_idx, j, ii - 1] = np.nan

            times[alpha_in_idx, j] = timeit.default_timer() - timebefore
    return (biAOut, zz1, zz2, zz3, times, snaps)


# performs the same task as alpBetSpectra but in a manner suitable to be spread
# across many processes
def alpBetSpecMpi(kt, alpha, beta, back, params, NB, nsnaps, tols, MTE):
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    side = np.size(alpha)
    Nbefore = NB
    rank = comm.Get_rank()
    size = comm.Get_size()
    num = side / size

    if float(side) / size != float(side // size):
        if rank == 0:
            print_warning(
                "warning! number size of alpha must be divisable by "
                "number of processes, exiting"
            )
        return (np.empty, np.empty, np.empty, np.empty, np.empty)
    else:
        alphaL = alpha[rank * num : rank * num + num]

        BzL, Pz1L, Pz2L, Pz3L, timesL, snaps = alpBetSpectra(
            kt, alphaL, beta, back, params, Nbefore, nsnaps, tols, MTE
        )

        if rank != 0:
            comm.send(Pz1L, dest=0)
            comm.send(Pz2L, dest=0)
            comm.send(Pz3L, dest=0)
            comm.send(BzL, dest=0)
            comm.send(timesL, dest=0)

        if rank == 0:
            Bztot = np.zeros(
                [np.size(alpha), np.size(beta), np.size(BzL[0, 0, :])]
            )
            Pz1tot = np.zeros(
                [np.size(alpha), np.size(beta), np.size(BzL[0, 0, :])]
            )
            Pz2tot = np.zeros(
                [np.size(alpha), np.size(beta), np.size(BzL[0, 0, :])]
            )
            Pz3tot = np.zeros(
                [np.size(alpha), np.size(beta), np.size(BzL[0, 0, :])]
            )
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

                Bztot[jj * num : jj * num + num, :, :] = BzL
                Pz1tot[jj * num : jj * num + num, :, :] = Pz1L
                Pz2tot[jj * num : jj * num + num, :, :] = Pz2L
                Pz3tot[jj * num : jj * num + num, :, :] = Pz3L
                timestot[jj * num : jj * num + num, :] = timesL
            return (Bztot, Pz1tot, Pz2tot, Pz3tot, timestot, snaps)
        else:
            return (np.empty, np.empty, np.empty, np.empty, np.empty, snaps)


def kexitN(
    Nexit: float, back: np.ndarray, params: np.ndarray, MTE, fit="nearest"
):
    """
    Computes horizon exit momenta at efold Nexit

    :param Nexit: Horizon exit time
    :param back: background evolution
    :param params: model parameters
    :param MTE: PyTransport model
    :param fit: method for fitting / finding k value. 'nearest' finds closest
                step to Nexit; 'spl' forms spline estimate
    :return: Horizon exit mode, k_{exit}
    """
    assert Nexit >= back[0][0], Nexit
    assert Nexit <= back[-1][0], Nexit

    # TODO: Read in spline tools, remove dep on depric

    if fit == "spl":
        back_exit = arr_eval_spl_val(Nexit, back, 0)
    elif fit == "nearest":
        back_exit = arr_eval_nearest_N(Nexit, back)
    else:
        raise KeyError(fit)

    Hexit = MTE.H(back_exit[1:], params)

    k = np.exp(Nexit + np.log(Hexit))

    return k


def kexitPhi(PhiExit, n, back, params, MTE):
    nF = np.size(back[0, 1:]) // 2
    backExit = np.zeros(2 * nF)
    positions = np.argsort(back[:, n])
    Nexit = interpolate.splev(
        PhiExit,
        interpolate.splrep(back[positions, n], back[positions, 0], s=1e-15),
        der=0,
    )

    for i in range(1, 2 * nF + 1):
        backExit[i - 1] = interpolate.splev(
            Nexit, interpolate.splrep(back[:, 0], back[:, i], s=1e-15), der=0
        )
    k = np.exp(Nexit) * MTE.H(backExit, params)
    return k


def compute_Nexit_for_matching(
    MTE, back: np.ndarray, params: np.ndarray, k=0.002
):
    const = 55.75 - np.log(k / 0.05) + np.log(243.5 * 3**0.25)

    Hend = MTE.H(back[-1][1:], params)

    Hk = None

    for ii in range(1, len(back)):
        jj = ii - 1

        N_ii = back[ii - 1][0]
        N_jj = back[ii][0]

        H_ii = MTE.H(back[ii][1:], params)
        H_jj = MTE.H(back[jj][1:], params)

        k_ii = np.exp(N_ii + np.log(H_ii))
        k_jj = np.exp(N_jj + np.log(H_jj))

        if k_jj > k:
            Hk = H_ii if abs(k_ii - k) < abs(k_jj - k) else H_jj

            break

    if Hk is None:
        return np.nan

    N_before_end = const + np.log(Hk) - np.log(Hend) / 2

    out = back.T[0][-1] - N_before_end

    if out < 0:
        return np.nan

    return out


def matchKExitN(back, params, MTE, k=0.002):  # Remove
    const = 55.75 - np.log(k / 0.05)

    # Mpl = 2.435 × 10^18 GeV
    # Tre = 10^16 Gev

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
        print_warning(
            "warning unable to compute pivot time Nk with " "background"
        )

        return np.nan

    try:
        HkSpl = UnivariateSpline(kExitArr, HArr)

        Hk = HkSpl(k)

        Hend = MTE.H(back[-1][1:], params)

        Nk = (
            const
            + np.log(243.5 * 3**0.25 * np.sqrt(Hk))
            + np.log(np.sqrt(Hk / Hend))
        )

    except ValueError:
        print_warning(
            "warning error fitting spline for pivot "
            "exit scale: {}".format(k)
        )

        return np.nan

    return Nk


def _eps_eta_mij_hub(back, params, Nexit, MTE, func):
    if func == "eps":
        METHOD = MTE.Epsilon
    elif func == "eta":
        METHOD = MTE.Eta
    elif func == "mij":
        METHOD = MTE.massMatrix
    elif func == "hub":
        METHOD = MTE.H
    else:
        raise KeyError(func)

    def _eval(step, params):
        return METHOD(step[1 : 2 * MTE.nF() + 1], params)

    back_exit = splt.approx_row_closest(Nexit, back, 0)
    back_end = back[-1]

    return np.array(
        [_eval(back_exit, params), _eval(back_end, params)], dtype=float
    )


def get_epsilon_data(back, params, Nexit, MTE):
    return _eps_eta_mij_hub(back, params, Nexit, MTE, "eps")


def get_eta_data(back, params, Nexit, MTE):
    return _eps_eta_mij_hub(back, params, Nexit, MTE, "eta")


def get_mass_data(back, params, Nexit, MTE):
    mij = _eps_eta_mij_hub(back, params, Nexit, MTE, "mij")

    eigs_exit = np.sort(np.linalg.eigvals(mij[0]))
    eigs_end = np.sort(np.linalg.eigvals(mij[1]))

    return np.array([eigs_exit, eigs_end], dtype=float)


def evolve_mass_spectrum(
    back, params, MTE, scale_eigs=False, hess_approx=False, covariant=False
):
    """
    Evolves the mass spectrum as a function of efolds,
    normalizing the result by H

    :param back: background evolution
    :param params: model params
    :param MTE: PyTransport model
    :param scale_eigs: if True, rescales evolution to absolute sqrt
    :param hess_approx: if True, hess aprox to matrix
    :param covariant: if True, returns index down down expression
    :return: eigs(M^i_j) / H^2
    """

    nF = MTE.nF()
    out = np.zeros((len(back), 1 + nF), dtype=float)

    for idx, step in enumerate(back):
        mij = MTE.massMatrix(step[1:], params, hess_approx, covariant)
        out[idx][0] = step[0]
        out[idx][1:] = np.sort(np.linalg.eigvals(mij))

    if scale_eigs:
        out[::, 1:] = np.sqrt(np.abs(out[::, 1:]))

    return out


def evolveMasses(
    back, params, MTE, scale_eigs=False, hess_approx=False, covariant=False
):
    """Computes the mass matrix along a given background evolution,
    M^{I}_{J}
    """

    # Build empty array to populate with efold number + eigenvalues of
    # mass-matrix at time step
    eigs = np.empty((len(back), 1 + MTE.nF()))

    for idx, step in enumerate(back):
        # Unpack efolding and fields as np arrays
        N, fieldsdotfields = step[0], step[1 : 1 + 2 * MTE.nF()]

        # Compute mass matrix
        Mij = MTE.massMatrix(fieldsdotfields, params, hess_approx, covariant)

        # Compute and sort eigenvalues
        masses = np.linalg.eigvals(Mij)
        masses = np.sort(masses)

        # If scale eigs m^2 -> m
        if scale_eigs:
            masses = np.sign(masses) * np.sqrt(np.abs(masses))

        # Assign row values for mass-matrix
        eigs[idx] = np.concatenate((np.array([N]), masses))

    return eigs


def compute_spectral_index(
    MTE,
    back: np.ndarray,
    params: np.ndarray,
    tols: np.ndarray,
    sub_evo: int or float,
    Nexit=None,
    tmax=None,
):
    """
    Compute ns, running and As

    :param MTE:
    :param back:
    :param params:
    :param tols:
    :param sub_evo:
    :param Nexit:
    :param tmax:
    :return:
    """
    if Nexit is None:
        Nexit = compute_Nexit_for_matching(MTE, back, params)

    if Nexit is None:
        assert 0, "FIX"

    kexit_arr = np.zeros(5, dtype=float)

    for idx in range(5):
        _NEXIT = Nexit - 0.5 + 0.25 * idx

        if _NEXIT < 0:
            return -51

        _k = kexitN(_NEXIT, back, params, MTE)

        if np.isnan(_k) or np.isinf(_k):
            return -51

        kexit_arr[idx] = _k

    kexit = kexit_arr[2]

    # Build power spectra values
    p_zeta_arr = np.zeros(5, dtype=float)

    Nend = back.T[0][-1]

    scalar_amplitude = None

    # Iterate over k-scales
    for idx, k in enumerate(kexit_arr):
        # Build initial conditions for mode based on massless condition
        Nstart, ICs = find_massless_condition(MTE, sub_evo, k, back, params)

        if not isinstance(ICs, np.ndarray) and np.isnan(ICs):
            return -50

        _2pf = MTE.sigEvolve(
            np.array([Nstart, Nend]), k, ICs, params, tols, False, tmax, True
        )

        if isinstance(_2pf, int):
            return _2pf

        p_zeta_arr[idx] = _2pf.T[1][-1]

        if idx == 2:
            scalar_amplitude = p_zeta_arr[2] * (k**3)

    assert scalar_amplitude is not None, scalar_amplitude

    #  Build log arrays in k and Pzeta
    log_k_arr = np.log(kexit_arr / kexit)
    log_pzeta_arr = np.log(p_zeta_arr)

    # Build cubic spline from log values
    pzeta_spl = UnivariateSpline(log_k_arr, log_pzeta_arr, k=3)

    # Differentiate for ns & running
    ns_spline = pzeta_spl.derivative()
    alpha_spline = ns_spline.derivative()

    # Define functions to map exit scale to observables subject to kPicot
    def ns(k):
        return ns_spline(np.log(k / kexit)) + 4.0

    def alpha(k):
        return float(alpha_spline(np.log(k / kexit)))

    return np.array(
        [ns(kexit), alpha(kexit), scalar_amplitude], dtype=np.float64
    )


def compute_fnl(
    MTE,
    back: np.ndarray,
    params: np.ndarray,
    tols: np.ndarray,
    sub_evo: int or float,
    Nexit=None,
    tmax=None,
    alpha=None,
    beta=None,
    eq=False,
    fo=False,
    sq=False,
):
    assert sum([eq, fo, sq]) > 0 or (alpha is not None and beta is not None), [
        alpha,
        beta,
        eq,
        fo,
        sq,
    ]

    Nend = back.T[0][-1]

    if Nexit is None:
        Nexit = compute_Nexit_for_matching(MTE, back, params)

    if Nexit is None:
        assert 0, "FIX"

    kexit = kexitN(Nexit, back, params, MTE)

    # Build Fourier triangle via Fergusson Shellard convention;
    def _k1(alpha, beta, k):
        return k * (1.0 - beta) / 2.0

    def _k2(alpha, beta, k):
        return k * (1.0 + alpha + beta) / 4.0

    def _k3(alpha, beta, k):
        return k * (1.0 - alpha + beta) / 4.0

    kt = 3 * kexit

    # standard template options
    if eq:
        alpha, beta = 0, 1.0 / 3.0
    elif fo:
        alpha, beta = -0.5, 0.5
    elif sq:
        alpha, beta = 0.9, 0.01
    else:
        pass

    # Compute momentum triangle
    k1 = _k1(alpha, beta, kt)
    k2 = _k2(alpha, beta, kt)
    k3 = _k3(alpha, beta, kt)

    for k in [k1, k2, k3]:
        # if momenta has become infinite, or subject to numerical error
        if np.isinf(k) or np.isnan(k):
            return -61

    del k

    # Find larges scale and compute ICs
    kmin = np.min([k1, k2, k3])

    Nstart, ICs = find_massless_condition(MTE, sub_evo, kmin, back, params)

    if not isinstance(ICs, np.ndarray) and np.isnan(ICs):
        return -60

    # Compute three-point function up until end of background
    _3pf = MTE.alphaEvolve(
        np.array([Nstart, Nend]),
        k1,
        k2,
        k3,
        ICs,
        params,
        tols,
        True,
        tmax,
        True,
    )

    # If flag has been returned when computing 3pf, return flag data
    if isinstance(_3pf, int):
        return _3pf

    # Compute amplitude of fNL at end of inflation
    pzeta_1, pzeta_2, pzeta_3, bzeta = [_3pf.T[i][-1] for i in range(1, 5)]
    fnl = (
        (5.0 / 6.0)
        * bzeta
        / (pzeta_1 * pzeta_2 + pzeta_2 * pzeta_3 + pzeta_1 * pzeta_3)
    )

    return np.array([fnl], dtype=np.float64)
