import os
import numpy as np
import pickle as pk

import pytsa.pytrans_scripts as py_scripts
from pytsa.cache_tools import hash_alpha_beta
from pytsa.sampler.setup_sampler import LatinSampler, APrioriSampler, SamplerMethods


class ProtoAttributes:
    # Prototype for typing

    methods: SamplerMethods
    sampler: APrioriSampler or LatinSampler
    index: int
    cache: str
    task_dict: dict
    n_samples: int


error_dict = {
    0: "",
    -30: "short",
    -31: "fpace",
    -32: "int_back",
    -33: "timeout_back",
    -40: "eps",
    -41: "eta",
    -42: "mij",
    -50: "ics_2pf",
    -51: "kexit_2pf",
    -52: "int_2pf",
    -53: "timeout_2pf",
    -60: "ics_3pf_{alpha_beta}",
    -61: "kexit_3pf_{alpha_beta}",
    -62: "int_3pf_{alpha_beta}",
    -63: "timeout_3pf_{alpha_beta}"
}


class SampleCore(object):

    # Track background quantities for observables, and log errors

    def __init__(self, index: int, err_number: int, cache_loc, back: np.ndarray or None = None,
                 back_raw: np.ndarray or None = None,
                 Nexit: float or None = None,
                 Nexit_raw: float or None = None,
                 params: np.ndarray or None = None):

        self.index = index
        self.cache_loc = cache_loc
        self.path = os.path.join(self.cache_loc, "sample.%06d" % self.index)

        self.params = params

        self.err_n = [err_number]
        self.err_s = [error_dict[err_number]]

        self.results = {}

        self.back = back
        self.back_raw = back_raw

        self.Nexit = Nexit
        self.Nexit_raw = Nexit_raw

        self.adiabatic = None

        self.completed = set()  # TODO: we can actually just look at the results dict instead of tracking twice !

        self.cache()

    def get_last_status(self):
        return self.index, self.err_n[-1]

    @classmethod
    def accepted(cls, index, cache_loc):
        return cls(index, 0, cache_loc)

    @classmethod
    def inflation_too_short(cls, index, cache_loc):
        return cls(index, -30, cache_loc)

    @classmethod
    def field_space_violation(cls, index, cache_loc):
        return cls(index, -31, cache_loc)

    @classmethod
    def int_background(cls, index, cache_loc):
        return cls(index, -32, cache_loc)

    def task_eps(self, failed: bool, value):

        self.completed.add("eps")

        if failed:
            self._update_sample(-40, result_key="eps", result_value=value)
        else:
            self._update_sample(0, "eps", result_key="eps", result_value=value)

    def task_eta(self, failed: bool, value=None):

        self.completed.add("eta")

        if failed:
            self._update_sample(-41, result_key="eta", result_value=value)
        else:
            self._update_sample(0, result_key="eta", result_value=value)

    def task_mij(self, failed: bool, adiabatic: bool, value=None):

        self.completed.add("mij")

        self.adiabatic = adiabatic

        if failed:
            self._update_sample(-42, result_key="mij", result_value=value)
        else:
            self._update_sample(0, result_key="mij", result_value=value)

    def task_2pf(self, code: int, value=None):

        self.completed.add("2pf")

        if -60 < code < -50:
            self._update_sample(code, result_key="2pf", result_value=value)
        else:
            assert code == 0, code
            self._update_sample(0, result_key="2pf", result_value=value)

    def task_3pf(self, code: int, alpha_beta, value=None):

        self.completed.add(alpha_beta)

        if -70 < code < -60:
            self._update_sample(code, alpha_beta, result_key=alpha_beta, result_value=value)
        else:
            self._update_sample(0, result_key=alpha_beta, result_value=value)

    def _update_sample(self, error_code: int, alpha_beta=None, result_value=None, result_key=None):

        self.err_n.append(error_code)

        es = error_dict[error_code]

        self.err_s.append(es if alpha_beta is None else es.format(alpha_beta=alpha_beta))

        assert result_key is not None

        self.results[result_key] = result_value

        self.cache()

    def cache(self):

        p = self.path

        if not os.path.exists(p):
            with open(p, "wb") as f:
                pk.dump(self, f)

        else:
            os.remove(p)
            self.cache()


def extract_core(data: ProtoAttributes) -> SampleCore:
    cache = data.cache

    sample_path = os.path.join(cache, "sample.%06d" % data.index)

    with open(sample_path, "rb") as f:
        sample_core: SampleCore = pk.load(f)

    return sample_core


def compute_background(data: ProtoAttributes):
    try:
        loaded = extract_core(data)
        return loaded.index, loaded.err_n[0]
    except FileNotFoundError:
        pass

    model = data.methods.model
    index = data.index
    cache = data.cache

    tols = np.array([*data.methods.tols], dtype=float)
    Nmin = data.methods.N_min
    ics, params = data.sampler.get_sample(index)

    Nend = model.findEndOfInflation(ics, params, tols)

    N_evo = np.linspace(0, Nend, int(100 * Nend))  # Dense N evo
    background = model.backEvolve(N_evo, ics, params, tols, False, -1, True)  # Raw background

    # Dictionary with keys for field indices, and arrays for end inflation
    terminate_with_fields = data.methods.fieldspace_terminate

    if len(terminate_with_fields) > 0:

        _idx = np.inf

        for idx, step in enumerate(background):

            for fidx, bounds in terminate_with_fields.items():

                if bounds[0] <= step[1 + fidx] <= bounds[1]:
                    _idx = np.min([_idx, idx])

        if np.isinf(_idx):
            return SampleCore.field_space_violation(index, cache).get_last_status()

        background = background[:_idx + 1]
        N_evo = background.T[0]
        Nend = N_evo[-1]

    reject_with_fields = data.methods.fieldspace_reject

    if len(reject_with_fields) > 0:

        for idx, step in enumerate(background):

            for fidx, bounds in reject_with_fields.items():

                if bounds[0] <= step[1 + fidx] <= bounds[1]:
                    SampleCore.field_space_violation(index, cache).get_last_status()

    if Nend < Nmin:
        return SampleCore.inflation_too_short(index, cache).get_last_status()

    if isinstance(Nend, tuple):
        return SampleCore.int_background(index, cache).get_last_status()

    if isinstance(background, tuple):
        return SampleCore.int_background(index, cache).get_last_status()

    back_adj = py_scripts.adjust_back(background)  # Adjusted background

    nexit_adj = py_scripts.compute_Nexit_for_matching(model, back_adj, params)
    nexit = py_scripts.compute_Nexit_for_matching(model, background, params)

    return SampleCore(index, 0, cache, back_adj, background, nexit_adj, nexit, params).get_last_status()


from numpy.lib import format as npf


def compute_epsilon(data: ProtoAttributes):
    sample_core = extract_core(data)

    if "eps" not in sample_core.completed:
        eps = py_scripts.get_epsilon_data(sample_core.back, sample_core.params, sample_core.Nexit, data.methods.model)
        sample_core.task_eps(False, eps)  # TODO: error


def compute_eta(data: ProtoAttributes):
    sample_core = extract_core(data)

    if "eta" not in sample_core.completed:
        eta = py_scripts.get_eta_data(sample_core.back, sample_core.params, sample_core.Nexit, data.methods.model)
        sample_core.task_eta(False, value=eta)  # TODO: error


def compute_mij(data: ProtoAttributes):
    sample_core = extract_core(data)

    if "mij" not in sample_core.completed:

        masses_exit, masses_end = py_scripts.get_mass_data(sample_core.back, sample_core.params, sample_core.Nexit,
                                                           data.methods.model)
        N_adi = data.methods.N_adiabitc

        back = sample_core.back

        NEND = back.T[0][-1]

        N_adi_check = np.array([N for N in back.T[0] if N > NEND - N_adi])

        # Check for adiabatic limit: Criteria here requires that m^2 >= 2 H^2 over the last efolds of inflation
        # for all but one tachyonic mode.
        # precise nember defined as N_adiabatic in setup

        adi = True
        for _ in N_adi_check:
            masses = py_scripts.get_mass_data(sample_core.back, sample_core.params, _, data.methods.model)

            n_tachyon = (masses < 0).sum()

            if n_tachyon != 1:
                adi = False
                break

            rest_heavy = np.all(np.sort(masses)[1:] > 2)

            if not rest_heavy:
                adi = False
                break

        sample_core.task_mij(False, adi, value=[masses_exit, masses_end])  # TODO: error


def compute_obs(data: ProtoAttributes):
    task_dict = data.task_dict

    if task_dict['task_2pt'] is True:
        compute_2pf(data)

    if task_dict['task_3pt_eq'] is True:
        compute_3pf(data, eq=True)

    if task_dict['task_3pt_fo'] is True:
        compute_3pf(data, fo=True)

    if task_dict['task_3pt_sq'] is True:
        compute_3pf(data, sq=True)

    if len(task_dict['alpha']) > 0:
        for alpha, beta in zip(task_dict['alpha'], ['beta']):
            compute_3pf(data, alpha=alpha, beta=beta)

    return 0


def compute_2pf(data: ProtoAttributes):
    sample_core = extract_core(data)

    if "2pf" not in sample_core.completed:
        result = py_scripts.compute_spectral_index(
            data.methods.model,
            sample_core.back,
            sample_core.params,
            np.array([*data.methods.tols], dtype=float),
            data.methods.N_sub_evo,
            Nexit=sample_core.Nexit,
            tmax=600  # 10 minute maximum, this should be plenty
        )

        sample_core.task_2pf(0, value=result)  # TODO, error


def compute_3pf(data: ProtoAttributes, **task_kwargs):
    sample_core = extract_core(data)

    if task_kwargs["eq"]:
        alpha_beta = "eq"
    elif task_kwargs["sq"]:
        alpha_beta = "sq"
    elif task_kwargs["fo"]:
        alpha_beta = "fo"
    else:
        alpha = task_kwargs['alpha']
        beta = task_kwargs['beta']
        alpha_beta = hash_alpha_beta(alpha, beta)

    if alpha_beta not in sample_core.completed:
        result = py_scripts.compute_fnl(
            data.methods.model,
            sample_core.back,
            sample_core.params,
            np.array([*data.methods.tols], dtype=float),
            data.methods.N_sub_evo,
            Nexit=sample_core.Nexit,
            tmax=600,  # 10 minute maximum, this should be plenty
            **task_kwargs
        )

        sample_core.task_3pf(0, alpha_beta, value=result)

#
# def buildICPs(modelnumber, rerun_model=False):
#     """ Builds initial conditions & parameters, or retrieves them if in rerun mode"""
#
#     if rerun_model is False:
#         # Generate new sample with parameter & field priori
#         n, fvals, vvals, pvals = genSample(modelnumber)
#
#     else:
#         # Load an existing set of parameters & model priori
#         model_path = os.path.join(saveloc, "samples", "{}.sample".format(modelnumber))
#         model_file = open(model_path, "rb")
#
#         with model_file as mf:
#             model = pk.load(mf)
#             fvals = model.fields
#             vvals = model.velocities
#             pvals = model.parameters
#
#             return np.concatenate(fvals, vvals), pvals
#
#     # If SR is in vvals, compute initial velocity value via slowroll equation
#     if "SR" in vvals:
#
#         # Compute intial velocities from approximation
#         vvalsSR = PyT.dotfieldsSR(fvals, pvals)
#
#         # Populate empty array with relevant values
#         vvals_ = np.zeros(PyT.nF())
#
#         # For each velocity def.
#         for jj in range(PyT.nF()):
#
#             #  If def. is SR, assign slow-roll velocity to float array
#             if vvals[jj] == "SR":
#                 vvals_[jj] = vvalsSR[jj]
#             else:
#                 # Otherwise force float-type to other array elements
#                 vvals_[jj] = float(vvals[jj])
#
#         # Redefine vvals
#         vvals = vvals_
#
#     # Concatenate field and velocity values into initial conditions array
#     return np.concatenate((fvals, vvals)), pvals
#
#
# def Initialize(modelnumber, rerun_model=False, initial=None, pvals=None, returnSample=False):
#     tmax_bg = int(os.environ['tmax_bg'])
#
#     if initial is None and pvals is None:
#         initial, pvals = buildICPs(modelnumber, rerun_model)
#
#     """ We have the following flags in place:
#
#     Timeout flags:
#     -10 : findEndOfInflation
#     -11 : backEvolve
#     -12 : sigEvolve
#     -13 : alphaEvolve
#
#     Integration error flags:
#     -20 : findEndOfInflation
#     -21 : beckEvolve
#     -22 : sigEvolve
#     -23 : alphaEvolve
#
#     Misc. PyT:
#     -30 : Nend < Nmin
#     -31 : k not found
#     -32 : No ICsBM (2pf)
#     -33 : No ICsBM (3pf)
#     -34 : Model violation
#     -35 : Eternal inflation
#
#     """
#
#     acceptIfAnyPath = os.path.join(pathLocalData, "acceptIfAny.localdata")
#     acceptIfAllPath = os.path.join(pathLocalData, "acceptIfAll.localdata")
#
#     rejectIfAnyPath = os.path.join(pathLocalData, "rejectIfAny.localdata")
#     rejectIfAllPath = os.path.join(pathLocalData, "rejectIfAll.localdata")
#
#     # We define non-canonical conditions to end inflation, and those which violate inflation
#     endAnyNC = os.path.exists(acceptIfAnyPath)  # NC end
#     endAllNC = os.path.exists(acceptIfAllPath)
#
#     exitAnyNC = os.path.exists(rejectIfAnyPath)  # NC violate
#     exitAllNC = os.path.exists(rejectIfAllPath)
#
#     canonicalEnd = not (endAnyNC or endAllNC)
#     canonicalExit = not (exitAnyNC or exitAllNC)
#
#     # if canonical end or not, suitable to start with cheap epsilon search
#     Nepsilon = PyT.findEndOfInflation(initial, pvals, tols, 0.0, 10000, tmax_bg, True)
#
#     # First we check models that end with SR violation
#     if canonicalEnd:
#
#         # Check for integration / time out error
#         if type(Nepsilon) is tuple:
#             return Nepsilon
#
#         # Check for enough inflation
#         elif Nepsilon < minN:
#             return (-30, Nepsilon)
#
#         # Looks good, go compute background
#         else:
#             back = PyT.backEvolve(
#                 np.linspace(0, Nepsilon, int(3 * Nepsilon)), initial, pvals, tols, True, tmax_bg, True)
#
#         Nend = Nepsilon
#
#     # For models that don't end with SR violation
#     else:
#
#         # Check for integration / time out error
#         if type(Nepsilon) is tuple:
#             # Instead, draw back farthest efold reached by 10% and try and compute fid. background
#             Nepsilon = 0.9 * Nepsilon[1]
#
#         # Try and obtain fiducial background up until epsilon = 1
#         backFid = PyT.backEvolve(
#             np.linspace(0, Nepsilon, 3 * (int(Nepsilon) + 1)), initial, pvals, tols, 0, tmax_bg, True)
#
#         """ We will now attempt to extend the background by up to 100 efolds passed the epsilon definition
#             Note that getting further than this is highly unlikely, since integration becomes hard
#             Note that we increase desired tols by an order of magnitude to try and achieve this """
#
#         backExt = PyT.backEvolve(
#             np.linspace(backFid[-1][0], backFid[-1][0] + 100, 300), backFid[-1][1:], pvals, tols * 1e-1, 0, tmax_bg,
#             True)
#
#         # Combine arrays, omiting the first row of the extension since this will match the last of fid.
#         back = np.vstack((backFid, backExt[1:]))
#
#         # We now will search through the background and look for Nend condition
#         endIndex = len(back)
#         foundEnd = False
#
#         # Load conditions
#         if endAnyNC:
#             with open(acceptIfAnyPath, "rb") as f: acceptIfAny = pk.load(f)
#             acceptIfAnyIdx = acceptIfAny['idx']
#             acceptIfAnyMin = acceptIfAny['min']
#             acceptIfAnyMax = acceptIfAny['max']
#
#         if endAllNC:
#             with open(acceptIfAllPath, "rb") as f: acceptIfAll = pk.load(f)
#             acceptIfAllIdx = acceptIfAll['idx']
#             acceptIfAllMin = acceptIfAll['min']
#             acceptIfAllMax = acceptIfAll['max']
#
#         for idx in range(len(back)):
#
#             row = back[idx]
#
#             if endAnyNC and np.any(
#                     np.logical_and(row[acceptIfAnyIdx] > acceptIfAnyMin, row[acceptIfAnyIdx] < acceptIfAnyMax)):
#                 endIndex = idx
#                 Nend = row[0]
#                 foundEnd = True
#                 break
#
#             if endAllNC and np.all(
#                     np.logical_and(row[acceptIfAllIdx] > acceptIfAllMin, row[acceptIfAllIdx] < acceptIfAllMax)):
#                 endIndex = idx
#                 Nend = row[0]
#                 foundEnd = True
#                 break
#
#         endIndex = min(endIndex + 1, len(back))
#
#         back = back[:endIndex]
#
#     if canonicalExit is False:
#
#         if exitAnyNC:
#             with open(rejectIfAnyPath, "rb") as f: rejectIfAnyPath = pk.load(f)
#             rejectIfAnyIdx = rejectIfAnyPath['idx']
#             rejectIfAnyMin = rejectIfAnyPath['min']
#             rejectIfAnyMax = rejectIfAnyPath['max']
#
#         if exitAllNC:
#             with open(rejectIfAllPath, "rb") as f: rejectIfAllPath = pk.load(f)
#             rejectIfAllIdx = rejectIfAllPath['idx']
#             rejectIfAllMin = rejectIfAllPath['min']
#             rejectIfAllMax = rejectIfAllPath['max']
#
#         for row in back:
#
#             if exitAnyNC and np.any(
#                     np.logical_and(row[rejectIfAnyIdx] > rejectIfAnyMin, row[rejectIfAnyIdx] < rejectIfAnyMax)):
#                 Nend = row[0]
#
#                 return -34, Nend
#
#             if exitAllNC and np.all(
#                     np.logical_and(row[rejectIfAllIdx] > rejectIfAllMin, row[rejectIfAllIdx] < rejectIfAllMax)):
#                 Nend = row[0]
#
#                 return -34, Nend
#
#     # Now perform the final check that the background
#     if canonicalEnd:
#         if Nend < minN:
#             return -30, Nend
#     else:
#
#         Nend = back[-1][0]
#
#         if Nend > minN and not foundEnd:
#             return -35, Nend
#
#         if Nend < minN and foundEnd:
#             return -30, Nend
#
#         if Nend < minN and not foundEnd:
#             return -21, Nend
#
#     # Track *actual* initial conditions prior to repositioning
#     initial_0 = initial
#
#     # Rescale background to start at Nend - minN to avoid exp. large k
#     back_adj = PyS.rescaleBack(back, Nr=minN + 10)
#
#     # Attempt computation of momenta at horizon crossing, the value used will be in the ballpark of requirements
#     # for 2pf and 3pf tasks. Hence, we discard the realization at this point if it's prospects look bad.
#
#     Npiv = PyS.matchKExitN(back_adj, pvals, PyT, 0.002)
#
#     if np.isnan(Npiv):
#         return -31
#
#     kExit = PyS.kexitN(Nend - Npiv, back_adj, pvals, PyT)
#
#     # Asses success by data type of momenta result
#     if type(kExit) not in [float, np.float, np.float32, np.float64] or np.isnan(kExit) or np.isinf(kExit):
#         return -31
#
#     # We now test for an adiabatic limit: We begin by assuming this is True, then test for violating conditions
#     adiabatic = True
#
#     # Define number of efolds from end of inflation which should be adiabatic TODO: Wrap this into config. file
#
#     # Find first efoldign to perform eigenvalue test
#     try:
#         adiabaticN_start = np.where(back.T[0] >= back_adj[-1][0] - adiabaticN)[0][0]
#     except IndexError:
#         return -21
#
#     # Compute mass-matrix evolution from this point
#     Mij_end = PyS.evolveMasses(back[adiabaticN_start:], pvals, PyT, scale_eigs=False, hess_approx=False,
#                                covariant=False)
#
#     # For each mass-matrix evolution step
#     for item in Mij_end:
#
#         # Determine if lightest mass is tachyonic
#         tachyon = item[1] < 0
#
#         # Determine number of hubble scale heavier masses
#         rest_large = sum([eig >= 1 for eig in item[2:]]) == PyT.nF() - 1
#
#         # If there is no tachyon or the remaining masses aren't large for the portion of evolution
#         if not (tachyon and rest_large):
#             # Change adiabatic property to False, and break out of search window
#             adiabatic = False
#
#             break
#
#     # Record all sample data
#     new_sample(modelnumber, initial_0[:PyT.nF()], initial_0[PyT.nF():], pvals, back_adj, adiabatic,
#                Nend, kExit, pathSamples, rerun_model)
#
#     print
#     "-- Generated sample: %05d" % modelnumber
#
#     assert type(modelnumber) == int, modelnumber
#
#     # Return model number, acts as flag to confirm successful sample
#     return modelnumber
#
#
# def updateStats(modelnumber, flagKey=None, timeEnd=None):
#     """ Rewrites stats object during demand sample routine """
#
#     # Get directory for sample stats log
#     sample_path = os.path.join(pathStats, "bg", "{}.bg".format(modelnumber))
#
#     if os.path.exists(sample_path):
#         with open(sample_path, "r") as f:
#             flagDict = pk.load(f)
#         pathExists = True
#
#     else:
#         pathExists = False
#
#         # Flag defs.
#         timeoutFlags = [
#             -10, -11, -12, -13  # end of inflation, background, 2pf, 3pf
#         ]
#
#         integratorFlags = [
#             -20, -21, -22, -23  # end of inflation, background, 2pf, 3pf
#         ]
#
#         samplerFlags = [
#             -30, -31, -32, -33, -34, -35, -36  # N < Nmin, k not found, no ICs 2pf, no ICs 3pf, model violation, eternal
#         ]
#
#         allFlags = timeoutFlags + integratorFlags + samplerFlags
#         flagDict = {str(f): 0 for f in allFlags}
#
#     if flagKey in flagDict:
#         flagDict[flagKey] += 1
#     elif flagKey is None and timeEnd is not None:
#         flagDict['time'] = timeEnd
#     else:
#         raise KeyError, flagKey
#
#     if pathExists:
#         os.remove(sample_path)
#
#     with open(sample_path, "w") as f:
#         pk.dump(flagDict, f)
#
#
# def DemandSample(modelnumber, rerun=False):
#     """ Repeat initialization until successful sample is found """
#
#     # # Get directory for sample stats log
#     # sample_path = os.path.join(pathStats, "bg", "{}.bg".format(modelnumber))
#     #
#     # Flag defs.
#     timeoutFlags = [
#         -10, -11, -12, -13  # end of inflation, background, 2pf, 3pf
#     ]
#
#     integratorFlags = [
#         -20, -21, -22, -23  # end of inflation, background, 2pf, 3pf
#     ]
#
#     samplerFlags = [
#         -30, -31, -32, -33, -34, -35, -36
#         # N < Nmin, k not found, no ICs 2pf, no ICs 3pf, model violation, eternal, mij eig fail
#     ]
#
#     allFlags = timeoutFlags + integratorFlags + samplerFlags
#     # flagDict = {str(f): 0 for f in allFlags}
#
#     # Start timer
#     tstart = time.clock()
#
#     ii = -1
#
#     # If successful sample generation, the model number is returned
#     while ii != modelnumber:
#
#         ii = Initialize(modelnumber, rerun_model=rerun)
#
#         # If fail flag received: log statistic
#         if type(ii) is tuple:
#             flagKey = str(ii[0])
#             timeEnd = None
#             # flagDict[flagKey] += 1
#
#         elif type(ii) is int and ii in allFlags:
#             flagKey = str(ii)
#             timeEnd = None
#             # flagDict[flagKey] += 1
#
#         # If model number, sum number of iterations
#         elif ii == modelnumber:
#             flagKey = None
#             timeEnd = time.clock() - tstart
#
#         else:
#             raise ValueError, modelnumber
#
#         updateStats(modelnumber, flagKey=flagKey, timeEnd=timeEnd)
#
#
# def DemandSample_rerun(modelnumber):
#     """
#     Repeat initialization until successful sample is found : using rerun samples
#     """
#     DemandSample(modelnumber, rerun=True)
#
#
# def Mij(modelnumber):
#     """
#     Computes the eigenvalues of the mass-square-matrix at horizon crossing
#     """
#
#     # Unload sample data
#     path = os.path.join(pathSamples, "{}.sample".format(modelnumber))
#     assert os.path.exists(path), "Unable to locate sample location: {}".format(path)
#     s = open(path, "rb")
#     with s:
#         sample = pk.load(s)
#
#     # Cross-check model number identifier
#     assert sample.modelnumber == modelnumber, "Pickled model data does not match! {m1} != {m2}".format(
#         m1=modelnumber, m2=sample.modelnumber
#     )
#
#     # Unpack background data
#     back = sample.background
#     Nend = sample.Nend
#     assert Nend is not None, "Nend = None -> Sample construction failure!"
#     pvals = np.asarray(sample.parameters)
#
#     # Transport background for col evolution rather than row evolution
#     backT = back.T
#
#     # Get number of fields
#     nF = PyT.nF()
#
#     # Begin timer
#     ti = time.clock()
#
#     # Define Nstar, i.e. exit time of interest
#     Nstar = Nend - Nexit
#
#     # Unpack efold evolution
#     Nevo = backT[0]
#
#     # If Nstar falls exactly in the background, simply load this background step
#     if Nstar in Nevo:
#         Nstar_idx = [Nevo.index(Nstar)]
#         back_star = back[Nstar_idx]
#
#     # Otherwise we build a spline about Nstar
#     else:
#
#         # Get background index values +/- 1 efold around the prescribed exit time
#         Nfit = 0
#         dN = 0.5
#         while Nfit < 3:
#             dN += 0.5
#             idx_above = np.where(np.logical_and(Nevo > Nstar, Nevo < Nstar + dN))[0]
#             idx_below = np.where(np.logical_and(Nevo < Nstar, Nevo > Nstar - dN))[0]
#
#             Nstar_idx = list(np.concatenate([idx_below, idx_above]))
#             Nfit = len(Nstar_idx)
#
#         # Define efold range to "smooth" over
#         Nsmooth = [Nevo[idx] for idx in Nstar_idx]
#
#         # We check for repeated efold numbers as this will throw the interpolation scheme
#         repLog = []
#         for ii in range(len(Nsmooth) - 1):
#             if Nsmooth[ii] == Nsmooth[ii + 1]:
#                 repLog.append(ii + 1)
#         for index in sorted(repLog, reverse=True):
#             del Nstar_idx[index]
#
#         Nsmooth = [Nevo[idx] for idx in Nstar_idx]
#
#         # BY default, we choose a cubic spline, but if this is not possible, we reduce the order
#         if len(Nstar_idx) <= 3:
#             k = len(Nstar_idx) - 1
#         else:
#             k = 3
#
#         # Construct splines for the field evolution over Nsmooth
#         fsplines = [
#             UnivariateSpline(
#                 Nsmooth, [backT[i][idx] for idx in Nstar_idx], k=k) for i in range(1, 1 + nF)
#         ]
#
#         # Construct splines for the velocity evolution over Nsmooth
#         vsplines = [
#             UnivariateSpline(
#                 Nsmooth, [backT[i][idx] for idx in Nstar_idx], k=k) for i in range(1 + nF, 1 + 2 * nF)
#         ]
#
#         # Evaluate the background at N_star = Nend - Nexit
#         back_star = np.concatenate([np.array([Nstar]), np.array([s(Nstar) for s in fsplines + vsplines])])
#
#     try:
#         # Compute mass matrix eigenvalues at the exit time
#         Mij = PyS.evolveMasses(
#             np.array([back_star]), pvals, PyT, scale_eigs=False, hess_approx=False, covariant=False)[0][1:]
#
#     except np.linalg.LinAlgError:
#
#         # Build dictionary of masses
#         Mij_eig = {}
#         for i in range(PyT.nF()):
#             Mij_eig['m{}'.format(i)] = None
#
#         Mij_eig['T_masses'] = None
#
#         sample.update_observables(Mij_eig)
#
#         return {"mn": modelnumber, "flag": -36, "ext": "mij"}
#
#     # Build dictionary of masses
#     Mij_eig = {}
#     for i in range(len(Mij)):
#         Mij_eig['m{}'.format(i)] = Mij[i]
#
#     # End timer
#     Mij_eig['T_masses'] = time.clock() - ti
#
#     # Update sample observables
#     sample.update_observables(Mij_eig)
#
#     return 0
#
#
# def SpectralIndex(modelnumber):
#     # Get timeout
#     tmax_2pf = int(os.environ['tmax_2pf'])
#
#     # Unload sample data
#     path = os.path.join(pathSamples, "{}.sample".format(modelnumber))
#     assert os.path.exists(path), "Unable to locate sample location: {}".format(path)
#     s = open(path, "rb")
#     with s:
#         sample = pk.load(s)
#
#     # Cross-check model number identifier
#     assert sample.modelnumber == modelnumber, "Pickled model data does not match! {m1} != {m2}".format(
#         m1=modelnumber, m2=sample.modelnumber
#     )
#
#     # def spectralIndex(back, pvals, Nexit, tols, subevo, MTE, kPivot=None, returnRunning=True, errorReturn=False, tmax=None):
#     ti = time.clock()
#     ns_alpha = PyS.spectralIndex(
#         sample.background, sample.parameters, Nexit, tols, subevo, PyT, errorReturn=True, tmax=tmax_2pf)
#
#     # Check for failed ns calculation
#     if isinstance(ns_alpha, tuple):
#         if ns_alpha[0] is ValueError:
#
#             errKey = ns_alpha[1]
#
#             # Define error dictionary to return
#             if errKey == "k":
#                 retDict = {"mn": modelnumber, "flag": -31, "ext": "2pf"}
#             elif errKey == "ics":
#                 retDict = {"mn": modelnumber, "flag": -32, "ext": "2pf"}
#             elif errKey == "2pf":
#                 retDict = {"mn": modelnumber, "flag": ns_alpha[2][0], "ext": "2pf"}
#             else:
#                 raise KeyError, "Unknown error: {}".format(errKey)
#
#             # Update sample with null result
#             sample.update_observables({'ns': None, 'alpha': None, 'T_2pf': None})
#
#             return retDict
#
#     # Unpack 2pt observables
#     ns, alpha = ns_alpha
#
#     # End calculation timer and build results dictionaryobs_
#     twoPt_dict = {'ns': ns, 'alpha': alpha, 'T_2pf': time.clock() - ti}
#
#     # Update sample observables
#     sample.update_observables(twoPt_dict)
#
#     return 0
#
#
# def fNL(modelnumber, configName):
#     tmax_3pf = int(os.environ['tmax_3pf'])
#
#     name = None
#     # Gather configuration data from configName
#     for d in fNLDict:
#         if d['name'] == configName:
#             name = configName
#             alpha = d['alpha']
#             beta = d['beta']
#             break
#
#     # Check that fNL definition has been found
#     assert name is not None, "Unable to locate fNL configuration in localdata: {}".format(name)
#
#     # Unload sample data
#     path = os.path.join(pathSamples, "{}.sample".format(modelnumber))
#     assert os.path.exists(path), "Unable to locate sample location: {}".format(path)
#     s = open(path, "rb")
#     with s:
#         sample = pk.load(s)
#
#     # Cross-check model number identifier
#     assert sample.modelnumber == modelnumber, "Pickled model data does not match! {m1} != {m2}".format(
#         m1=modelnumber, m2=sample.modelnumber
#     )
#     # Begin timer
#     ti = time.clock()
#
#     # Compute fNL
#     fNL = PyS.fNL(
#         sample.background, sample.parameters, Nexit, tols, subevo, PyT, alpha=alpha, beta=beta, stdConfig=name,
#         errorReturn=True, tmax=tmax_3pf)
#
#     # Check for failed fNL calculation calculation
#     if isinstance(fNL, tuple):
#         if fNL[0] is ValueError:
#
#             errKey = fNL[1]
#
#             # Define error dictionary to return
#             if errKey == "k":
#                 retDict = {"mn": modelnumber, "flag": -31, "ext": name}
#             elif errKey == "ics":
#                 retDict = {"mn": modelnumber, "flag": -33, "ext": name}
#             elif errKey == "3pf":
#                 retDict = {"mn": modelnumber, "flag": fNL[2][0], "ext": name}
#             else:
#                 raise KeyError, "Unknown error: {}".format(errKey)
#
#             # Update sample with null result
#             sample.update_observables({'{}'.format(name): None, 'T_{}'.format(configName): None})
#
#             return retDict
#
#     # Build result dictionary
#     threePt_Dict = {'{}'.format(name): fNL, 'T_{}'.format(configName): time.clock() - ti}
#
#     # Update sample observables
#     sample.update_observables(threePt_Dict)
#
#     return 0
#
#
# def computations(mn_calc):
#     """ Simple handler for computations called from sampler """
#
#     # Unpack model number and calculation type
#     modelnumber, calculation = mn_calc
#
#     # Prevent computation fails from hanging via catch_warnings
#     with warnings.catch_warnings(record=True):
#
#         if calculation == "masses":
#             print
#             "-- Start MAb: model %06d" % modelnumber
#             r = Mij(modelnumber)
#             status = "" if r == 0 else ", failed"
#             print
#             "--   End MAb: model %06d" % modelnumber + status
#
#         elif calculation == "2pf":
#             print
#             "-- Start 2pf: model %06d" % modelnumber
#             r = SpectralIndex(modelnumber)
#             status = "" if r == 0 else ", failed"
#             print
#             "--   End 2pf: model %06d" % modelnumber + status
#
#         elif calculation not in ["masses", "2pf"]:
#             print
#             "-- Start fNL: model %06d, config {}".format(calculation) % modelnumber
#             r = fNL(modelnumber, calculation)
#             status = "" if r == 0 else ", failed"
#             print
#             "--   End fNL: model %06d, config {}".format(calculation) % modelnumber + status
#
#         else:
#             raise ValueError, "Undefined calculation: {}".format(calculation)
#
#         if r != 0:
#
#             # Determine subdirectory for error records
#             # For now, we will store masses fails in 2pf dir to avoid restructuring
#             if calculation == "masses":
#                 subdir = "mij"
#             elif calculation == "2pf":
#                 subdir = "2pf"
#             else:
#                 subdir = "3pf"
#
#             # Get / remove extension key from dictionary, used to discriminate 3pf varieties
#             pkEXT = r.pop('ext')
#
#             # Dump dictionary object
#             pkPath = os.path.join(pathStats, subdir, "{}.{}".format(modelnumber, pkEXT))
#
#             pkFile = open(pkPath, "wb")
#             with pkFile:
#                 pk.dump(r, pkFile)
