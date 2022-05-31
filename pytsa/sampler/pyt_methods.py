import os
import numpy as np
import pickle as pk
import dill

from pytsa import pytrans_scripts as py_scripts
from ..cache_tools import hash_alpha_beta
from .setup_sampler import LatinSampler, APrioriSampler, SamplerMethods
from pytsa.sampler.mpi_helpers import rank

try:
    env_verbose = os.environ['pytsa_VERBOSE']
    verbose = bool(int(env_verbose))
except KeyError:
    env_verbose = None
    verbose = False

try:
    env_sampler = os.environ['pytsa_SAMPLER']
    with open(env_sampler, "rb") as f:
        sampler: LatinSampler or APrioriSampler = dill.load(f)

    root_dir = os.path.split(env_sampler)[0]
    samples_dir = os.path.join(root_dir, "samples_core")

    model = sampler.hyp_pars.model
    n_fields = sampler.hyp_pars.nF
    n_params = sampler.hyp_pars.nP

    fspace_reject = sampler.hyp_pars.fieldspace_reject
    fspace_terminate = sampler.hyp_pars.fieldspace_terminate

    Nadi = sampler.hyp_pars.N_adiabatic
    Nmin = sampler.hyp_pars.N_min
    Nsub = sampler.hyp_pars.N_sub_evo
    tols = np.array([*sampler.hyp_pars.tols], dtype=float)
    step_density = sampler.hyp_pars.step_density
    tmag_bg = sampler.hyp_pars.tmax_bg
    tmax_2pf = sampler.hyp_pars.tmax_2pf
    tmax_3pf = sampler.hyp_pars.tmax_3pf

except KeyError:
    env_sampler = None
    sampler = None
    print("-- Unable to load sampler data, read-only mode")

error_dict = {
    -00: "",
    -30: "short",
    -31: "fpace",
    -32: "int_back",
    -33: "timeout_back",
    -34: "nexit",
    -35: "eternal",
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

_BAD = np.nan


class SampleCore(object):

    # Track background quantities for observables, and log errors

    def __init__(self, index: int, back_status: int, back: np.ndarray or None = None,
                 back_raw: np.ndarray or None = None,
                 Nexit: float or None = None,
                 Nexit_raw: float or None = None,
                 params: np.ndarray or None = None):

        self.index = index
        self.path = os.path.join(samples_dir, "sample.%06d" % index)

        self.params = params

        self.err_n = [back_status]
        self.err_s = [error_dict[back_status]]

        self.results = {}

        self.back = back
        self.back_raw = back_raw

        self.Nexit = Nexit
        self.Nexit_raw = Nexit_raw

        self.adiabatic = None

        self.cache()

    def get_last_status(self):
        return self.index, self.err_n[-1]

    @classmethod
    def accepted(cls, index):
        return cls(index, 0)

    @classmethod
    def inflation_too_short(cls, index):
        return cls(index, -30)

    @classmethod
    def field_space_violation(cls, index):
        return cls(index, -31)

    @classmethod
    def int_background(cls, index):
        return cls(index, -32)

    @classmethod
    def timeout_background(cls, index):
        return cls(index, -33)

    @classmethod
    def nexit_background(cls, index):
        return cls(index, -34)

    def task_eps(self, failed: bool, value):

        if failed:
            self._update_sample(-40, result_key="eps", result_value=value)
        else:
            self._update_sample(0, "eps", result_key="eps", result_value=value)

    def task_eta(self, failed: bool, value):

        if failed:
            self._update_sample(-41, result_key="eta", result_value=value)
        else:
            self._update_sample(0, result_key="eta", result_value=value)

    def task_mij(self, failed: bool, adiabatic: bool, value):

        self.adiabatic = adiabatic

        if failed:
            self._update_sample(-42, result_key="mij", result_value=value)
        else:
            self._update_sample(0, result_key="mij", result_value=value)

    def task_2pf(self, code: int, value):

        if -60 <= code <= -50:
            self._update_sample(code, result_key="2pf", result_value=value)
        else:
            assert code == 0, [code, value]
            self._update_sample(0, result_key="2pf", result_value=value)

    def task_3pf(self, code: int, alpha_beta, value):

        if -70 <= code <= -60:
            self._update_sample(code, alpha_beta, result_key=alpha_beta, result_value=value)
        else:
            assert code == 0, [code, value]
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

    def get_row_data(self, *obs_keys):

        if self.err_n[0] != 0:  # Did not find background
            return _BAD

        results = np.array([], dtype=np.float64)

        for k in ['eps', 'eta']:
            results = np.concatenate((results, self.results[k]))

        for masses in self.results['mij']:
            results = np.concatenate((results, masses))

        for obs in obs_keys:
            results = np.concatenate((results, self.results[obs]))

        return results


def extract_core(index: int) -> SampleCore:

    sample_path = os.path.join(samples_dir, "sample.%06d" % index)

    with open(sample_path, "rb") as f:
        sample_core: SampleCore = pk.load(f)

    return sample_core


def compute_background(index: int):
    try:
        loaded = extract_core(index)
        return loaded.index, loaded.err_n[0]
    except FileNotFoundError:
        pass

    if verbose:
        print(f"-- computing background {index} @ {rank}")

    ics, params = sampler.get_sample(index)

    N_evo = np.linspace(0, 3000, step_density * 3000)
    background = model.backEvolve(
        N_evo,
        ics,
        params,
        tols,
        True,
        tmag_bg,
    )  # Raw background

    if isinstance(background, int):

        if background == -32:
            return SampleCore.int_background(index).get_last_status()

        elif background == -33:
            return SampleCore.timeout_background(index).get_last_status()

        else:
            raise ValueError(background)

    # Update evolution to correspond to epsilon > 1 end
    N_evo = background.T[0]
    Nend = N_evo[-1]

    # check for exit condition in field space constraints
    if len(fspace_terminate) > 0:

        _idx = np.inf

        for idx, step in enumerate(background):

            for fidx, bounds in fspace_terminate.items():

                if bounds[0] <= step[1 + fidx] <= bounds[1]:
                    _idx = np.min([_idx, idx])

        if np.isinf(_idx):
            return SampleCore.field_space_violation(index).get_last_status()

        background = background[:_idx + 1]
        N_evo = background.T[0]
        Nend = N_evo[-1]

    if len(fspace_reject) > 0:

        for idx, step in enumerate(background):

            for fidx, bounds in fspace_reject.items():

                if bounds[0] <= step[1 + fidx] <= bounds[1]:
                    SampleCore.field_space_violation(index).get_last_status()

    if Nend < Nmin:
        return SampleCore.inflation_too_short(index).get_last_status()

    if isinstance(Nend, tuple):
        return SampleCore.int_background(index).get_last_status()

    if isinstance(background, tuple):
        return SampleCore.int_background(index).get_last_status()

    back_adj = py_scripts.adjust_back(background)  # Adjusted background

    nexit_adj = py_scripts.compute_Nexit_for_matching(model, back_adj, params)
    nexit = py_scripts.compute_Nexit_for_matching(model, background, params)

    if np.isnan(nexit_adj):
        return SampleCore.nexit_background(index).get_last_status()

    return SampleCore(index, 0, back_adj, background, nexit_adj, nexit, params).get_last_status()


def any_nan_inf(*vals):
    arr = np.array([*vals])

    return np.any(np.isnan(arr)) or np.any(np.isinf(arr))


class _HintTasks:
    index: int
    task_dict: dict


def compute_epsilon(index: int):
    sample_core = extract_core(index)

    if "eps" not in sample_core.results:
        if verbose:
            print(f"-- computing epsilon {sample_core.index} @ {rank}")
        eps = py_scripts.get_epsilon_data(sample_core.back, sample_core.params, sample_core.Nexit, model)

        failed = any_nan_inf(eps)

        sample_core.task_eps(failed, eps)


def compute_eta(index: int):
    sample_core = extract_core(index)

    if "eta" not in sample_core.results:
        if verbose:
            print(f"-- computing eta {sample_core.index} @ {rank}")
        eta = py_scripts.get_eta_data(sample_core.back, sample_core.params, sample_core.Nexit, model)

        failed = any_nan_inf(eta)

        sample_core.task_eta(failed, value=eta)


def compute_mij(index: int):
    sample_core = extract_core(index)

    if "mij" not in sample_core.results:

        if verbose:
            print(f"-- computing masses {sample_core.index} @ {rank}")

        masses_exit, masses_end = py_scripts.get_mass_data(sample_core.back, sample_core.params, sample_core.Nexit,
                                                           model)
        back = sample_core.back

        NEND = back.T[0][-1]

        N_adi_check = np.array([N for N in back.T[0] if N > NEND - Nadi])

        # Check for adiabatic limit: Criteria here requires that m^2 >= 2 H^2 over the last efolds of inflation
        # for all but one tachyonic mode.
        # precise nember defined as N_adiabatic in setup

        adi = True
        for _ in N_adi_check:
            masses = py_scripts.get_mass_data(sample_core.back, sample_core.params, _, model)

            n_tachyon = (masses < 0).sum()

            if n_tachyon != 1:
                adi = False
                break

            rest_heavy = np.all(np.sort(masses)[1:] > 2)

            if not rest_heavy or any_nan_inf(rest_heavy):
                adi = False
                break

        failed = any_nan_inf(masses_exit, masses_end)

        sample_core.task_mij(failed, adi, value=[masses_exit, masses_end])


def compute_obs(data: _HintTasks):
    task_dict = data.task_dict

    if task_dict['task_2pt'] is True:
        compute_2pf(data)

    if task_dict['task_3pt_eq'] is True:
        compute_3pf(data, eq=True, fo=False, sq=False)

    if task_dict['task_3pt_fo'] is True:
        compute_3pf(data, eq=False, fo=True, sq=False)

    if task_dict['task_3pt_sq'] is True:
        compute_3pf(data, eq=False, fo=False, sq=True)

    if len(task_dict['alpha']) > 0:
        for alpha, beta in zip(task_dict['alpha'], ['beta']):
            compute_3pf(data, alpha=alpha, beta=beta, eq=False, fo=False, sq=False)

    return 0


def compute_2pf(data: _HintTasks):
    sample_core = extract_core(data.index)

    if "2pf" not in sample_core.results:

        if verbose:
            print(f"-- computing 2pf {sample_core.index} @ {rank}")

        result = py_scripts.compute_spectral_index(
            model,
            sample_core.back,
            sample_core.params,
            tols,
            Nsub,
            Nexit=sample_core.Nexit,
            tmax=tmax_2pf
        )

        if isinstance(result, int):
            status = result
            result = np.array([_BAD, _BAD, _BAD], dtype=np.float64)
        else:
            status = 0

        sample_core.task_2pf(status, value=result)


def compute_3pf(data: _HintTasks, **task_kwargs):
    sample_core = extract_core(data.index)

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

    if alpha_beta not in sample_core.results:

        if verbose:
            print(f"-- computing fnl {alpha_beta} {sample_core.index} @ {rank}")

        result = py_scripts.compute_fnl(
            model,
            sample_core.back,
            sample_core.params,
            tols,
            Nsub,
            Nexit=sample_core.Nexit,
            tmax=tmax_3pf,
            **task_kwargs
        )

        if isinstance(result, int):
            status = result
            result = np.array([_BAD], dtype=np.float64)
        else:
            status = 0

        sample_core.task_3pf(status, alpha_beta, value=result)
