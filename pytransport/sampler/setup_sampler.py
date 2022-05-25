import os
import shutil
import dill
import pickle as pk
import numpy as np
import scipy.stats
from pytransport.cache_tools import hash_alpha_beta
from pytransport.sampler.rng_states import RandomStates
from pytransport.sampler.mpi_helpers import make_dir, single_proc_task, barrier, rank

default_cache = os.path.abspath(os.path.join(os.path.dirname(__file__), "", "..", "..", "samplers"))

assert os.path.exists(default_cache)


class Constant(object):

    def __init__(self, val):
        """
        Mimics methods from scipy.stats for constructing samples via rvs or ppf. Returns constant in either case
        for any value

        :param val: constant value to return
        """
        self.val = val
        self.random_state = None

    def rvs(self, n=1):
        return np.array([self.val]).repeat(n) if n > 1 else self.val

    def ppf(self, *x):
        return self.val

    @classmethod
    def slow_roll(cls):
        return cls("sr")


def _dist2pars(d):
    """
    Translate scipy distributions into hashable params

    :param d: scipy.stats._xxx_distns.yyy_gen or Constant instance
    :return: [dist_name, loc, scale] or [const]
    """
    if isinstance(d, Constant):
        return [d.rvs()]

    assert hasattr(d, "dist"), d
    assert hasattr(d, "args"), d

    pars = [str(d.dist).split("_gen")[0].split(".")[-1]]
    for a in d.args:
        pars.append(a)

    return pars


def _check_method(m):
    """
    Checks methods are available for dists called by sampling mode

    :param m: method
    """
    assert hasattr(m, "rvs"), m
    assert hasattr(m, "ppf"), m


class _SamplingParameters(object):

    def __init__(self, pyt_model):
        """
        Container for sampling parameters

        :param pyt_model: pytransport module
        """

        self.model = pyt_model
        self.nF = pyt_model.nF()
        self.nP = pyt_model.nP()

        self.fields = {f: None for f in range(self.nF)}
        self.dot_fields = {f: None for f in range(self.nF)}
        self.params = {p: None for p in range(self.nP)}

        self.latex_f = {f: None for f in range(self.nF)}
        self.latex_df = {f: None for f in range(self.nF)}
        self.latex_p = {p: None for p in range(self.nP)}

    def _check(self):
        for f in range(self.nF):
            assert self.fields[f] is not None, f
            assert self.dot_fields[f] is not None, f
        for p in range(self.nP):
            assert self.params[p] is not None, p

    def set_field(self, *indices, method, latex: str or None = None, record=True):
        """
        Set initial values for fields

        :param indices: field index
        :param method: distribution to sample
        :param latex: latex definition
        :param record: if True, records IC in results
        """
        if method == "sr" or method == "SR":
            method = Constant.slow_roll()

        _check_method(method)

        for idx in indices:
            self.fields[idx] = method
            self.latex_f[idx] = (latex, record)

    def set_dot_field(self, *indices, method, latex: str or None = None, record=True):
        """
        Set initial values for initial field derivatives w.r.t. cosmic time

        :param indices: field index
        :param method: distribution to sample
        :param latex: latex definition
        :param record: if True, records IC in results
        """
        if method == "sr" or method == "SR":
            method = Constant.slow_roll()

        _check_method(method)

        for idx in indices:
            self.dot_fields[idx] = method
            self.latex_df[idx] = (latex, record)

    def set_param(self, *indices, method, latex: str or None = None, record=True):
        """
        Set model parameter values

        :param indices: parameter index
        :param method: distribution to sample
        :param latex: latex definition
        :param record: if True, records IC in results
        """
        _check_method(method)

        for idx in indices:
            self.params[idx] = method
            self.latex_p[idx] = (latex, record)

    def _get_hash_data(self):

        hd = {}

        for idx in range(self.nF):
            hd[f"f_{idx}"] = _dist2pars(self.fields[idx])
            hd[f"v_{idx}"] = _dist2pars(self.dot_fields[idx])

        for idx in range(self.nF):
            hd[f"p_{idx}"] = _dist2pars(self.params[idx])

        return hd


class SamplerMethods(_SamplingParameters):

    def __init__(self, pyt_model, sampler_name, cache_loc=default_cache):
        """
        Setup routine for sampler

        :param pyt_model: pytransport module
        """
        super(SamplerMethods, self).__init__(pyt_model)

        self.cache_loc = os.path.join(cache_loc, sampler_name)
        self.sampler_path = os.path.join(self.cache_loc, "sampler.methods")

        self.fieldspace_reject = {}
        self.dotfieldspace_reject = {}

        self.fieldspace_end_inflation = {}
        self.dotfieldspace_end_inflation = {}

        self._analysis_pars_set = False

    def add_fieldspace_constraint(self, fidx, fmin=-np.inf, fmax=np.inf, dot_field=False):
        """
        Define region of fieldspace where inflation is permitted. If deviation occurs, the sample is rejected.

        :param fidx: Field index
        :param fmin: minimum allowed field value. By default, -inf, i.e. no lower bound
        :param fmax: maximum allowed field value. By default, +inf, i.e. no upper bound
        :param dot_field: if True, constraint applies to field velocities. Else field value
        """

        if dot_field:
            if fidx in self.dotfieldspace_reject:
                raise AttributeError("Attribute already defined for field number {}".format(fidx))
            self.dotfieldspace_reject[fidx] = np.array([fmin, fmax])

        else:
            if fidx in self.fieldspace_reject:
                raise AttributeError("Attribute already defined for field number {}".format(fidx))
            self.fieldspace_reject[fidx] = np.array([fmin, fmax])

    def add_end_inflation_constraint(self, fidx, fmin, fmax, dot_field=False):
        """
        Define region of fieldspace where inflation instantaneously ends. E.g. in d-brane inflation.
        The first efold N when f \in [fmin, fmax] will terminate inflation. Sample(s) will be accepted
        if this condition is met, in conjunction with a sufficient number of efoldings.

        :param fidx: Field index
        :param fmin: minimum allowed field value. By default,
        :param fmax: maximum allowed field value. By default, +inf, i.e. no upper bound
        :param dot_field: if True, constraint applies to field velocities. Else field value
        """

        if dot_field:
            if fidx in self.dotfieldspace_end_inflation:
                raise AttributeError("Attribute already defined for field number {}".format(fidx))
            self.dotfieldspace_end_inflation[fidx] = np.array([fmin, fmax])

        else:
            if fidx in self.fieldspace_end_inflation:
                raise AttributeError("Attribute already defined for field number {}".format(fidx))
            self.fieldspace_end_inflation[fidx] = np.array([fmin, fmax])

    def add_mutual_constraint(self, fields=None, dotfields=None, **kwargs):
        assert 0, "currently unsupported"

    def set_analysis_params(self, N_sub_evo=6., tols=None, N_adiabatic=1., N_min=60.):
        """
        Set integration parameters and efoldings for ICs, adiabicity and minimal inflation

        :param N_sub_evo: efoldings of sub-horizon evolution to track
        :param tols: integration tols (abs, rel)
        :param N_adiabatic: Number of efoldings to probe the adiabaitic limit
        :param N_min: min efolds of inflation
        """
        if tols is None:
            tols = [1e-8, 1e-8]
        self.N_sub_evo = N_sub_evo
        self.tols = tols
        self.N_adiabitc = N_adiabatic
        self.N_min = N_min
        self._analysis_pars_set = True

    def get_hash_data(self) -> dict:

        assert self._analysis_pars_set, "Analysis parameters not set!"

        base = self._get_hash_data()

        prefixes = ['fields_rej', 'dotfields_rej', 'fields_end', 'dotfields_end']

        for conditions, p in zip([self.fieldspace_reject, self.dotfieldspace_reject,
                                  self.fieldspace_end_inflation, self.dotfieldspace_end_inflation], prefixes):
            for key in conditions:
                base[p + f"_{key}"] = conditions[key]

        keys = ['Nsub', 'tols', 'Nadi', 'Nmin']

        for value, key in zip([self.N_sub_evo, self.tols, self.N_adiabitc, self.N_min], keys):
            base[key] = value

        return base

    def build_sampler(self):
        """
        Execute construction of sampling module
        """

        self._check()

        if not self._analysis_pars_set:
            self.set_analysis_params()  # use defaults

        cache_loc = self.cache_loc

        make_dir(cache_loc)

        self.hash_check()

        self._dump_self()

        dest = os.path.join(cache_loc, "run.py")
        if rank == 0 and not os.path.exists(dest):
            src = os.path.join(os.path.abspath(os.path.dirname(__file__)), "run.py")
            shutil.copy(src, dest)

        barrier()

    def _dump_self(self):
        """
        Dump binary file of sampler
        """

        path = os.path.join(self.sampler_path)

        if not os.path.exists(path) and rank == 0:
            with open(path, "wb") as f:
                print("-- sampler configuration data @ {}".format(path))
                dill.dump(self, f)

        barrier()

    def hash_check(self):
        """
        If construction is attempted for a sampler in an existing location,
        then we compare the keys and only throw error if there is a mismatch

        """

        hash_path = os.path.join(self.cache_loc, "sampler.pars")

        current_data = self.get_hash_data()

        if rank == 0:

            if not os.path.exists(hash_path):
                with open(hash_path, "wb") as f:
                    pk.dump(current_data, f)
            else:
                with open(hash_path, "rb") as f:
                    existing_data: dict = pk.load(f)

                # check for key data
                assert current_data.keys() == existing_data.keys(), [current_data.keys(), existing_data.keys()]

                # check values match
                for k in current_data.keys():
                    cur_val = current_data[k]  # Current value
                    exi_val = existing_data[k]  # Existing value
                    assert current_data[k] == existing_data[k], f"KEY ERROR @ {k}:\n{cur_val} != {exi_val}"

        barrier()

    def sr2val(self, idx: int or np.ndarray, fvals: np.ndarray, pvals: np.ndarray):
        """
        Get slow roll value for field index/indices
        :param idx:
        :param fvals:
        :param pvals:
        :return:
        """

        V = self.model.V(fvals.astype(float), pvals)
        dV = self.model.dV(fvals.astype(float), pvals)

        return -dV[idx] / np.sqrt(3 * V)


class _SamplerRoutine(object):

    def __init__(self, random_states: RandomStates, field_methods: list, dotfield_methods: list, param_methods: list):
        """
        Base class for samplers.

        :param seed: if given, seed is used for prng before 'get_samples' call in derived function.
                        Note: in LatinHypercube method, this is called earlier, during the distribution
                                of grid squares.
        """
        self.random_states = random_states

        self.n_field = len(field_methods)
        self.field_methods = field_methods

        assert len(dotfield_methods) == self.n_field, [len(dotfield_methods), self.n_field]
        self.dotfield_methods = dotfield_methods

        self.n_params = len(param_methods)
        self.param_methods = param_methods

        if isinstance(self, LatinHypercube):
            print("-- I am latin")
            for m in self.field_methods + self.dotfield_methods + self.param_methods:
                assert hasattr(m, "ppf"), m

        elif isinstance(self, APriori):
            print("-- I am apriori")
            for m in self.field_methods + self.dotfield_methods + self.param_methods:
                assert hasattr(m, "rvs"), m

        else:
            raise TypeError(self)

    def get_sample(self, *args, **kwargs) -> tuple:
        pass


class APriori(_SamplerRoutine):

    def __init__(self, random_states: RandomStates, field_methods: list, dotfield_methods: list,
                 param_methods: list):
        """
        Draws values apriori from statistical distributions defined with paramters

        :param seed: integer or None, defining a random (or absence of) random seed
        """
        super(APriori, self).__init__(random_states, field_methods, dotfield_methods, param_methods)

    def get_sample(self, sample_index: int):
        """
        Get samples for parameters as defined by corresponding statistical distributions

        :param n_samples: number of samples
        :return: array of samples with shape (n_samples, n_params)
        """

        dotfields_out = np.zeros(self.n_field * 2, dtype=object)
        params_out = np.zeros(self.n_params, dtype=float)

        state = self.random_states.get_state(sample_index)

        for idx in range(self.n_field):
            m = self.field_methods[idx]
            m.random_state = state
            dotfields_out[idx] = m.rvs()

            m = self.dotfield_methods[idx]
            m.random_state = state
            dotfields_out[idx + self.n_field] = m.rvs()

        for idx in range(self.n_params):
            m = self.param_methods[idx]
            m.random_state = state
            params_out[idx] = m.rvs()

        return dotfields_out, params_out


class LatinHypercube(_SamplerRoutine):

    def __init__(self, random_states: RandomStates, field_methods: list, dotfield_methods: list,
                 param_methods: list):
        """
        Constructs latin hypercube for inverse sampling the CDF of a given distribution

        :param n_cells: number of cells (along one dimension) to define sample regions
        :param seed: prng seed for constructing samples
        :param cube_seed: prng seed for constructing hypercube
        """
        super(LatinHypercube, self).__init__(random_states, field_methods, dotfield_methods, param_methods)
        self.n_samples = random_states.n_states - 1
        self.grid = self.build_grid()
        g = np.linspace(0, 1, self.n_samples + 1)
        self.unif = list(zip(g, np.roll(g, -1)))[:-1]  # const density intervals

    def coords_1d(self, grid_state):
        """
        Construct non-repeating list of integers on the interval [0, n_cells-1], that represent
        sampling intervals for a single parameter of a unique sample

        :return: interval coords
        """

        n_samples = self.n_samples

        avail = list(range(n_samples))
        coords = np.zeros(n_samples, dtype=int)

        for ii in range(n_samples):
            m = scipy.stats.randint(0, len(avail))
            m.random_state = grid_state
            coords[ii] = avail.pop(m.rvs())

        return coords, grid_state

    def build_grid(self):
        """
        Construct grid of coordinates representing sampling intervals for the n-dimensional parameter space

        :return: nxn coordinate grid
        """

        n_tot = self.n_params + self.n_field * 2

        coords = np.zeros((n_tot, self.n_samples), dtype=int)

        grid_state = self.random_states.get_state(0)

        for ii in range(n_tot):
            coords[ii], grid_state = self.coords_1d(grid_state)

        return coords.T

    def get_sample(self, sample_index: int):
        """
        Get samples for parameters as defined by corresponding statistical distributions

        :param verbose: print statements if True
        :return: array of samples with shape (n_cells, n_params)
        """

        # +1 to state index since we reserve the zeroth for grid building
        state = self.random_states.get_state(sample_index + 1)

        fdf_out = np.zeros(2 * self.n_field, dtype=object)
        par_out = np.zeros(self.n_params, dtype=float)

        grid_indices = self.grid[sample_index]
        unif_regions = self.unif

        for idx in range(self.n_field):
            lb, ub = unif_regions[grid_indices[idx]]  # gets boundaries for subnterval between 0 and 1
            m = scipy.stats.uniform(loc=lb, scale=ub - lb)  # initialize rng for region
            m.random_state = state  # update random state to match for current sample
            unif_rv = m.rvs()  # random point from cdf on [0, 1]
            fdf_out[idx] = self.field_methods[idx].ppf(unif_rv)  # Get sampled value

            lb, ub = unif_regions[grid_indices[idx + self.n_field]]  # gets boundaries for subnterval between 0 and 1
            m = scipy.stats.uniform(loc=lb, scale=ub - lb)  # initialize rng for region
            m.random_state = state  # update random state to match for current sample
            unif_rv = m.rvs()  # random point from cdf on [0, 1]
            fdf_out[idx + self.n_field] = self.dotfield_methods[idx].ppf(unif_rv)  # Get sampled value

        for idx in range(self.n_params):
            lb, ub = self.unif[grid_indices[idx + self.n_field * 2]]
            m = scipy.stats.uniform(loc=lb, scale=ub - lb)  # initialize rng for region
            m.random_state = state  # update random state to match for current sample
            unif_rv = m.rvs()  # random point from cdf on [0, 1]
            par_out[idx] = self.param_methods[idx].ppf(unif_rv)  # Get sampled value

        return fdf_out, par_out


class _XSampler:

    def __init__(self, is_apriori: bool, random_states: RandomStates, hyp_pars: SamplerMethods):
        self.random_states = random_states
        self.hyp_pars = hyp_pars

        fields = [hyp_pars.fields[idx] for idx in range(hyp_pars.nF)]
        dotfields = [hyp_pars.dot_fields[idx] for idx in range(hyp_pars.nF)]
        params = [hyp_pars.params[idx] for idx in range(hyp_pars.nP)]

        if is_apriori:
            self.sampler = APriori(random_states, fields, dotfields, params)
        else:
            self.sampler = LatinHypercube(random_states, fields, dotfields, params)

    def get_sample(self, idx: int):
        _fdf, pars = self.sampler.get_sample(idx)

        if np.any(_fdf == "sr"):
            sr_idx = np.where(_fdf == "sr")[0]
            _fdf[sr_idx] = self.hyp_pars.sr2val(sr_idx - self.hyp_pars.nF, _fdf[:self.hyp_pars.nF], pars)

        # TODO: This is an inefficient way of type-casting the fdf results from object to float
        #       NOTE: WE have object type to accomodate the 'sr' key; maybe define this as a floating
        #             point constant such that we don't have to do this?
        fdf = np.zeros(self.hyp_pars.nF * 2, dtype=float)
        for idx, val in enumerate(_fdf):
            fdf[idx] = val

        return fdf, pars


class APrioriSampler(_XSampler):

    def __init__(self, random_states: RandomStates, hyp_pars: SamplerMethods):
        super(APrioriSampler, self).__init__(True, random_states, hyp_pars)


class LatinSampler(_XSampler):

    def __init__(self, random_states: RandomStates, hyp_pars: SamplerMethods):
        super(LatinSampler, self).__init__(False, random_states, hyp_pars)


def build_results_template(cache, n_samples, n_vals):
    path = os.path.join(cache, "data.npy")
    if not os.path.exists(path) and rank == 0:
        np.save(path, np.zeros((n_samples, n_vals), dtype=np.float64))

    barrier()


def job_config(task_pars: dict):
    name = task_pars['name']
    n_samples = task_pars['n_samples']
    entropy = task_pars['entropy']
    apriori = task_pars['apriori']
    sampler_cwd = task_pars['cwd']

    n_states = n_samples if apriori else n_samples + 1

    root_dir = os.path.abspath(os.path.join(sampler_cwd, name))

    make_dir(root_dir)
    barrier()

    tasks_path = os.path.join(root_dir, "sampler.tasks")

    if rank == 0:

        if os.path.exists(tasks_path):

            with open(tasks_path, "rb") as f:

                cached_task_pars: dict = pk.load(f)

            assert cached_task_pars.keys() == task_pars.keys(), [cached_task_pars.keys(), task_pars.keys()]

            for k in cached_task_pars.keys():
                parsed_par = task_pars[k]
                cached_par = cached_task_pars[k]

                assert parsed_par == cached_par or k == "n_procs", f"Parameter error @ {k}: {parsed_par} != {cached_par}"

        else:

            with open(tasks_path, "wb") as f:

                pk.dump(task_pars, f)

        random_states = RandomStates.from_cache(root_dir, entropy=entropy, n_states=n_states)

        sampler_path = os.path.join(root_dir, "sampler.run")

        if not os.path.exists(sampler_path):
            methods_path = os.path.abspath(os.path.join(sampler_cwd, "sampler.methods"))

            with open(methods_path, "rb") as f:
                sampler_methods = dill.load(f)  # DILL??

            _m = APrioriSampler if apriori else LatinSampler

            sampler_run = _m(random_states, sampler_methods)

            with open(sampler_path, "wb") as f:
                dill.dump(sampler_run, f)

        samples_core_dir = os.path.join(root_dir, "samples_core")
        make_dir(samples_core_dir)

        result_dirs = [(x, 2) for x in ['mij', 'epsilon', 'eta']]

        if 'task_2pt' in task_pars:
            result_dirs.append(('ns_alpha', 2))

        for ext in ['eq', 'fo', 'sq']:
            task_name = f'task_3pt_{ext}'
            if task_name in task_pars:
                result_dirs.append((f'fnl_{ext}', 1))

        for a, b in zip(task_pars['alpha'], task_pars['beta']):
            result_dirs.append((f'fnl_{hash_alpha_beta(a, b)}', 1))

        for rd in result_dirs:
            cache = os.path.join(root_dir, rd[0])

            if not os.path.exists(cache):
                os.makedirs(cache)

    barrier()


if __name__ == "__main__":

    import pyt_dquad_euclidean as model
    import scipy.stats as stats

    sampler_setup = SamplerMethods(model, "dquad_test")

    sampler_setup.set_analysis_params(tols=[1e-5, 1e-5])
    sampler_setup.set_field(0, 1, method=stats.uniform(-20, 40))
    sampler_setup.set_dot_field(0, method="sr")
    sampler_setup.set_dot_field(1, method="sr")
    sampler_setup.set_param(0, 1, method=stats.loguniform(1e-6, 1e-3))
    sampler_setup.build_sampler()

    make_demo_fig = False

    if make_demo_fig:

        import matplotlib.pyplot as plt

        N = 2

        n_samples = 1000

        states_lh = RandomStates(n_samples + 1, entropy=438473848392)
        states_ap = RandomStates(n_samples, entropy=438473848392)

        m1 = scipy.stats.uniform(loc=-20, scale=40)
        m2 = scipy.stats.norm(1e-3)
        m3 = scipy.stats.lognorm(1e-6, 1e-3)

        lh = LatinHypercube(states_lh, [m1] * N, [m2] * N, [m3] * N)
        ap = APriori(states_ap, [m1] * N, [m2] * N, [m3] * N)

        lh_fdfs = np.zeros((n_samples, N * 2), dtype=float)
        lh_pars = np.zeros((n_samples, N), dtype=float)

        ap_fdfs = lh_fdfs.copy()
        ap_pars = lh_pars.copy()

        for ii in range(n_samples):
            lh_fdfs[ii], lh_pars[ii] = lh.get_sample(ii)
            ap_fdfs[ii], ap_pars[ii] = ap.get_sample(ii)

        fig, axs = plt.subplots(2, 3, figsize=(10, 5))

        for idx, ax in enumerate(axs.flatten()[:4]):
            ax.hist(lh_fdfs.T[idx], density=True, bins="auto", alpha=0.3)
            ax.hist(ap_fdfs.T[idx], density=True, bins="auto", alpha=0.3)

        for idx, ax in enumerate(axs.flatten()[4:]):
            ax.hist(lh_pars.T[idx], density=True, bins="auto", alpha=0.3)
            ax.hist(ap_pars.T[idx], density=True, bins="auto", alpha=0.3)

        plt.show()
        plt.tight_layout()
