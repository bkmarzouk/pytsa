import os
import shutil
import dill
import pickle as pk
import numpy as np

from pytransport.cache_tools import hash_pars
from pytransport.sampler.methods import samplers

apriori = "apriori"
latin = "latin"
default_cache = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "samplers"))

assert os.path.exists(default_cache)


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

        self.PyT = pyt_model
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


class Setup(_SamplingParameters):

    def __init__(self, pyt_model, sampler_name, cache_loc=default_cache):
        """
        Setup routine for sampler

        :param pyt_model: pytransport module
        """
        super(Setup, self).__init__(pyt_model)

        self.cache_loc = os.path.join(cache_loc, sampler_name)
        self.sampler_path = os.path.join(self.cache_loc, "sampler")

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

        if not os.path.exists(cache_loc):
            os.makedirs(cache_loc)

        self.hash_check()

        self._dump_self()
        src = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..", "methods", "run.py")
        dest = os.path.join(cache_loc, "run.py")
        shutil.copy(src, dest)

    def _dump_self(self):
        """
        Dump binary file of sampler
        """

        path = os.path.join(self.sampler_path)

        with open(path, "wb") as f:
            print("-- sampler configuration data @ {}".format(path))
            dill.dump(self, f)

    def hash_check(self):
        """
        If construction is attempted for a sampler in an existing location,
        then we compare the keys and only throw error if there is a mismatch

        """

        hash_path = os.path.join(self.cache_loc, "pars.pk")

        current_data = self.get_hash_data()

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


def build_catalogue(s: Setup, parsed_args, path_only=False):  # TODO: Remove....
    pars = [s.N_min, s.N_adiabitc, s.N_sub_evo, s.tols]
    pars += [parsed_args.seed, parsed_args.grid_seed, parsed_args.n_samples]

    hash = hash_pars(*pars)

    dir_name = "latin_" if parsed_args.latin else "apriori_"

    dir_name += hash

    dir_path = s.path

    path = os.path.join(dir_path, dir_name)

    if not os.path.exists(path):
        os.makedirs(path)

    samples_path = os.path.join(path, "catalogue.npy")

    if path_only:
        return samples_path

    if not os.path.exists(samples_path):
        if parsed_args.apriori:
            x = samplers.APriori(parsed_args.seed)
        else:
            x = samplers.LatinHypercube(parsed_args.n_samples, seed=parsed_args.seed,
                                        cube_seed=parsed_args.grid_seed)

        print("-- Constructing catalogue for ICs & Params")

        for f in range(s.nF):
            x.add_param(s.fields[f])
        for f in range(s.nF):
            x.add_param(s.dot_fields[f])
        for p in range(s.nP):
            x.add_param(s.params[p])

        if parsed_args.apriori:
            samples = x.get_samples(parsed_args.n_samples).T
        else:
            samples = x.get_samples().T

        for row_idx, row in enumerate(samples):

            if "sr" in row:
                f = np.array([*row[:s.nF]], dtype=np.float64)
                p = np.array([*row[-s.nP:]], dtype=np.float64)
                V = s.PyT.V(f, p)
                dV = s.PyT.dV(f, p)

                sr_ic = -dV / np.sqrt(3 * V)

                for idx in range(s.nF):
                    if samples[row_idx][s.nF + idx] == "sr":
                        samples[row_idx][s.nF + idx] = sr_ic[idx]

        np.save(samples_path, np.asarray(samples, dtype=np.float64))

    print("-- Sample ICs & Params @ {}".format(samples_path))


if __name__ == "__main__":
    import pyt_dquad_euclidean as model
    from pytransport.sampler.methods.samplers import Constant
    import scipy.stats as stats

    sampler_setup = Setup(model, "dquad_test")

    sampler_setup.set_analysis_params(tols=[1e-7, 1e-7])
    sampler_setup.set_field(0, 1, method=stats.uniform(-20, 20))
    sampler_setup.set_dot_field(0, method="sr")
    sampler_setup.set_dot_field(1, method=stats.norm(0, 1e-7))
    sampler_setup.set_param(0, 1, method=stats.loguniform(1e-6, 1e-3))
    sampler_setup.build_sampler()
