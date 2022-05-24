import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

from pytransport.sampler.configs.rng_states import RandomStates


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


class _Sampler(object):

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


class APriori(_Sampler):

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
        params_out = np.zeros(self.n_params, dtype=object)

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


class LatinHypercube(_Sampler):

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

    def v2s(self, *values):
        """
        Translates uniformly sampled values on [0, 1] to inverse CDF values for each parameter

        :param values: uniformly distributed values from sampled intervals
        :return: samples of each parameter
        """

        out = np.zeros(self.n_params, dtype=object)

        for ii in range(self.n_params):
            out[ii] = self.dists[ii].ppf(values[ii])
        return out

    def get_sample(self, sample_index: int):
        """
        Get samples for parameters as defined by corresponding statistical distributions

        :param verbose: print statements if True
        :return: array of samples with shape (n_cells, n_params)
        """

        # +1 to state index since we reserve the zeroth for grid building
        state = self.random_states.get_state(sample_index + 1)

        fdf_out = np.zeros(2 * self.n_field, dtype=object)
        par_out = np.zeros(self.n_params, dtype=object)

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


if __name__ == "__main__":

    make_demo_fig = True

    if make_demo_fig:

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
