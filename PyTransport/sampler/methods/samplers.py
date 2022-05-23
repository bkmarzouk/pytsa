import numpy as np
import scipy.stats
from pytransport.sampler.methods.rng import RNGAPriori, RNGLatin


class Constant(object):

    def __init__(self, val):
        """
        Mimics methods from scipy.stats for constructing samples via rvs or ppf. Returns constant in either case
        for any value

        :param val: constant value to return
        """
        self.val = val

    def rvs(self, n=1):
        return np.array([self.val]).repeat(n) if n > 1 else self.val

    def ppf(self, *x):
        return self.val

    @classmethod
    def slow_roll(cls):
        return cls("sr")


class _Sampler(object):

    def __init__(self, seed=None):
        """
        Base class for samplers.

        :param seed: if given, seed is used for prng before 'get_samples' call in derived function.
                        Note: in LatinHypercube method, this is called earlier, during the distribution
                                of grid squares.
        """
        self.seed = seed
        self.n_params = 0
        self.dists = []

        self.n_consts = 0
        self.consts = []

        self.order = {}

    def add_param(self, dist=scipy.stats.uniform(0, 1)):
        """
        Adds parameter to be sampled via a specified distribution.

        :param dist: distribution object. Must have methods 'rvs' for Apriori sampling, or 'ppf' for hypercube.
                     These conditions are satisfied by calling the standard scipy.stats modules.
        """

        dist = Constant.slow_roll() if dist == "sr" else dist

        if issubclass(_Sampler, LatinHypercube):
            assert hasattr(dist, "ppf"), "Must have inverse cdf, defined as 'ppf' method."
        if issubclass(_Sampler, APriori):
            assert hasattr(dist, "rvs"), "Must have random variables generator, as defined as 'rvs' method."

        if isinstance(dist, Constant):
            self.consts.append(dist)
            self.n_consts += 1
            self.order[len(self.order)] = "c"
        else:
            self.dists.append(dist)
            self.n_params += 1
            self.order[len(self.order)] = "p"

    def get_samples(self):
        assert 0, "override method"


class APriori(_Sampler):

    def __init__(self, seed=None):
        """
        Draws values apriori from statistical distributions defined with paramters

        :param seed: integer or None, defining a random (or absence of) random seed
        """
        super(APriori, self).__init__(seed)

    def get_samples(self, n_samples=1):
        """
        Get samples for parameters as defined by corresponding statistical distributions

        :param n_samples: number of samples
        :return: array of samples with shape (n_samples, n_params)
        """
        out = np.zeros((self.n_params + self.n_consts, n_samples), dtype=object)

        p_samples = np.zeros((self.n_params, n_samples))
        c_samples = np.array([c.rvs() for c in self.consts]).repeat(n_samples).reshape((self.n_consts, n_samples))

        np.random.seed(self.seed)

        for idx, d in enumerate(self.dists):
            p_samples[idx] = d.rvs(n_samples)

        p_count = 0
        c_count = 0
        for idx in range(len(self.order)):
            if self.order[idx] == "p":
                out[idx] = p_samples[p_count]
                p_count += 1
            else:
                out[idx] = c_samples[c_count]
                c_count += 1

        return out


class LatinHypercube(_Sampler):

    def __init__(self, n_cells=10, seed=None, cube_seed=None):
        """
        Constructs latin hypercube for inverse sampling the CDF of a given distribution

        :param n_cells: number of cells (along one dimension) to define sample regions
        :param seed: prng seed for constructing samples
        :param cube_seed: prng seed for constructing hypercube
        """
        super(LatinHypercube, self).__init__(seed)
        self.n_cells = n_cells
        self.cube_seed = cube_seed

    def coords_1d(self):
        """
        Construct non-repeating list of integers on the interval [0, n_cells-1], that represent
        sampling intervals for a single parameter of a unique sample

        :return: interval coords
        """

        n_cells = self.n_cells

        avail = list(range(n_cells))
        coords = np.zeros(n_cells, dtype=int)

        for ii in range(n_cells):
            coords[ii] = avail.pop(np.random.randint(n_cells - ii))

        return coords

    def build_grid(self):
        """
        Construct grid of coordinates representing sampling intervals for the n-dimensional parameter space

        :return: nxn coordinate grid
        """
        np.random.seed(self.cube_seed)

        coords = np.zeros((self.n_params, self.n_cells), dtype=int)

        for ii in range(self.n_params):
            coords[ii] = self.coords_1d()

        return coords

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

    def uniform_domains(self):
        """
        Splits the interval [0, 1] into equal sized subsets, s.t. Union([s1, s2], ..., [sn-1, sn]) = [0, 1]

        :return: uniform domains to sample from
        """
        g = np.linspace(0, 1, self.n_cells + 1)
        return list(zip(g, np.roll(g, -1)))[:-1]

    def get_samples(self, verbose=False):
        """
        Get samples for parameters as defined by corresponding statistical distributions

        :param verbose: print statements if True
        :return: array of samples with shape (n_cells, n_params)
        """

        if verbose:
            print("-- building grid")
        coords = self.build_grid()

        if verbose:
            print("-- building domains")
        domains = self.uniform_domains()

        unifs = np.zeros((self.n_cells, self.n_params), dtype=float)

        if verbose:
            print("-- building uniform samples")
        for idx, domain_indices in enumerate(coords.T):
            unifs[idx] = [np.random.uniform(*domains[d]) for d in domain_indices]

        if verbose:
            print("-- building samples")

        p_samples = np.zeros((self.n_cells, self.n_params), dtype=object)
        for ii in range(self.n_cells):
            p_samples[ii] = self.v2s(*unifs[ii])

        if verbose:
            print("-- done")

        p_samples = p_samples.T

        c_samples = np.zeros((self.n_consts, self.n_cells), dtype=object)

        for idx, c in enumerate(self.consts):
            val = c.rvs()

            if val == "sr":
                pass

            c_samples[idx] = val

        out = np.zeros((self.n_params + self.n_consts, self.n_cells), dtype=object)

        p_count = 0
        c_count = 0
        for idx in range(len(self.order)):
            if self.order[idx] == "p":
                out[idx] = p_samples[p_count]
                p_count += 1
            else:
                out[idx] = c_samples[c_count]
                c_count += 1

        return out


if __name__ == "__main__":

    make_demo_fig = True

    if make_demo_fig:

        import os
        import matplotlib.pyplot as plt

        n_samples = 5000

        a = APriori()
        l = LatinHypercube(n_cells=n_samples)

        dists = [scipy.stats.uniform(0, 1), scipy.stats.norm(0, 1), scipy.stats.betaprime(5, 6)]

        for d in dists:
            a.add_param(d)
            l.add_param(d)

        samples_a = a.get_samples(n_samples)
        samples_l = l.get_samples()

        fig, axs = plt.subplots(1, 3, figsize=(12, 4))

        dist_lims = [[0, 1], [-5, 5], [0, 5]]

        for ax, sa, sl, d, l, t in zip(axs, samples_a, samples_l, dists, dist_lims,
                                       ['Uniform', 'Normal', '$\\beta$-prime']):

            ax.hist(sa, density=True, bins=40, alpha=0.33, label="Apriori")
            ax.hist(sl, density=True, bins=40, alpha=0.33, label="Latin")
            v = np.linspace(l[0], l[1], 1000)
            ax.plot(v, d.pdf(v), ls="--", c="k", lw=2, label="PDF")
            ax.set_title(t)
            ax.set_xlabel("$x$", size=13)
            ax.set_ylabel("$\mathbb{P}(x)$", size=13)

        axs[-1].legend()

        plt.tight_layout()

        plt.savefig(os.path.join("/home/kareem/cosmo-share", "apriori_latin_demo.pdf"), bbox_inches="tight")

        plt.show()
