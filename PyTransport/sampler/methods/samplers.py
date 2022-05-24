import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

from pytransport.sampler.configs.rng_states import RandomStates





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
