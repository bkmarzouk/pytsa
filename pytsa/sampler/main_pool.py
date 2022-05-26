import os.path

import numpy as np
import dill

from pytsa.sampler import pyt_methods
from pytsa.sampler.setup_sampler import SamplerMethods, APrioriSampler, LatinSampler


class TrackMethods:

    def __init__(self, methods: SamplerMethods, sampler: APrioriSampler or LatinSampler, index, cache):
        self.methods = methods
        self.sampler = sampler
        self.index = index
        self.cache = cache


def _dill_load(path):
    with open(path, "rb") as f:
        return dill.load(f)


def build_back_pool(loc: str, n_samples: int):
    sampler_path = os.path.join(loc, "sampler.run")
    sampler: APrioriSampler or LatinSampler = _dill_load(sampler_path)

    methods_path = os.path.abspath(os.path.join(loc, "..", "sampler.methods"))
    methods: SamplerMethods = _dill_load(methods_path)

    pool_data = np.empty(n_samples, dtype=TrackMethods)

    samples_dir = os.path.join(loc, "samples_core")

    for idx in range(n_samples):
        pool_data[idx] = TrackMethods(methods, sampler, idx, samples_dir)

    return pool_data


def build_obs_pool(loc: str, status_dict: list):
    sampler_path = os.path.join(loc, "sampler.run")
    sampler: APrioriSampler or LatinSampler = _dill_load(sampler_path)

    methods_path = os.path.abspath(os.path.join(loc, "..", "sampler.methods"))
    methods: SamplerMethods = _dill_load(methods_path)

    indices = {sd[0] for sd in status_dict if sd[1] == 0}

    n_samples = len(indices)

    pool_data = np.empty(n_samples, dtype=TrackMethods)

    samples_dir = os.path.join(loc, "samples_core")

    for _, idx in enumerate(indices):
        pool_data[_] = TrackMethods(methods, sampler, idx, samples_dir)

    return pool_data


def main(pool, args_dict: dict):
    n_samples = args_dict['n_samples']

    # NOTE: list(...) call required for serial pools
    back_pool = build_back_pool(os.path.join(args_dict['cwd'], args_dict['name']), n_samples)

    back_status = list(pool.map(pyt_methods.compute_background, back_pool))

    obs_pool = build_obs_pool(os.path.join(args_dict['cwd'], args_dict['name']), back_status)

    list(pool.map(pyt_methods.compute_epsilon, obs_pool))
    list(pool.map(pyt_methods.compute_eta, obs_pool))
    list(pool.map(pyt_methods.compute_mij, obs_pool))
    list(pool.map(pyt_methods.compute_2pf, obs_pool))
    list(pool.map(pyt_methods.compute_3pf, obs_pool))

    return 0
    #
    # build_catalogue(setup, args)
    #
    # # Print configuration to command line
    # print("-- Starting sampling routine: {}".format(setup.PyT))
    #
    # pool_data = []
    #
    # # Construct initial numpy array of samples (required to support latin sampling)
    # catalogue_path = build_catalogue(setup, args, path_only=True)
    #
    # for idx in range(n_samples):
    #     pd = _Data()
    #     pd.PyT = setup.PyT.__name__
    #     pd.model_number = idx
    #     pd.path = catalogue_path
    #     pd.tols = setup.tols
    #     pd.N_min = setup.N_min
    #     pd.N_adi = setup.N_adiabitc
    #     pd.N_sub = setup.N_sub_evo
    #     pool_data.append(pd)
    #
    # # get results from background (returns data objs. with updated ics, params and back vals
    # r = pool.map(pyt_methods.compute_background, pool_data)
    #
    # # TODO: cache and reload pool data to following list
    #
    # # redefine pool data with successful trajectories
    # pool_data = [dat for dat in r if isinstance(dat.background, np.ndarray)]
