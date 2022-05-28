import os.path

import numpy as np
import dill

from pytsa.sampler import pyt_methods
from pytsa.sampler.setup_sampler import SamplerMethods, APrioriSampler, LatinSampler


class TrackMethods:

    def __init__(self, methods: SamplerMethods, sampler: APrioriSampler or LatinSampler, index, cache, task_dict=None):
        self.methods = methods
        self.sampler = sampler
        self.index = index
        self.cache = cache
        self.task_dict = task_dict


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


def build_obs_pool(loc: str, status_dict: list, task_dict: dict):
    sampler_path = os.path.join(loc, "sampler.run")
    sampler: APrioriSampler or LatinSampler = _dill_load(sampler_path)

    methods_path = os.path.abspath(os.path.join(loc, "..", "sampler.methods"))
    methods: SamplerMethods = _dill_load(methods_path)

    indices = {sd[0] for sd in status_dict if sd[1] == 0}

    n_samples = len(indices)

    pool_data = np.empty(n_samples, dtype=TrackMethods)

    samples_dir = os.path.join(loc, "samples_core")

    for _, idx in enumerate(indices):
        pool_data[_] = TrackMethods(methods, sampler, idx, samples_dir, task_dict=task_dict)

    return pool_data


def print_out(s: str):
    print(f"\n-- {s}\n")


def main(pool, args_dict: dict):
    n_samples = args_dict['n_samples']

    sampler_run_dir = os.path.join(args_dict['cwd'], args_dict['name'])

    back_pool = build_back_pool(sampler_run_dir, n_samples)

    print_out("Computing background trajectories")
    back_status = list(pool.map(pyt_methods.compute_background, back_pool))
    print_out("Background complete.")

    # Gather successful trajectories
    obs_pool = build_obs_pool(sampler_run_dir, back_status, args_dict)

    # Compute further background data for successful trajectories

    print_out("Computing epsilon data")
    list(pool.map(pyt_methods.compute_epsilon, obs_pool))
    print_out("Epsilon complete.")

    print_out("Computing eta data")
    list(pool.map(pyt_methods.compute_eta, obs_pool))
    print_out("Eta complete.")

    print_out("Computing mass data")
    list(pool.map(pyt_methods.compute_mij, obs_pool))
    print_out("Masses complete.")

    print_out("Computing observables")
    list(pool.map(pyt_methods.compute_obs, obs_pool))
    print_out("Observables complete.")
