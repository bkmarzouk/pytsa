import os.path

import numpy as np
import pandas as pd
import dill

from pytsa.sampler import pyt_methods
from pytsa.sampler.setup_sampler import SamplerMethods, APrioriSampler, LatinSampler
from pytsa.cache_tools import hash_alpha_beta

LINE_SIZE = None

_BAD = np.nan


class TrackMethods:

    def __init__(self, n_samples: int, methods: SamplerMethods, sampler: APrioriSampler or LatinSampler, index, cache,
                 task_dict=None):
        self.n_samples = n_samples
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
        pool_data[idx] = TrackMethods(n_samples, methods, sampler, idx, samples_dir)

    return pool_data


def build_obs_pool(loc: str, status_dict: list, task_dict: dict):
    sampler_path = os.path.join(loc, "sampler.run")
    sampler: APrioriSampler or LatinSampler = _dill_load(sampler_path)

    methods_path = os.path.abspath(os.path.join(loc, "..", "sampler.methods"))
    methods: SamplerMethods = _dill_load(methods_path)

    indices = {sd[0] for sd in status_dict if sd[1] == 0}

    n_total = task_dict['n_samples']

    n_samples = len(indices)

    pool_data = np.empty(n_samples, dtype=TrackMethods)

    samples_dir = os.path.join(loc, "samples_core")

    for _, idx in enumerate(indices):
        pool_data[_] = TrackMethods(n_total, methods, sampler, idx, samples_dir, task_dict=task_dict)

    return pool_data


def print_out(s: str):
    print(f"\n-- {s}\n")


def write_results(args_dict: dict, sample_pool: np.ndarray):
    nF = sample_pool[0].methods.nF
    nP = sample_pool[0].methods.nP
    latexs = []

    for ii in range(nF):
        latexs.append(sample_pool[0].methods.latex_f[ii])
    for ii in range(nF):
        latexs.append(sample_pool[0].methods.latex_df[ii])
    for ii in range(nP):
        latexs.append(sample_pool[0].methods.latex_p[ii])

    latexs.append(r"\epsilon_*")
    latexs.append(r"\epsilon")

    latexs.append(r"\eta_*")
    latexs.append(r"\eta")

    for ii in range(nF):
        latexs.append(r"M_{*," + str(ii) + "}^2/H_*^2")
    for ii in range(nF):
        latexs.append(r"M_{ii}^2/H^2".format(ii=ii))

    # headers for field ics, dot field ics and params
    headers = [f'f_{idx}' for idx in range(nF)]
    headers += [f'v_{idx}' for idx in range(nF)]
    headers += [f'p_{idx}' for idx in range(nP)]

    # headers for sr pars
    headers += ['eps_exit', 'eps_end', 'eta_exit', 'eta_end']

    # headers for mass eigenvalues
    headers += [f'm_exit_{n}' for n in range(nF)]
    headers += [f'm_end_{n}' for n in range(nF)]

    # get relevant keys for calling results in order
    obs_keys = []

    # headers for 2pf
    if args_dict['task_2pt']:
        obs_keys.append("2pf")

        headers.append("ns")
        headers.append("running")
        headers.append("As")

        latexs.append("n_s")
        latexs.append("d n_s / dk")
        latexs.append("A_s")

    # headers for template 3pf
    for ext in ['eq', 'fo', 'sq']:
        full = f"task_3pt_{ext}"
        if args_dict[full]:
            obs_keys.append(ext)
            headers.append(f"fnl_{ext}")

            latexs.append(r"F_\mathrm{NL}^{" + ext + "}")

    # headers for custom 3pf
    for alpha, beta in zip(args_dict['alpha'], args_dict['beta']):
        key = hash_alpha_beta(alpha, beta)
        obs_keys.append(key)
        headers.append(f"fnl_{key}")
        latexs.append(r"F_\mathrm{NL}(" + f"{alpha}, {beta}" + ")")

    # initialize as zeros array
    raw = np.zeros((args_dict['n_samples'], len(headers)), dtype=np.float64)

    for item in sample_pool:
        sample_path = os.path.join(item.cache, "sample.%06d" % item.index)

        with open(sample_path, "rb") as f:
            sample = dill.load(f)

        row_data = sample.get_row_data(*obs_keys)

        if isinstance(row_data, float) and np.isnan(row_data):
            raw[sample.index][:] = _BAD

        else:
            ics, pars = item.sampler.get_sample(item.index)
            raw[sample.index] = np.concatenate((ics, pars, row_data))

    df = pd.DataFrame(raw, columns=headers)

    sampler_run_dir = os.path.join(args_dict['cwd'], args_dict['name'])
    results_path = os.path.join(sampler_run_dir, "pandas_{}.df".format(args_dict['name']))

    if os.path.exists(results_path):
        os.remove(results_path)

    df.to_pickle(results_path)

    res_getdist = np.ones((raw.shape[0], raw.shape[1] + 2), dtype=np.float64)

    res_getdist[:, 2:] = raw

    Nrows = len(res_getdist)

    # Filter any rows that contain bad data values that would corrupt getdist analysis
    for _idx in range(len(res_getdist)):
        idx = Nrows - 1 - _idx
        r = res_getdist[idx]
        if np.any(np.isinf(r)) or np.any(np.isnan(r)) or np.any(r == _BAD):
            res_getdist = np.delete(res_getdist, idx, 0)

    name = args_dict['name']

    getdist_r_path = os.path.join(sampler_run_dir, f"getdist_{name}.txt")

    with open(getdist_r_path, "w") as f:
        for row_data in res_getdist:
            f.write(" ".join(map(str, row_data)) + "\n")

    getdist_p_path = os.path.join(sampler_run_dir, f"getdist_{name}.paramnames")
    with open(getdist_p_path, "w") as f:
        for h, l in zip(headers, latexs):
            f.write(f"{h} {l}\n")


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

    print("Writing results")
    write_results(args_dict, back_pool)
    print("All done.")
