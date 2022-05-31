import os.path

import numpy as np
import pandas as pd
import dill

from . import pyt_methods
from .setup_sampler import SamplerMethods, APrioriSampler, LatinSampler
from ..cache_tools import hash_alpha_beta

_BAD = np.nan

reps = [
    ['short', "Nend < Nmin"],
    ['fspace', "Violated field space condition"],
    ['int_back', "Background integration error"],
    ['timeout_back', "Background integration timeout"],
    ['nexit', "Failed to compute acceptable horizon exit time"],
    ['eternal', "Could not find end of inflation"],
    ['eps', "Failed to compute epsilon slow-roll data"],
    ['eta', "Failed to compute eta slow-roll data"],
    ['mij', "Failed to compute Mij eigenvalue data"],
    ['ics_2pf', "Unable to find initial conditions for 2pf"],
    ['kexit_2pf', "Unable to find horizon exit mode for 2pf"],
    ['int_2pf', "2pf integration error"],
    ['timeout_2pf', "2pf integration timeout"],
    ["ics_3pf_", "Unable to find initial conditions for 3pf, "],
    ["kexit_3pf_", "Unable to find horizon exit mode for 3pf, "],
    ["int_3pf_", "3pf integration error, "],
    ["timeout_3pf_", "3pf integration timeout, "]
]


class SampleData:

    def __init__(self, n_samples: int, index: int, task_dict: dict or None = None):
        self.n_samples = n_samples
        self.index = index
        self.task_dict = task_dict


def _dill_load(path):
    with open(path, "rb") as f:
        return dill.load(f)


def build_obs_pool(indices: list, task_dict: dict):
    pool_data = np.empty(len(indices), dtype=SampleData)

    n_total = task_dict['n_samples']

    for _, idx in enumerate(indices):
        pool_data[_] = SampleData(n_total, idx, task_dict=task_dict)

    return pool_data


def print_out(s: str):
    print(f"\n-- {s}\n")


def write_results(args_dict: dict, indices: np.ndarray):

    path = os.path.join(args_dict['cwd'], args_dict['name'], 'sampler.run')
    samples_dir = os.path.join(args_dict['cwd'], args_dict['name'], 'samples_core')

    with open(path, "rb") as f:
        sampler: APrioriSampler or LatinSampler = dill.load(f)

    nF = sampler.hyp_pars.nF
    nP = sampler.hyp_pars.nP
    latexs = []

    for ii in range(nF):
        latexs.append(sampler.hyp_pars.latex_f[ii])
    for ii in range(nF):
        latexs.append(sampler.hyp_pars.latex_df[ii])
    for ii in range(nP):
        latexs.append(sampler.hyp_pars.latex_p[ii])

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

    errors_data = {}

    for idx in indices:
        sample_path = os.path.join(samples_dir, "sample.%06d" % idx)

        with open(sample_path, "rb") as f:
            sample = dill.load(f)

        row_data = sample.get_row_data(*obs_keys)

        for n, s in zip(sample.err_n, sample.err_s):

            if n != 0:

                if s in errors_data:
                    errors_data[s] += 1
                else:
                    errors_data[s] = 1

        if isinstance(row_data, float) and np.isnan(row_data):
            raw[sample.index][:] = _BAD

        else:
            ics, pars = sampler.get_sample(idx)
            raw[idx] = np.concatenate((ics, pars, row_data))

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

    errors_path = os.path.join(sampler_run_dir, f"errors_{name}.txt")

    error_lines = []

    for key, txt in reps:
        if key in errors_data:
            error_lines.append("{:55s} {}\n".format(txt, errors_data[key]))

    keys_3pf = [k for k in errors_data.keys() if "3pf" in k]

    for k3 in keys_3pf:
        for pattern, rep in reps[::-1]:

            if k3.startswith(pattern):
                txt = k3.replace(pattern, rep)

                error_lines.append("{:55s} {}\n".format(txt, errors_data[k3]))

    efficiency = round(100 * len(res_getdist) / len(raw), 2)

    error_lines.append("\n-> {:55s} {}%\n".format("Sampling efficiency", efficiency))

    with open(errors_path, "w") as f:
        f.writelines(error_lines)


def main(pool, args_dict: dict):
    n_samples = args_dict['n_samples']

    indices = np.arange(n_samples)

    print_out("Computing background trajectories")
    back_status = list(pool.map(pyt_methods.compute_background, indices))
    print_out("Background complete.")

    # exit code ok
    status_0s = list(filter(lambda x: x[1] == 0, back_status))
    status_0s = [item[0] for item in status_0s]

    print_out("Computing epsilon data")
    list(pool.map(pyt_methods.compute_epsilon, status_0s))
    print_out("Epsilon complete.")

    print_out("Computing eta data")
    list(pool.map(pyt_methods.compute_eta, status_0s))
    print_out("Eta complete.")

    print_out("Computing mass data")
    list(pool.map(pyt_methods.compute_mij, status_0s))
    print_out("Masses complete.")

    # Gather successful trajectories
    obs_pool = build_obs_pool(status_0s, args_dict)

    print_out("Computing observables")
    list(pool.map(pyt_methods.compute_obs, obs_pool))
    print_out("Observables complete.")

    pool.close()

    print("Writing results")
    write_results(args_dict, indices)
    print("All done.")
