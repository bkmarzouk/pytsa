import argparse
import os
import dill as pk
import numpy as np

from PyTransport.Sampler.configs.setup_sampler import Setup
from PyTransport.Sampler.methods import samplers
from PyTransport.cache_tools import hash_pars

parser = argparse.ArgumentParser()

parser.add_argument("--n_samples", type=int, dest="n_samples", required=True)
parser.add_argument("--seed", type=int, dest="seed", required=False, default=None)
parser.add_argument("--grid_seed", type=int, dest="grid_seed", required=False, default=None)

mode = parser.add_mutually_exclusive_group(required=True)
mode.add_argument("--apriori", action="store_true", dest="apriori")
mode.add_argument("--latin", action="store_true", dest="latin")

args = parser.parse_args()

pk_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "sampler.pk"))

assert os.path.exists(pk_path), "sampler file not found! {}".format(pk_path)

with open(pk_path, "rb") as f:
    _s = pk.load(f)
    assert isinstance(_s, Setup), "{} should be a Setup instance for the PyTransport sampler".format(_s)


def build_catalogue(s: Setup, parsed_args):
    pars = [s.N_min, s.N_adiabitc, s.N_sub_evo, s.tols]
    pars += [parsed_args.seed, parsed_args.grid_seed, parsed_args.n_samples]

    hash = hash_pars(*pars)

    dir_name = "latin_" if parsed_args.latin else "apriori_"

    dir_name += hash

    dir_path = s.path

    path = os.path.join(dir_path, dir_name)

    if not os.path.exists(path):
        os.makedirs(path)

    if parsed_args.apriori:
        x = samplers.APriori(parsed_args.seed)
    else:
        x = samplers.LatinHypercube(parsed_args.n_samples, seed=parsed_args.seed, cube_seed=parsed_args.grid_seed)

    samples_path = os.path.join(path, "catalogue.npy")

    if not os.path.exists(samples_path):

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
                f = np.array([*row[1:1 + s.nF]], dtype=float)
                p = np.array([*row[-s.nP:]], dtype=float)
                V = s.PyT.V(f, p)
                dV = s.PyT.dV(f, p)

                sr_ic = -dV / np.sqrt(3 * V)

                for idx in range(s.nF):
                    if samples[row_idx][s.nF + idx] == "sr":
                        samples[row_idx][s.nF + idx] = sr_ic[idx]

        np.save(samples_path, np.asarray(samples, dtype=float))

    print("-- Sample ICs & Params @ {}".format(samples_path))
#
#
# build_catalogue(_s, args)
