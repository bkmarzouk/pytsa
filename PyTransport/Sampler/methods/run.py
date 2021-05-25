import argparse
import os
import dill as pk

from PyTransport.Sampler.configs.setup_sampler import build_catalogue, Setup  # Import Setup for dill handling
from PyTransport.Sampler.methods import pyt_methods


class _Data(object):

    def __init__(self):
        pass


def main(pool, sampler, args):
    n_samples = args.n_samples

    # Determine tasks for sampler run
    tasks = []

    build_catalogue(sampler, args)

    # Print configuration to command line
    print("-- Starting sampling routine: {}".format(sampler.PyT))

    pool_data = []

    catalogue_path = os.path.join(sampler.path, "catalogue.npy")

    for idx in range(n_samples):
        pd = _Data()
        pd.PyT = sampler.PyT.__name__
        pd.model_number = idx
        pd.path = catalogue_path
        pool_data.append(pd)

    pool.map(pyt_methods.compute_background, pool_data)


#
#
# # write summary stats for background
# if use_samples is True:
#     if not os.path.exists(os.path.join(pathStats, "summary_bg.pk")):
#         w.bg_summary()
# else:
#     w.bg_summary()
#
# # Generate task pool(s) from specs. in configuration file
# taskpools = {'masses': []} if use_samples is False else {}  # Horizon crossing masses (Generated automatically)
# if run_2pf: taskpools['2pf'] = []  # 2pt function data
# if run_3pf:  # 3pt function data
#     for d in fNLDict: taskpools[d['name']] = []  # -> 3pt configurations
#
# # Each task is defined via a model number and key for task execution
# for i in sample_range:
#     for key in taskpools:
#         task = i, key
#         taskpools[key].append(task)
#
# # Now we map individual task pools to processors, updating the sample objects on the fly
# for key in taskpools:
#     print
#     "\n\n-- START ensemble tasks: {}\n\n".format(key)
#     pool.map(pytm.computations, taskpools[key])
#     print
#     "\n\n-- END ensemble tasks: {}\n\n".format(key)
#
# # Close pool
# pool.close()
#
# # Write results file(s)
# w.write_error_report()
# w.write_results(use_samples)
# print
# "\n\n-- All Complete.\n\n"

if __name__ == "__main__":
    # Import schwimmbad for MPI pool processes
    import schwimmbad

    # Import argument parser for cmd args
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Configure PyTransport Sampler Routine")

    # Build mutually exclusive group for schwimmbad processes
    sb_group = parser.add_mutually_exclusive_group()

    # Number of cores to use
    sb_group.add_argument("--ncores", dest="n_cores", default=1, type=int,
                          help="Number of processes (uses multiprocessing).")

    # Use mpi (required by schwimmbad routine)
    sb_group.add_argument("--mpi", dest="mpi", default=False,
                          action="store_true", help="Run with MPI.")

    sampler_group = parser.add_argument_group()

    sampler_group.add_argument("--n_samples", type=int, dest="n_samples", required=True)
    sampler_group.add_argument("--seed", type=int, dest="seed", required=False, default=None)
    sampler_group.add_argument("--grid_seed", type=int, dest="grid_seed", required=False, default=None)

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--apriori", action="store_true", dest="apriori")
    mode.add_argument("--latin", action="store_true", dest="latin")

    args = parser.parse_args()

    pk_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "sampler.pk"))

    assert os.path.exists(pk_path), "sampler file not found! {}".format(pk_path)

    with open(pk_path, "rb") as f:
        s = pk.load(f)

    pool = schwimmbad.choose_pool(mpi=args.mpi, processes=args.n_cores)

    main(pool, s, args)
