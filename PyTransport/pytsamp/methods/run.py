import os
import dill as dill

from PyTransport.sampler.configs import setup_sampler
from PyTransport.sampler.methods import main_pool

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

    path = os.path.abspath(os.path.join(os.path.dirname(__file__), "sampler"))

    assert os.path.exists(path), "sampler file not found! {}".format(path)

    with open(path, "rb") as f:
        s = dill.load(f)

    pool = schwimmbad.choose_pool(mpi=args.mpi, processes=args.n_cores)

    main_pool.main(pool, s, args)
