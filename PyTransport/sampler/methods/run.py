import os
import dill as dill

from pytransport.sampler.configs.setup_sampler import Setup
from pytransport.sampler.methods import main_pool

if __name__ == "__main__":
    os.environ['SAMPLER_MODE'] = "TRUE"

    # Import schwimmbad for MPI pool processes
    import schwimmbad

    # Import argument parser for cmd args
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Configure PyTransport sampler Routine")

    parser.add_argument("name", metavar="name", type=str,
                        help="Name for sampling routine: defines subdir")

    parser.add_argument("--nproc", dest="n_procs", default=1, type=int,
                        help="Number of processes, intializes mpi pool")

    sampler_group = parser.add_argument_group()

    sampler_group.add_argument("--n_samples", type=int, dest="n_samples", required=True)
    sampler_group.add_argument("--entropy", type=int, dest="entropy", required=False, default=None)

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--apriori", action="store_true", dest="apriori")
    mode.add_argument("--latin", action="store_true", dest="latin")

    observables_group = parser.add_argument_group()

    # TODO: docs
    sampler_group.add_argument("--ns", action="store_true", dest="run_2pt")
    sampler_group.add_argument("--eq", action="store_true", dest="run_3pt_eq")
    sampler_group.add_argument("--fo", action="store_true", dest="run_3pt_fo")
    sampler_group.add_argument("--sq", action="store_true", dest="run_3pt_sq")

    extra_group = parser.add_argument_group()

    extra_group.add_argument("--alpha", nargs="+", help="Additional 3pt tasks: alpha definitions", required=False)
    extra_group.add_argument("--beta", nargs="+", help="Additional 3pt tasks: beta definitions", required=False)

    args = parser.parse_args()

    path = os.path.abspath(os.path.join(os.path.dirname(__file__), "sampler"))

    assert os.path.exists(path), "sampler file not found! {}".format(path)

    with open(path, "rb") as f:
        s: Setup = dill.load(f)

    assert 0, args

    pool = schwimmbad.choose_pool(mpi=args.n_procs > 1, processes=args.n_procs)

    main_pool.main(pool, s, args)
