import os
from pytsa.sampler import main_pool
from pytsa.sampler.setup_sampler import job_config

if __name__ == "__main__":
    os.environ['SAMPLER_MODE'] = "TRUE"

    # Import schwimmbad for MPI pool processes
    import schwimmbad
    from pytsa.sampler.mpi_helpers import size as n_proc

    # Import argument parser for cmd args
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Configure PyTransport sampler Routine")

    parser.add_argument("name", metavar="name", type=str,
                        help="Name for sampling routine: defines subdir")

    sampler_group = parser.add_argument_group()

    sampler_group.add_argument("--n_samples", type=int, dest="n_samples", required=True,
                               help="Number of samples to compute")
    sampler_group.add_argument("--entropy", type=int, dest="entropy", required=False, default=None,
                               help="System entropy value: Enables reproducible results")

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--apriori", action="store_true", dest="apriori",
                      help="If flagged, then run sampler in apriori mode")
    mode.add_argument("--latin", action="store_true", dest="latin",
                      help="If flagged, then run sampler in latin hypercube mode")

    observables_group = parser.add_argument_group()

    sampler_group.add_argument("--ns", action="store_true", dest="task_2pt", help="If flagged, computes 2-point data")
    sampler_group.add_argument("--eq", action="store_true", dest="task_3pt_eq",
                               help="If flagged, computes 3-point data in equilateral mode configuration")
    sampler_group.add_argument("--fo", action="store_true", dest="task_3pt_fo",
                               help="If flagged, computes 3-point data in folded mode configuration")
    sampler_group.add_argument("--sq", action="store_true", dest="task_3pt_sq",
                               help="If flagged, computes 3-point data in squeezed mode configutarion")

    extra_group = parser.add_argument_group()

    extra_group.add_argument("--alpha", nargs="+", required=False,
                             help="Pass values of alpha (Fergusson Shellard convention) to compute custom fnls. "
                                  "NOTE: Requires equal number of beta definitions.")
    extra_group.add_argument("--beta", nargs="+", required=False,
                             help="Pass values of beta (Fergusson Shellard convention) to compute custom fnls. "
                                  "NOTE: Requires equal number of alpha definitions.")

    parsed_args = parser.parse_args()

    args_dict = vars(parsed_args)

    for k in ['alpha', 'beta']:
        if args_dict[k] is None:
            args_dict[k] = []

    assert len(parsed_args.alpha) == len(parsed_args.beta), [parsed_args.alpha, parsed_args.beta]

    args_dict['cwd'] = os.path.abspath(os.path.dirname(__file__))

    job_config(args_dict)

    pool = schwimmbad.choose_pool(mpi=n_proc > 1, processes=n_proc, use_dill=True)

    main_pool.main(pool, args_dict)
