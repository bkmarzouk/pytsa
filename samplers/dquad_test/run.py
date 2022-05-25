import os
import dill as dill
import pickle as pk

from pytransport.sampler.configs.setup_sampler import SamplerMethods, APrioriSampler, LatinSampler
from pytransport.sampler.methods import main_pool
from pytransport.sampler.configs.rng_states import RandomStates


def job_config(task_pars: dict):
    name = task_pars['name']
    n_samples = task_pars['n_samples']
    entropy = task_pars['entropy']
    apriori = task_pars['apriori']

    n_states = n_samples if apriori else n_samples + 1

    root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), name))

    if not os.path.exists(root_dir):
        os.makedirs(root_dir)

    tasks_path = os.path.join(root_dir, "sampler.tasks")

    if os.path.exists(tasks_path):

        with open(tasks_path, "rb") as f:

            cached_task_pars: dict = pk.load(f)

        assert cached_task_pars.keys() == task_pars.keys(), [cached_task_pars.keys(), task_pars.keys()]

        for k in cached_task_pars.keys():
            parsed_par = task_pars[k]
            cached_par = cached_task_pars[k]

            assert parsed_par == cached_par, f"Parameter error @ {k}: {parsed_par} != {cached_par}"

    else:

        with open(tasks_path, "wb") as f:

            pk.dump(task_pars, f)

    random_states = RandomStates.from_cache(root_dir, entropy=entropy, n_states=n_states)

    sampler_path = os.path.join(root_dir, "sampler.run")

    if not os.path.exists(sampler_path):
        methods_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "sampler.methods"))

        with open(methods_path, "rb") as f:
            sampler_methods = dill.load(f)  # DILL??

        _m = APrioriSampler if apriori else LatinSampler

        sampler_run = _m(random_states, sampler_methods)

        with open(sampler_path, "wb") as f:
            dill.dump(sampler_run, f)


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

    parsed_args = parser.parse_args()

    args_dict = vars(parsed_args)

    job_config(args_dict)

    assert 0, args_dict

    pool = schwimmbad.choose_pool(mpi=args.n_procs > 1, processes=args.n_procs)

    main_pool.main(pool, s, args)
