# Import system tools
import os, sys, importlib, pickle as pk

# Import numerical tools
import numpy as np
from scipy.interpolate import UnivariateSpline

# Import other useful tools
import time, warnings


# Navigate to local files
runDir = os.path.dirname(os.path.abspath(__file__))
pathLocalData = os.path.join(runDir, ".localdata")
sys.path.append(pathLocalData)
import generator


# Get environment data
envPath = os.path.join(pathLocalData, "env.localdata")
envFile = open(envPath, "rb")
with envFile: envDict = pk.load(envFile)


# Define set of requried keys, export to environment & define system paths
envKeys = [
    'PyTS_pathPyT', 'PyTS_pathRoot', 'PyTS_pathSamples',
    'PyTS_pathStats', 'PyTS_pathClasses', 'PyTS_pathMethods'
]

for eK in envKeys + ['PyTS_pathLocalData']: os.environ[eK] = envDict[eK]

pathPyT, pathRoot, \
pathSamples, pathStats, pathClasses, pathMethods = [envDict[eK] for eK in envKeys]


# Get fNL data
fNLPath = os.path.join(pathLocalData, "fNL.localdata")
fNLFile = open(fNLPath, "rb")
with fNLFile: fNLDict = pk.load(fNLFile)


# get transport data
transportPath = os.path.join(pathLocalData, "transport.localdata")
transportFile = open(transportPath, "rb")
with transportFile: transportDict = pk.load(transportFile)


# Add paths to additional tools and objects
# sys.path.append(pathClasses)
sys.path.append(pathMethods)


# Import PyTransport and writer tools
import pyt_methods as pytm
import writer as w


# Finally configure internal PyTransport workings
sys.path.append(pathPyT)
import PyTransSetup
PyTransSetup.pathSet()


# Load PyTransport module installation and scripts
PyT = importlib.import_module(transportDict['transportModule'])
import PyTransScripts as PyS


def main(pool, n_samples, run_2pf, run_3pf, use_samples):

    # Determine tasks for sampler run
    tasks = []

    # Get tasks
    if run_2pf: tasks.append("2pf")
    if run_3pf: tasks.append("3pf")
    

    # Print configuration to command line
    print "\n\n-- Starting sampling routine: {}".format(transportDict['transportModule'])
    if run_2pf: print "-- Computing 2pt data"
    if run_3pf: print "-- Computing 3pt data"


    print "\n\n\n\n-- Initializing ensemble"


    # Run existing sample data
    if use_samples is True:
        print "----- Loading ensemble\n\n"
        
        # Need to undo compression to use samples again
        w.compress_samples(undo=True)
        
        sample_range = [int(os.path.splitext(item)[0]) for item in os.listdir(pathSamples)]
        
    # Build ensemble from scratch
    else:
        print "----- Building ensemble\n\n"
        sample_range = range(n_samples)
        pool.map(pytm.DemandSample, sample_range)
        

    # write summary stats for background
    if use_samples is True:
        if not os.path.exists(os.path.join(pathStats, "summary_bg.pk")):
            w.bg_summary()
    else: w.bg_summary()


    # Generate task pool(s) from specs. in configuration file
    taskpools = {'masses': []} if use_samples is False else {}  # Horizon crossing masses (Generated automatically)
    if run_2pf: taskpools['2pf'] = []                           # 2pt function data
    if run_3pf:                                                 # 3pt function data
        for d in fNLDict: taskpools[d['name']] = []             # -> 3pt configurations

    # Each task is defined via a model number and key for task execution
    for i in sample_range:
        for key in taskpools:
            task = i, key
            taskpools[key].append(task)
    
    # Now we map individual task pools to processors, updating the sample objects on the fly
    for key in taskpools:
        
        print "\n\n-- START ensemble tasks: {}\n\n".format(key)
        pool.map(pytm.computations, taskpools[key])
        print "\n\n-- END ensemble tasks: {}\n\n".format(key)


    # Close pool
    pool.close()

    # Write results file(s)
    print "\n\n-- Writing result files\n\n"
    w.write_error_report()
    w.write_results()
    print "\n\n-- All Complete.\n\n"


if __name__ == "__main__":

    # Import schwimmbad for MPI pool processes
    import schwimmbad

    # Import argument parser for cmd args
    from argparse import ArgumentParser
    
    parser = ArgumentParser(description="Configure PyTransport Sampler Routine")

    """ Args1: Set schwimmbad mpi """

    # Build mutually exclusive group for schwimmbad processes
    sb_group = parser.add_mutually_exclusive_group()

    # Number of cores to use
    sb_group.add_argument("--ncores", dest="n_cores", default=1, type=int,
                       help="Number of processes (uses multiprocessing).")

    # Use mpi (required by schwimmbad routine)
    sb_group.add_argument("--mpi", dest="mpi", default=False,
                       action="store_true", help="Run with MPI.")

    """ Args2: Set sampling routine: Some number of samples / use existing sample data / rerun existing samples """

    # Add group for samples
    sampler_group = parser.add_mutually_exclusive_group(required=True)

    sampler_group.add_argument("--n_samples", dest="n_samples", type=int,
                        help="Specify number of samples to obtain")

    # Flag to rerun computations from existing saved samples
    sampler_group.add_argument("--use_samples", dest="use_samples", default=False,
                       action="store_true", help="Use existing sample data")

    """ Args3: Set which calculations to perform """

    # Add group for observables & timeout conditions for integrator
    obs_group = parser.add_argument_group()

    # Compute spectral index / running
    obs_group.add_argument("--run_2pf", dest="run_2pf", default=False,
                       action="store_true", help="Compute 2pt function observables")

    # Compute fNL
    obs_group.add_argument("--run_3pf", dest="run_3pf", default=False,
                       action="store_true", help="Compute 3pt function observables")

    # Max number of seconds that background tasks are allowed to run for (default 5 mins)
    obs_group.add_argument("--tmax_bg", dest="tmax_bg", default=300, type=int,
                        help="Number of seconds background tasks may run for")

    # Max number of seconds that background tasks are allowed to run for (default 10 mins)
    obs_group.add_argument("--tmax_2pf", dest="tmax_2pf", default=600, type=int,
                        help="Number of seconds 2-point function tasks may run for")

    # Max number of seconds that background tasks are allowed to run for (default 30 mins)
    obs_group.add_argument("--tmax_3pf", dest="tmax_3pf", default=1800, type=int,
                        help="Number of seconds 3-point function tasks may run for")


    """ Parse arguments"""

    args = parser.parse_args()

    # Export max execution times as environment variables
    os.environ["tmax_bg"]  = str(args.tmax_bg)
    os.environ["tmax_2pf"] = str(args.tmax_2pf)
    os.environ["tmax_3pf"] = str(args.tmax_3pf)

    # If using ensemble of samples, we want to allow dictionary key rewriting; hence we export flag
    os.environ['PyTS_rewriteObs'] = str(args.use_samples)

    # Build pool via schwimmbad
    pool = schwimmbad.choose_pool(mpi=args.mpi, processes=args.n_cores)

    # Execute main
    main(pool, args.n_samples, args.run_2pf, args.run_3pf, args.use_samples)
