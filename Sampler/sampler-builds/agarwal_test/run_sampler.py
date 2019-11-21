# Import numerical tools
import numpy as np
from scipy.interpolate import UnivariateSpline


# Import system tools
import os
import sys
import importlib
import pickle as pk
import time
import warnings


# Load configuration file, set PyTransport paths and import relevant modules
import config as cfg
pytpath = cfg.system['pytpath']  # PyTransport installation
saveloc = cfg.system['saveloc']  # Save location for sampler outputs
smppath = cfg.system['smppath']  # Sampler files directory


# Export environment variables for external functions
os.environ['PyTS_pytpath'] = pytpath
os.environ['PyTS_saveloc'] = saveloc
os.environ['PyTS_smppath'] = smppath
os.environ['PyTS_logpath'] = os.path.join(saveloc, "sampler_stats")


# Set Paths
rootpath = os.getcwd()
os.environ['PyTSamplerRoot'] = rootpath


path_2pf = os.path.join(saveloc, "2pf")
path_3pf = os.path.join(saveloc, "3pf")
Mij_path = os.path.join(saveloc, "Mij")
cfg_path = os.path.join(saveloc, "config.py")
smp_path = os.path.join(saveloc, "samples")


# Add methods to system path
sys.path.append(os.path.join(smppath, "sampler-methods"))


# Import PyTransport and writer tools
import pyt_methods as pytm
import writer as w


# Set PyTransport internal paths
sys.path.append(pytpath)
import PyTransSetup
PyTransSetup.pathSet()


# Load PyTransport module installation and scripts
PyT = importlib.import_module(cfg.sampler['PyTransportModule'])
import PyTransScripts as PyS



def main(pool, use_existing_samples, rerun_samples):

    print "\n... Beginning main...\n"

    # Determine tasks for sampler run
    tasks = []

    # Get 2pt tasks
    compute_2pf = cfg.computations['2pf']
    if compute_2pf: tasks.append('2pf')

    # Get 3pt tasks
    compute_3pf = cfg.computations['3pf']
    if compute_3pf:
        tasks.append('3pf')
        which3pf = cfg.which3pf

    # Get mass-matrix tasks
    compute_Mij = cfg.computations['Mij']
    if compute_Mij:
        tasks.append('Mij')

    # Get macro-parameters for sampler run
    nsamples = cfg.sampler['NSamples']
    pytmod   = cfg.sampler['PyTransportModule']
    
    # Print sampler configuration
    if sum([compute_2pf, compute_3pf, compute_Mij])==0:
        print "| Begin sampling routing: {}".format(pytmod)
        print "| -- Background only"

    else:
        print "| Begin sampling routing: {}".format(pytmod)
        print "| -- 2pf: {}".format(compute_2pf)
        print "| -- 3pf: {}".format(compute_3pf)
        if compute_3pf:
            for item in which3pf:
                print "| --------fNL {}".format(item['config_name'])
        print "| -- Mij: {}".format(compute_Mij)

    # Get number of fields and model parameters
    nF, nP = PyT.nF(), PyT.nP()


    print "\n-- Initializing ensemble"
    
    
    # Run samples with their associated background (and params)
    if use_existing_samples is True and rerun_samples is False:
        print "----- Loading ensemble"
        sample_loc = os.path.join(cfg.system['saveloc'], "samples")
        sample_range = [int(os.path.splitext(item)[0]) for item in os.listdir(sample_loc)]
    
    # Recompute samples from scratch (maintaining params only)
    elif use_existing_samples is True and rerun_samples is True:
        print "----- Re-running ensemble"
        sample_loc = os.path.join(cfg.system['saveloc'], "samples")
        sample_range = [int(os.path.splitext(item)[0]) for item in os.listdir(sample_loc)]
        pool.map(pytm.DemandSample_rerun, sample_range)
    
    # Build ensemble from scratch
    else:
        print "----- Building ensemble"
        sample_range = range(nsamples)
        pool.map(pytm.DemandSample, sample_range)

    # Generate task pool from specs. in configuration file
    taskpool = []
    
    # Each task is defined via a model number and key for task execution
    for i in sample_range:

        if compute_2pf is True:
            task = i, "2pf"
            taskpool.append(task)

        if compute_Mij is True:
            task = i, "Mij"
            taskpool.append(task)

        if compute_3pf is True:
            for config in which3pf:
                task = i, config['config_name']
                taskpool.append(task)

    # Map task pool to computation handler
    pool.map(pytm.computations, taskpool)

    # Update samples with results from task pool
    print "\n-- Updating sample objects\n"
    pool.map(w.update_samples, sample_range)

    # Close pool
    pool.close()

    # Write results file(s)
    print "\n-- Writing result files\n"
    w.write_results(nF)
    print "\n-- All Complete.\n"


if __name__ == "__main__":

    # Import schwimmbad for MPI pool processes
    import schwimmbad

    # Import argument parser for cmd args
    from argparse import ArgumentParser
    
    # Build argument parser for cmd line
    parser = ArgumentParser(description="Configure PyTransport Sampler Routine")
    
    # Build mutually exclusive group for Schwimmbad processes
    sb_group = parser.add_mutually_exclusive_group()

    sb_group.add_argument("--ncores", dest="n_cores", default=1, type=int,
                       help="Number of processes (uses multiprocessing).")
    
    sb_group.add_argument("--mpi", dest="mpi", default=False,
                       action="store_true", help="Run with MPI.")
    
    # Add group for PyTransport sampler
    pyt_group = parser.add_argument_group()
    
    # Flag to rerun computations from existing saved samples
    pyt_group.add_argument("--use_samples", dest="use_samples", default=False,
                       action="store_true", help="Use existing sample data")

    # Flag to rerun samples with their associated parameters from scratch
    pyt_group.add_argument("--rerun_samples", dest="rerun_samples", default=False,
                       action="store_true", help="Re-run with existing sample core data")
    
    args = parser.parse_args()
    
    # Build pool via schwimmbad
    pool = schwimmbad.choose_pool(mpi=args.mpi, processes=args.n_cores)

    # Execute main
    main(pool, args.use_samples, args.rerun_samples)
