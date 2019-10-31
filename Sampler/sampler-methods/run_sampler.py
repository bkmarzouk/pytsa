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


# Set Paths
rootpath = os.getcwd(); os.environ['PyTSamplerRoot'] = rootpath
path_2pf = os.path.join(saveloc, "2pf")
path_3pf = os.path.join(saveloc, "3pf")
Mij_path = os.path.join(saveloc, "Mij")
cfg_path = os.path.join(saveloc, "config.py")
smp_path = os.path.join(saveloc, "samples")

# Add methods to system path
sys.path.append(os.path.join(smppath, "sampler-methods"))

# Get additional tools
import pyt_methods as pytm
import writer as w

# Set PyTransport internal paths
sys.path.append(pytpath)
import PyTransSetup
PyTransSetup.pathSet()

# Load PyTransport module installation and scripts
PyT = importlib.import_module(cfg.sampler['PyTransportModule'])
import PyTransScripts as PyS

# Begin main function
def main(pool, use_existing_samples):

    print "\n... Beginning main...\n"

    # Determine which computations
    tasks = []

    compute_2pf = cfg.computations['2pf']
    if compute_2pf: tasks.append('2pf')

    compute_3pf = cfg.computations['3pf']
    if compute_3pf:
        tasks.append('3pf')
        which3pf = cfg.which3pf

    compute_Mij = cfg.computations['Mij']
    if compute_Mij: tasks.append('Mij')

    # Get sampler confiuration
    nsamples = cfg.sampler['NSamples']
    pytmod   = cfg.sampler['PyTransportModule']
    
    if sum([compute_2pf, compute_3pf, compute_Mij])==0:
        print "| Begin sampling routing: {}".format(pytmod)
        print "| -- Background only"

    else:
        print "| Begin sampling routing: {}".format(pytmod)
        print "| -- 2pf: {}".format(compute_2pf)
        print "| -- 3pf: {}".format(compute_3pf)
        if compute_3pf:
            for item in which3pf:
                print "|         fNL {}".format(item['config_name'])
        print "| -- Mij: {}".format(compute_Mij)

    nF, nP = PyT.nF(), PyT.nP()

    # Generate ensemble of models that support sufficient inflation
    print "\n-- Initializing ensemble\n"
    
    # If use existing sample flag, get sample numbers from sample directory
    if use_existing_samples is True:
        print "----- Loading ensemble"
        sample_loc = os.path.join(cfg.system['saveloc'], "samples")
        sample_range = [int(os.path.splitext(item)[0]) for item in os.listdir(sample_loc)]
    
    else:
        sample_range = range(nsamples)
        pool.map(pytm.DemandSample, sample_range)

    # Generate ensemble of tasks based on user specifications
    taskpool = []
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

    pool.map(pytm.computations, taskpool)

    print "\n-- Updating sample objects\n"
    pool.map(w.update_samples, sample_range)

    pool.close()

    print "\n-- Writing result files\n"
    w.write_results(nF)
    print "\n-- All Complete.\n"


# MPI Pool settings
if __name__ == "__main__":

    # Import schwimmbad for MPI pool processes
    import schwimmbad

    # Import argument parser for cmd args
    from argparse import ArgumentParser
    sb_parser = ArgumentParser(description="")

    # Add MPI standard args
    group = sb_parser.add_mutually_exclusive_group()
    group.add_argument("--ncores", dest="n_cores", default=1, type=int,
                       help="Number of processes (uses multiprocessing).")
    group.add_argument("--mpi", dest="mpi", default=False,
                       action="store_true", help="Run with MPI.")
    
    # Add flag for existing samples use
    group.add_argument("--use_samples", dest="use_samples", default=False,
                       action="store_true", help="Use existing sample data")
    
    args = sb_parser.parse_args()

    pool = schwimmbad.choose_pool(mpi=args.mpi, processes=args.n_cores)

    # Execute main program
    main(pool, args.use_samples)
