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

# Set Paths
rootpath = os.getcwd(); os.environ['PyTSamplerRoot'] = rootpath
path_2pf = os.path.join(rootpath, "2pf")
path_3pf = os.path.join(rootpath, "3pf")
Mij_path = os.path.join(rootpath, "Mij")
cfg_path = os.path.join(rootpath, "config.py")
smp_path = os.path.join(rootpath, "samples")
tool_path = os.path.abspath(os.path.join(rootpath, "..", "..", "samplerfiles"))
sys.path.append(tool_path)

import pyt-methods as oscript
import Writer as w

# Load configuration file, set PyTransport paths and import relevant modules
import config as cfg
sys.path.append(cfg.pytpath)
import PyTransSetup
PyTransSetup.pathSet()

PyT = importlib.import_module(cfg.sampler['PyTransportModule'])
import PyTransScripts as PyS


# Begin main function
def main(pool):

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


    assert compute_2pf or compute_3pf or compute_Mij, "No calculations found."

    # Get sampler confiuration
    nsamples = cfg.sampler['NSamples']
    pytmod   = cfg.sampler['PyTransportModule']

    print "| Begin sampling routing: {}".format(pytmod)
    print "| -- 2pf: {}".format(compute_2pf)
    print "| -- 3pf: {}".format(compute_3pf)
    if compute_3pf:
        for item in which3pf:
            print "|         fNL {}".format(item['config_name'])
    print "| -- Mij: {}".format(compute_Mij)

    nF, nP = PyT.nF(), PyT.nP()
    #
    # # Setup writer for results
    # result_cols = ()
    # if compute_2pf:
    #     result_cols+=("ns", "alpha", "2pf_t")
    #
    # if compute_3pf:
    #     for item in which3pf:result_cols+=(item['config_name'], "3pf_{}_t".format(item['config_name']))
    #
    # if compute_Mij:
    #     for f in range(nF): result_cols+=("m_{}".format(f),)
    #
    # for p in cfg.parameter_values:
    #     pname = p['ParameterName']
    #     if pname is not None:result_cols+=(pname,)


    # Generate ensemble of models that support sufficient inflation
    print "\n-- Initializing ensemble\n"
    pool.map(oscript.DemandSample, range(nsamples))

    # Generate ensemble of tasks based on user specifications
    taskpool = []
    for i in range(nsamples):

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

    pool.map(oscript.computations, taskpool)

    print "\n-- Updating sample objects\n"
    pool.map(w.update_samples, range(nsamples))

    pool.close()

    print "\n-- Writing result files\n"
    w.write_results(nF)

    #
    # print '\n-- Beginning MapReduce job: Background evolution\n'
    #
    #
    # # Calculate background cosmology for samples
    # BGcosmo = pool.map(oscript.ModelSetup, samples)
    #
    # print '\n-- finished MapReduce job: Background evolution\n'
    #
    # # For active observable types, build job pool
    # jobpool = []
    # for item in BGcosmo:
    #
    #     if item is not None:
    #
    #         if ns_active:
    #             job = item, "ns"
    #             jobpool.append(job)
    #
    #         for s in shapes:
    #             job = item, s
    #             jobpool.append(job)
    #
    #
    # # MapReduce observables
    # print '\n-- beginning MapReduce job: Observables\n'
    #
    # # Perform MapReduce on observables.
    # # Note: We don't define variable, as all data is saved as pickle file and converted to results.txt file
    # pool.map(taskhandler.observables, jobpool)
    # pool.close()
    #
    # print '\n-- finished MapReduce job: Observables\n'
    #
    # print '\n-- begin writing phase\n'
    #
    # # Call writer to write initial ensemble file and results in GetDist format
    # writer.writegetdist()
    #
    # print '\n-- writing phase complete.\n'
    #
    # print '\n-- Job complete.'


# MPI Pool settings
if __name__ == "__main__":

    import schwimmbad

    # Import argument parser for handling command line arguments, setup parser for Schwimmbad
    from argparse import ArgumentParser
    sb_parser = ArgumentParser(description="")

    # We define the mutually exclusive group for MPI standard procedures
    group = sb_parser.add_mutually_exclusive_group()
    group.add_argument("--ncores", dest="n_cores", default=1, type=int,
                       help="Number of processes (uses multiprocessing).")
    group.add_argument("--mpi", dest="mpi", default=False, action="store_true", help="Run with MPI.")

    args = sb_parser.parse_args()

    pool = schwimmbad.choose_pool(mpi=args.mpi, processes=args.n_cores)

    # Execute main program
    main(pool)
