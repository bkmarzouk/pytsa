# Import system tools
import os
import sys
import pickle as pk
import csv
import shutil
import numpy as np

# Retrieve paths from environment variables
envKeys = [
    'PyTS_pathPyT', 'PyTS_pathRoot', 'PyTS_pathSamples',
    'PyTS_pathStats', 'PyTS_pathClasses', 'PyTS_pathMethods', 'PyTS_pathLocalData'
]

pathPyT, pathRoot, \
pathSamples, pathStats, pathClasses, pathMethods, pathLocalData = [os.environ[eK] for eK in envKeys]


# Get fNL data
fNLPath = os.path.join(pathLocalData, "fNL.localdata")
fNLFile = open(fNLPath, "rb")
with fNLFile: fNLDict = pk.load(fNLFile)


def load_obs(pk_path):
    f = open(pk_path, "rb")
    with f: obs = pk.load(f)
    return obs


def load_pk(pk_path):
    """ Wrapper for safe-load of binary files """
    f = open(pk_path, "rb")
    with f: obj = pk.load(f)
    return obj


def rewrite_pk(pk_obj, pk_path):
    """ Wrapper for safe-rewrite of binary files """
    assert os.path.exists(pk_path), "Cannot rewrite non-existing file: {}".format(pk_path)
    os.remove(pk_path)
    f = open(pk_path, "wb")
    with f: pk.dump(pk_obj, pk_path)


def update_sample(modelnumber, obsDict, timerSubscript):
    
    pathFullSample = os.path.join(pathSamples, "{}.sample".format(modelnumber))
    
    sample = load_pk(pathFullSample)
    
    assert sample.modelnumber == modelnumber
    
    sample.update_observables(obsDict, timerSubscript)
    
    rewrite_pk(sample, pathFullSample)





def bg_summary():
    
    # Build paths to pickled minor stats
    bgd = os.path.join(pathStats, "bg") # bg dir
    
    mdir = os.path.join(pathStats, "mij")
    
    if not os.path.exists(mdir):
        os.makedirs(mdir)
    
    stat_paths = [os.path.join(bgd, item) for item in os.listdir(bgd)]
    
    # Set counters dictionary (will sum minor stats)
    # Flag defs.
    timeoutFlags = [
        -10, -11, -12, -13  # end of inflation, background, 2pf, 3pf
    ]

    integratorFlags = [
        -20, -21, -22, -23  # end of inflation, background, 2pf, 3pf
    ]

    samplerFlags = [
        -30, -31, -32, -33, -34, -35, -36  # N < Nmin, k not found, no ICs 2pf, no ICs 3pf, model violation, eternal, mij
    ]

    allFlags = timeoutFlags + integratorFlags + samplerFlags
    allKeys = [str(f) for f in allFlags] + ["time"]
    totDict = {k: 0 for k in allKeys}
    
    print "\n\n-- Writing background statistics\n\n"
    
    for sp in stat_paths:
    
        # unload binary file
        f = open(sp, "rb")
        with f:
            try:
                stats = pk.load(f)
            except:
                raise OSError, sp
        
        # add corresponding dictionary keys
        for k in allKeys:

            try:

                totDict[k] += stats[k]
    
            except KeyError:

                if k == 'time':
                    print 'Sample data unavailable since no "time key": {}'.format(sp)

                else:
                    raise KeyError, "key={} not in {}".format(k, sp)

    for sp in stat_paths:
        # remove minor stats data
        os.remove(sp)
        
    # Write single
    f = open(os.path.join(pathStats, "summary_bg.pk"), "wb")
    
    with f: pk.dump(totDict, f)
    
    
def write_error_report():
    
    # Build paths to pickled minor stats
    mijDir = os.path.join(pathStats, "mij")
    twoptDir = os.path.join(pathStats, "2pf")
    threeptDir = os.path.join(pathStats, "3pf")

    mij_paths = []
    twopt_paths = []
    threept_paths = []

    if os.path.exists(mijDir):    
        mij_paths = [os.path.join(mijDir, p) for p in os.listdir(mijDir)]
    if os.path.exists(twoptDir):
        twopt_paths = [os.path.join(twoptDir, p) for p in os.listdir(twoptDir)]
    if os.path.exists(threeptDir):
        threept_paths = [os.path.join(threeptDir, p) for p in os.listdir(threeptDir)]
    
    # Flag defs.
    
    timeoutFlags = [
        -10, -11, -12, -13  # end of inflation, background, 2pf, 3pf
    ]
    
    integratorFlags = [
        -20, -21, -22, -23  # end of inflation, background, 2pf, 3pf
    ]
    
    samplerFlags = [
        -30, -31, -32, -33, -34, -35, -36  # N < Nmin, k not found, no ICs 2pf, no ICs 3pf, model violation, eternal
    ]
    
    allFlags = timeoutFlags + integratorFlags + samplerFlags
    allKeys = [str(f) for f in allFlags]
    
    sumBG = os.path.join(pathStats, "summary_bg.pk")
    sumObs = os.path.join(pathStats, "summary_obs.pk")

    if os.path.exists(sumObs):
        
        with open(sumObs, "rb") as f:
            
            totDict = pk.load(f)
    
    else:
        totDict = {k: 0 for k in allKeys}


    # TODO: Report will be wrong in instances of rerunning observables
    # since error calculations will be additive to previous run !
    # Find some means of tracking these errors better ...


    # Combine 2pt error stats
    for p in mij_paths + twopt_paths + threept_paths:
        
        f = open(p, "rb")
        with f: s = pk.load(f)
        
        k = str(s['flag'])
        totDict[k] += 1

    for p in mij_paths + twopt_paths + threept_paths:
        os.remove(p)
    
    f = open(sumObs, "wb")
    with f: pk.dump(totDict, f)

    g = open(sumBG, "rb")
    with g:
        bgStats = pk.load(g)
    
    bgFailCounter = 0
    totFailCounter = 0
    
    n_obtained = len(os.listdir(pathSamples))
    
    for k in allKeys:
        if k in bgStats:
            bgFailCounter += bgStats[k]
            totDict[k]    += bgStats[k]
        if k in totDict:
            totFailCounter += totDict[k]
    
    averageSampleTime = float(bgStats['time']) / float(n_obtained)
    
    sampleEfficiency = float(n_obtained) / float(n_obtained + bgFailCounter)
    
    lines = [
    
        "---- Sample rejection summary:\n\n",
        
        "-- Integration timeout:\n"
        "Find end of inflation:             {}\n".format(totDict['-10']),
        "Background evolution:              {}\n".format(totDict['-11']),
        "Two-point function:                {}\n".format(totDict['-12']),
        "Three-point function:              {}\n\n".format(totDict['-13']),
        
        "-- Integration error:\n"
        "Find end of inflation:             {}\n".format(totDict['-20']),
        "Background evolution:              {}\n".format(totDict['-21']),
        "Two-point function:                {}\n".format(totDict['-22']),
        "Three-point function:              {}\n\n".format(totDict['-23']),
        
        "-- Miscellaneous:\n",
        "Inflation too short:               {}\n".format(totDict['-30']),
        "Unable to find end of inflation:   {}\n".format(totDict['-35']),
        "Unable to find momenta (all):      {}\n".format(totDict['-31']),
        "Unable to find ICs (2pf):          {}\n".format(totDict['-32']),
        "Unable to find ICs (3pf):          {}\n".format(totDict['-33']),
        "Model violation:                   {}\n".format(totDict['-34']),
        "Eigs linalg. error:                {}\n\n".format(totDict['-36']),
        
        "Total number of rejections:        {}\n\n".format(totFailCounter),
        
        "---- Probability of inflating (background-level):\n\n",
        
        "Sample efficiency:                 {}\n".format(sampleEfficiency),
        "Total time to build ensemble:      {}\n".format(bgStats['time']),
        "Average time to obtain sample:     {}\n".format(averageSampleTime)
        
    ]

    summary_file = open(os.path.join(pathRoot, "outputs", "sampling_summary.txt"), "w")
    
    with summary_file:
        summary_file.writelines(lines)


def compress_samples(undo=False):
    """ Wraps all sample objects into single-large binary file.
    Large storage is usually preferred on scratch systems. """
    
    # If undo is True, the reverse procedure will be performed
    if undo is False:
        
        # Get full path definitions for samples
        sample_paths = [
            os.path.join(pathSamples, item) for item in os.listdir(pathSamples)
        ]
        
        # Open file to dump all samples to single binary object
        ensemblePickles = open(os.path.join(pathSamples, "ensemble.samples"), "wb")
        
        # With file handler
        with ensemblePickles:
            # For each item in the sample paths
            for item in sample_paths:
                # Open the sample file
                samplePickle = open(item, "rb")
                
                # Load the sample
                with samplePickle:
                    sample = pk.load(samplePickle)
                
                # Dump to the ensemble
                pk.dump(sample, ensemblePickles)
                
                # Remove single sample object
                os.remove(item)

    else:
        
        ensemblePaths = os.listdir(pathSamples)
        
        if len(ensemblePaths) == 1 and ensemblePaths[0] == "ensemble.samples":
            
            # open ensemble pickle file
            ensemblePickles = open(os.path.join(pathSamples, "ensemble.samples"), "rb")
            
            # Start sample count
            sampleCounter = 0
            
            try:
                
                print "\n\n-- Unpacking ensemble samples\n\n"
                
                with ensemblePickles:
                    
                    # Keep unloading files until an EOFError is raise
                    while True:
                        
                        s = pk.load(ensemblePickles)
                        
                        modelnumber = s.modelnumber
                        
                        # Write with standard syntax as before
                        samplePickle = open(os.path.join(pathSamples, "{}.sample").format(modelnumber), "wb")
                        
                        with samplePickle: pk.dump(s, samplePickle)
                        
                        sampleCounter += 1
            
            except EOFError: # Raised when ran out of binary objects in file
                print "\n\n-- Unpack complete: {} samples\n\n".format(sampleCounter)
                os.remove(os.path.join(pathSamples, "ensemble.samples"))

        elif len(os.listdir(pathSamples)) == 0:  # No samples found
            raise OSError, "No samples found in {}".format(pathSamples)
        
        elif np.all([item.endswith(".sample") for item in os.listdir(pathSamples)]): # If previous run didn't reach unpack stage
            print "\n\n-- Samples already uncompressed \n\n"
            
        else: # Something else...
            raise OSError, "Unrecognized files in {}".format(pathSamples)

def write_results(rerun=False):

    headers = ('weight', 'like')
    paths   = [ os.path.join(pathSamples, m) for m in os.listdir(pathSamples)]
    labels  = []
    latexs  = []

    ncols = 2
    
    latexFile = open(os.path.join(pathLocalData, "latex.localdata"), "rb")
    
    with latexFile: latexDefs = pk.load(latexFile)
    
    fkeys = []
    vkeys = []
    pkeys = []
    for k in sorted(latexDefs.keys()):
        if   k.startswith("f"): fkeys.append(k)
        elif k.startswith("v"): vkeys.append(k)
        elif k.startswith("p"): pkeys.append(k)
        else: raise KeyError, k
    
    sortedKeys = fkeys + vkeys + pkeys
    
    for key in sortedKeys:
        headers += (key,)
        labels.append(key)
        ncols += 1
        latexs.append(latexDefs[key])
    
    # The dictionary of observables associated with any gifven sample is expected to be representative
    f = open(paths[0], "rb")
    with f: s0 = pk.load(f)
    
    obsKeys = s0.observables.keys()
    
    nF = len(s0.fields)
    
    # Define all possible keys
    possibleKeys  = ["m{}".format(ii) for ii in range(nF)] + ["T_masses"]
    possibleLaTeX = ["m_{}^2/H^2".format(ii) for ii in range(nF)] + ["T_{M^A_B}"]

    possibleKeys  += ["ns", "alpha", "T_2pf"]
    possibleLaTeX += ["n_s", "\\alpha", "T_{n_s}"]
    
    fNLFile = open(os.path.join(pathLocalData, "fNL.localdata"), "rb")
    with fNLFile:
        configDefs = pk.load(fNLFile)
        
    for cD in configDefs:
        possibleKeys.append(cD['name'])
        possibleKeys.append("T_{}".format(cD['name']))
        
        possibleLaTeX.append(cD['latex'])
        possibleLaTeX.append("T_" + "{f_{NL}^{" + cD['name'] + "}}")
    
    for k in possibleKeys:
        if k in obsKeys:
            headers += (k,)
            labels.append(k)
            ncols += 1
            latexs.append(possibleLaTeX[possibleKeys.index(k)])

    # We will write 3 results files that contain all, evolving and adiabatic samples
    dicts_all       = []
    dicts_evolving  = []
    dicts_adiabatic = []
    failed_samples  = []

    # Iterate over all sample paths
    for p in paths:
        
        if rerun:
            with open(p, "rb") as f:
                s = pk.load(f)
            s.check_reject()
        
        # Open sample
        f = open(p, "rb")
        
        # With handler
        with f:
            
            # load binary file and get sample line from method
            s = pk.load(f)
            
            d = s.line_dict()

            # If the sample hasn't been rejected on account of failed computation(s)
            if s.reject is False:
                
                # Add to 'all' results
                dicts_all.append(d)
                
                # If adiabatic, add to adiabatic results
                if s.adiabatic is True:
                    dicts_adiabatic.append(d)
                # Otherwise add to evolving results
                else:
                    dicts_evolving.append(d)

            # Record failed sample model number
            else:
                failed_samples.append(p+"\n")

    p = open(os.path.join(pathRoot, "outputs", "adiabatic.txt"), "w")
    q = open(os.path.join(pathRoot, "outputs", "evolving.txt" ), "w")
    r = open(os.path.join(pathRoot, "outputs", "all.txt"), "w")

    for h_r in zip([p, q, r], [dicts_adiabatic, dicts_evolving, dicts_all]):
        
        handler, dicts_results = h_r
        
        with handler:
            w = csv.DictWriter(handler, headers, delimiter=' ')
            for d in dicts_results:
                w.writerow(d)

    for pfile in [open(os.path.join(pathRoot, "outputs", "{}.paramnames").format(fname), "w")
                for fname in ["adiabatic", "evolving", "all"]]:
        
        with pfile:
            w = csv.DictWriter(pfile, ('label', 'LaTeX'), delimiter = ' ')
            for lab_lat in zip(labels, latexs):
                lab, lat = lab_lat
                w.writerow({'label': lab, 'LaTeX': lat})
    

    pfile1 = os.path.join(pathRoot, "outputs", "adiabatic.paramnames")
    pfile2 = os.path.join(pathRoot, "outputs", "evolving.paramnames")
    h = open(pfile1, "w")
    headers = ('label', 'LaTeX')
    with h:
        x = csv.DictWriter(h, headers, delimiter=' ')
        for row in zip(labels, latexs):
            lab, lat = row
            x.writerow({'label': lab, 'LaTeX': lat})
    shutil.copyfile(pfile1, pfile2)

    i = open(os.path.join(pathRoot, "outputs", "failed.txt"), "w")
    with i:
        i.writelines(failed_samples)

    compress_samples()
