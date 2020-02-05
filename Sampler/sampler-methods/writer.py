# Import system tools
import os
import sys
import pickle as pk
import csv
import shutil

# Retrieve paths from environment variables
envKeys = [
    'PyTS_pathPyT', 'PyTS_pathRoot', 'PyTS_path2pf', 'PYTS_path3pf', 'PyTS_pathMasses', 'PyTS_pathSamples',
    'PyTS_pathStats', 'PyTS_pathClasses', 'PyTS_pathMethods', 'PyTS_pathLocalData'
]

pathPyT, pathRoot, path2pf, path3pf, pathMasses, \
pathSamples, pathStats, pathClasses, pathMethods, pathLocalData = [os.environ[eK] for eK in envKeys]


# Get fNL data
fNLPath = os.path.join(pathLocalData, "fNL.localdata")
fNLFile = open(fNLPath, "rb")
with fNLFile: fNLDict = pk.load(fNLFile)


def load_obs(pk_path):
    f = open(pk_path, "rb")
    with f: obs = pk.load(f)
    return obs


def bg_summary():
    
    # Build paths to pickled minor stats
    bgd = os.path.join(pathStats, "bg")
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
        -30, -31, -32, -33, -34, -35  # N < Nmin, k not found, no ICs 2pf, no ICs 3pf, model violation, eternal
    ]

    allFlags = timeoutFlags + integratorFlags + samplerFlags
    allKeys = [str(f) for f in allFlags] + ["time"]
    totDict = {k: 0 for k in allKeys}
    
    print "\n-- Writing background statistics\n"
    
    rej_tot = 0
    
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
            totDict[k] += stats[k]
    
    for sp in stat_paths:
        # remove minor stats data
        os.remove(sp)
        
    f = open(os.path.join(pathStats, "bg", "summary.pk"), "rb")
    
    with f: pk.dump(totDict, f)
    
    #
    # stats_file = open(os.path.join(pathStats, "summary.txt"), "w")
    #
    # error    = "-- Rejections summary:\n"
    # short    = "Insufficient inflation:                {}\n".format(counters['short'])
    # kexit    = "Failed to compute kExit:               {}\n".format(counters['kexit'])
    # feoi     = "Integration error, findEndOfInflation: {}\n".format(counters['feoi'])
    # back     = "Integration error, backEvolve:         {}\n".format(counters['back'])
    # violated = "Model violation:                       {}\n".format(counters['violated'])
    # timeout  = "Integration timout:                    {}\n".format(counters['timeout'])
    # rejects  = "Total rejected samples:                {}\n\n".format(rej_tot)
    #
    # success  = "-- Ensemble overview:\n"
    # tsa      = "Total number of attempted samples:     {}\n".format(counters['end'])
    # tss      = "Total number of successful samples:    {}\n".format(counters['samples'])
    # p_N      = "Probability of inflation:              {}\n".format(float(counters['samples']) / float(counters['end']))
    # ttotal   = "Total time to build ensemble:          {} seconds\n".format(float(counters['time']))
    # taverage = "Av. time to obtain successful sample:  {}\n".format(float(counters['time']/ float(counters['samples'])))
    #
    # lines = [error, short, kexit, feoi, back, violated, timeout, rejects, success, tsa, tss, p_N, ttotal, taverage]
    #
    # with stats_file as f:
    #     for line in lines:
    #         f.write(line)
    #
    

def update_samples(modelnumber):

    # Load sample
    save_path = os.path.join(pathSamples, "{}.sample".format(modelnumber))
    f=open(save_path, "rb")
    with f:
        sample = pk.load(f)
    assert sample.modelnumber == modelnumber

    if len(os.listdir(path2pf)) > 0:
        twoPt_path = os.path.join(path2pf, "{}.2pf".format(modelnumber))
        obs = load_obs(twoPt_path)
        sample.update_observables(obs, "2pf")

    if len(os.listdir(path3pf)) > 0:
        for d in fNLDict:
            cname = d['name']
            cpath = os.path.join(path3pf, "{m}.{c}".format(m=modelnumber, c=cname))
            obs = load_obs(cpath)
            sample.update_observables(obs, cname)

    if len(os.listdir(pathMasses)) > 0:
        Mij_obspath = os.path.join(pathMasses, "{}.masses".format(modelnumber))
        obs = load_obs(Mij_obspath)
        sample.update_observables(obs, "masses")

    os.remove(save_path)
    f = open(save_path, "wb")
    with f: pk.dump(sample, f)

def write_results(nF):

    headers = ('weight', 'like')
    paths   = [ os.path.join(pathSamples, m) for m in os.listdir(pathSamples)]
    labels  = []
    latexs = []

    ncols = 2
    
    latexFile = open(os.path.join(pathLocalData, "latex.localdata"), "rb")
    
    with latexFile: latexDefs = pk.load(latexFile)
    
    fkeys = []
    vkeys = []
    pkeys = []
    for k in sorted(latexDefs.keys()):
        if k.startswith("f"): fkeys.append(k)
        elif k.startswith("v"): vkeys.append(k)
        elif k.startswith("p"): pkeys.append(k)
        else: raise KeyError, k
    
    sortedKeys = fkeys + vkeys + pkeys
    
    for key in sortedKeys:
        headers += (key,)
        labels.append(key)
        ncols += 1
        latexs.append(latexDefs[key])
        
    # Since dictionary items will be read in arbitrary order, resort list based on param number
    

    if len(os.listdir(path2pf)) > 0:
        labs = ["ns", "alpha", "2pf_t"]
        lats = ["n_s", "\\alpha", "T_{2pf}"]

        for item in labs:
            headers+=(item,)
            labels.append(item)
            ncols+=1

        for item in lats: latexs.append(item)


    if len(os.listdir(path3pf)) > 0:
        
        for d in fNLDict:
            labs = [d['name'], "{}_t".format(d['name'])]
            lats = ["f_{NL}"+"^{}".format(d['name']),
                    "T_{}".format("{"+d['name']+"}") ]

            for item in labs:
                headers+=(item,)
                labels.append(item)
                ncols+=1

            for item in lats: latexs.append(item)

    if len(os.listdir(pathMasses)) > 0:
        
        for i in range(nF):
            headers+=("m{}".format(i), )
            labels.append("m{}".format(i))
            latexs.append("m_{}^2/H^2".format(i))
            ncols+=1
            
        headers+=("masses_t",)
        labels.append("masses_t")
        latexs.append("T_{M}")
        ncols+=1


    dicts_evolving  = []
    dicts_adiabatic = []
    failed_samples  = []

    for p in paths:
        f = open(p, "rb")
        with f:
            s = pk.load(f)
            d = s.line_dict()

            if len(d) == ncols:
                if s.adiabatic is True:
                    dicts_adiabatic.append(d)
                else:
                    dicts_evolving.append(d)

            else:
                failed_samples.append(p+"\n")

    f = open(os.path.join(pathRoot, "outputs", "r_adiabatic.txt"), "w")
    g = open(os.path.join(pathRoot, "outputs", "r_evolving.txt" ), "w")

    with f:
        w = csv.DictWriter(f, headers, delimiter=' ')
        for d in dicts_adiabatic:
            w.writerow(d)

    with g:
        v = csv.DictWriter(g, headers, delimiter=' ')
        for d in dicts_evolving:
            v.writerow(d)

    pfile1 = os.path.join(pathRoot, "outputs", "r_adiabatic.paramnames")
    pfile2 = os.path.join(pathRoot, "outputs", "r_evolving.paramnames")
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
