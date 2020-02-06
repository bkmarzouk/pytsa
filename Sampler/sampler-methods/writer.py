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


def bg_summary():
    
    # Build paths to pickled minor stats
    bgd = os.path.join(pathStats, "bg") # bg dir
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
        
    # Write single
    f = open(os.path.join(pathStats, "summary_bg.pk"), "wb")
    
    with f: pk.dump(totDict, f)
    
    
def write_error_report():
    
    # Build paths to pickled minor stats
    twoptDir = os.path.join(pathStats, "2pf")
    threeptDir = os.path.join(pathStats, "3pf")
    
    twopt_paths = [os.path.join(twoptDir, p) for p in os.listdir(twoptDir)]
    threept_paths = [os.path.join(threeptDir, p) for p in os.listdir(threeptDir)]
    
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
    allKeys = [str(f) for f in allFlags]
    
    totDict   = {k: 0 for k in allKeys}

    # Combine 2pt error stats
    for p in twopt_paths + threept_paths:
    
        f = open(p, "rb")
        with f: s = pk.load(f)
    
        k = str(s['flag'])
        totDict[k] += 1

    for p in twopt_paths + threept_paths: os.remove(p)

    f = open(os.path.join(pathStats, "summary_2pf_3pf.pk"), "wb")
    with f: pk.dump(totDict, f)

    g = open(os.path.join(pathStats, "summary_bg.pk"), "rb")
    with g:
        bgStats = pk.load(g)
    
    bgFailCounter = 0
    totFailCounter = 0
    
    n_obtained = len(os.listdir(pathSamples))
    
    for k in allKeys:
        bgFailCounter += bgStats[k]
        totDict[k]    += bgStats[k]
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
        "Unable to find momenta:            {}\n".format(totDict['-31']),
        "Unable to find ICs (2pf):          {}\n".format(totDict['-32']),
        "Unable to find ICs (3pf):          {}\n".format(totDict['-33']),
        "Model violation:                   {}\n\n".format(totDict['-34']),
        
        "Total number of rejections:        {}\n\n".format(totFailCounter),
        
        "---- Probability of inflating (background-level):\n\n",
        
        "Sample efficiency:                 {}\n".format(sampleEfficiency),
        "Total time to build ensemble:      {}\n".format(bgStats['time']),
        "Average time to obtain sample:     {}\n".format(averageSampleTime)
        
    ]

    summary_file = open(os.path.join(pathRoot, "outputs", "sampling_summary.txt"), "w")
    
    with summary_file as f:
        summary_file.writelines(lines)
    

def write_results(nF):

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
            
            print dir(s)
            
            d = s.line_dict()

            if s.reject is False:
                if s.adiabatic is True:
                    dicts_adiabatic.append(d)
                else:
                    dicts_evolving.append(d)

            else:
                print "FAILED SAMPLE: {}".format(p)
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
