# Import system tools
import os
import sys
import pickle as pk
import csv
import shutil

rootpath = os.environ['PyTSamplerRoot']; import config as cfg
path_2pf = os.path.join(rootpath, "2pf")
path_3pf = os.path.join(rootpath, "3pf")
Mij_path = os.path.join(rootpath, "Mij")
cfg_path = os.path.join(rootpath, "config.py")
smp_path = os.path.join(rootpath, "samples")

def load_obs(pk_path):
    f = open(pk_path, "rb")
    with f: obs = pk.load(f)
    return obs

def update_samples(modelnumber):

    # Load sample
    save_path = os.path.join(smp_path, "{}.sample".format(modelnumber))
    f=open(save_path, "rb")
    with f: sample = pk.load(f); assert sample.modelnumber == modelnumber

    # Get information about observables
    observables = cfg.computations

    if observables['2pf'] is True:
        twoPt_path = os.path.join(path_2pf, "{}.2pf".format(modelnumber))
        obs = load_obs(twoPt_path)
        sample.update_observables(obs, "2pf")

    if observables['3pf'] is True:
        which3pf = cfg.which3pf
        for c in which3pf:
            cname = c['config_name']
            cpath = os.path.join(path_3pf, "{m}.{c}".format(m=modelnumber, c=cname))
            obs = load_obs(cpath)
            sample.update_observables(obs, cname)

    if observables['Mij'] is True:
        Mij_obspath = os.path.join(Mij_path, "{}.Mij".format(modelnumber))
        obs = load_obs(Mij_obspath)
        sample.update_observables(obs, "Mij")

    os.remove(save_path)
    f = open(save_path, "wb")
    with f: pk.dump(sample, f)

def write_results(nF):

    headers = ('weight', 'like')
    paths   = [ os.path.join(smp_path, m) for m in os.listdir(smp_path)]
    labels  = []
    latexs = []

    ncols = 2

    # Get information about observables
    observables = cfg.computations

    if observables['2pf'] is True:
        labs = ["ns", "alpha", "2pf_t"]
        lats = ["n_s", "\\alpha", "T_{2pf}"]

        for item in labs:
            headers+=(item,)
            labels.append(item)
            ncols+=1

        for item in lats: latexs.append(item)


    if observables['3pf'] is True:
        which3pf = cfg.which3pf
        for c in which3pf:
            labs = [c['config_name'], "{}_t".format(c['config_name'])]
            lats = ["f_{NL}"+"^{}".format(c['config_name']),
                    "T_{}".format("{"+c['config_name']+"}") ]

            for item in labs:
                headers+=(item,)
                labels.append(item)
                ncols+=1

            for item in lats: latexs.append(item)

    if observables['Mij'] is True:
        for i in range(nF):
            headers+=("m{}".format(i), )
            labels.append("m{}".format(i))
            latexs.append("m_{}^2/H^2".format(i))
            ncols+=1
        headers+=("Mij_t",)
        labels.append("Mij_t")
        latexs.append("T_{M}")
        ncols+=1



    pinfo = cfg.parameter_values
    for item in pinfo:
        if item['LaTeX'] is not None:
            pnumb = item['ParameterNumber']
            pname = "p{}".format(pnumb)
            headers+=(pname,)
            labels.append(pname)
            latexs.append(item['LaTeX'])
            ncols+=1

    dicts_evolving  = []
    dicts_adiabatic = []
    failed_samples  = []

    for p in paths:
        f = open(p, "rb")
        with f:
            s = pk.load(f)
            d = s.line_dict()

            print s.observables

            if len(d) == ncols:
                if s.adiabatic is True:
                    dicts_adiabatic.append(d)
                else:
                    dicts_evolving.append(d)

            else:
                failed_samples.append(p+"\n")

    f = open(os.path.join(rootpath, "outputs", "r_adiabatic.txt"), "w")
    g = open(os.path.join(rootpath, "outputs", "r_evolving.txt" ), "w")

    with f:
        w = csv.DictWriter(f, headers, delimiter=' ')
        for d in dicts_adiabatic:
            w.writerow(d)

    with g:
        v = csv.DictWriter(g, headers, delimiter=' ')
        for d in dicts_evolving:
            v.writerow(d)

    pfile1 = os.path.join(rootpath, "outputs", "r_adiabatic.paramnames")
    pfile2 = os.path.join(rootpath, "outputs", "r_evolving.paramnames")
    h = open(pfile1, "w")
    headers = ('label', 'LaTeX')
    with h:
        x = csv.DictWriter(h, headers, delimiter=' ')
        for row in zip(labels, latexs):
            lab, lat = row
            x.writerow({'label': lab, 'LaTeX': lat})
    shutil.copyfile(pfile1, pfile2)

    i = open(os.path.join(rootpath, "outputs", "failed.txt"), "w")
    with i:
        i.writelines(failed_samples)
