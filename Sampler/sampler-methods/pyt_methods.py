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

rootpath = os.environ['PyTSamplerRoot']
path_2pf = os.path.join(rootpath, "2pf")
path_3pf = os.path.join(rootpath, "3pf")
Mij_path = os.path.join(rootpath, "Mij")
cfg_path = os.path.join(rootpath, "config.py")
smp_path = os.path.join(rootpath, "samples")

# Load configuration file, set PyTransport paths and import relevant modules
import config as cfg
pytpath = cfg.system['pytpath']  # PyTransport installation
saveloc = cfg.system['saveloc']  # Save location for sampler outputs
smppath = cfg.system['smppath']  # Sampler build directory

sys.path.append(pytpath)
import PyTransSetup

PyTransSetup.pathSet()

PyT = importlib.import_module(cfg.sampler['PyTransportModule'])
import PyTransScripts as PyS

# Yield configuration data
subevo = cfg.sampler['SubEvolution']
Nexit = cfg.sampler['ExitTime']
minN = cfg.accept_criteria['MinimumEfolds']
tols = np.asarray(cfg.sampler['tols'])
canonical = cfg.end_inflation['Canonical']
if canonical is False:
    conditions = cfg.end_conditions
else:
    conditions = None

from generator import gen_sample
from realization import realization
import record_stats

"""

Model Setup function

For a given Sample and configuration module, we compute the end of inflation subject to appropriate conditions

"""

def Initialize(modelnumber, rerun_model=False):
    
    if rerun_model is False:
        # Generate new sample with parameter & field priori
        n, fvals, vvals, pvals = gen_sample(modelnumber)
        fvals = np.asarray(fvals)
        vvals = np.asarray(vvals)
        pvals = np.asarray(pvals)
    else:
        # Load an existing set of parameters & model priori
        model_path = os.path.join(saveloc, "samples", "{}.sample".format(modelnumber))
        model_file = open(model_path, "rb")
        with model_file as mf:
            model = pk.load(mf)
            fvals = model.fields
            vvals = model.velocities
            pvals = model.parameters

    # Cross-check the number of fields needed for the potential
    nF = PyT.nF()
    assert nF == len(fvals), "PyTransport module defined nF = {nF}, nF = {lF} != {nF}".format(nF=nF, lF=len(fvals))

    # Compute potential based on the fields and model parameters
    V = PyT.V(fvals, pvals)

    # Iterate over velocity definitions
    if "SlowRoll" in vvals:
        assert all(
            v == "SlowRoll" for v in vvals), "If one field velocity is defined via the slow roll equation, all must be"

        # Compute derivatices of the potential
        dV = PyT.dV(fvals, pvals)

        # Concatenate field values with the slow-roll velocities
        initial = np.concatenate((fvals, -dV / np.sqrt(3 * V)))

    # If specified initial velocities: Read directly from configuration module
    else:
        initial = np.concatenate((fvals, vvals))

    # If we choose to end inflation with eps = 1
    if canonical is True:

        Nend = PyT.findEndOfInflation(initial, pvals, tols, 0.0, 10000)

        # Integration failure
        if type(Nend) == type('S32'):
            return "feoi_fail"
        
        # If no end to inflation is found
        elif Nend is None:
            return "eternal"
        
        # Attempt computation of background
        else:
            back = PyT.backEvolve(np.linspace(0, Nend, 1000), initial, pvals, tols, False)
            
            # If not numpy array, background computation failed
            if type(back) != np.ndarray:
                return "back_fail"

    # Inflation ends due to exotic conditioning: Integrate as far as possible to maximise background data
    else:

        # Compute extended background evolution
        back_ext = PyS.ExtendedBackEvolve(initial, pvals, PyT)

        # If integration / background failure / eternal model is reported, return this information
        if back_ext ["feoi_fail", "back_fail", "eternal"]:
            return back_ext

        # Otherwise unpack background and efolding where epsilon = 1
        back, Nepsilon = back_ext

        # Iterate over conditions list
        for item in conditions:

            # Unpack conditions for fields values subject to background data
            field_number = item['FieldNumber']
            field_value = item['FieldValue']
            from_above = item['FromAbove']
            successful = item['Successful']

            # Transpose evolution to extract independent evolution for field in question
            Nevo = back.T[0]
            Fevo = back.T[1 + field_number]
            Vevo = back.T[1 + field_number + PyT.nF()] # TODO: Consider velocity conditions on fields

            # Try to find instances of field conditions within background evolution
            try:
                # If conditional field value is approached from above
                if from_above is True:
                    Nloc = np.where(Fevo <= field_value)[0][0] # Find first instance of condition being True
                    
                # If conditional field value is approached from below
                elif from_above is False:
                    Nloc = np.where(Fevo >= field_value)[0][0] # Fine first instance of condition being True
                
                else:
                    raise ValueError, "Unrecognized value assignment 'from above': {}".format(from_above)

                # If condition successfully ends inflation, redefine Nend to this point
                if successful is True:
                    Nend = Nevo[Nloc]
                    
                # If condition unsuccessfully terminates inflation, redefine Nend as "Violated" model
                elif successful is False:
                    Nend = "Violated"
                    
                else:
                    raise ValueError, "Unrecognized value assignment: {}".format(successful)

            # Index error is raised if field value condition cannot be foudn
            except IndexError:

                # If the condition corresponded to a successful end to the evolution, return "NotFound"
                if successful is True:
                    Nend = "NotFound"
                
                # If the condition was supposed to terminate the evolution, pass (everything is OK)
                elif successful is False:
                    pass
                
                else:
                    raise ValueError, "Unrecognized value assignment: {}".format(successful)

    # We infer whether the background evolution was successful and whether it was too short subject to definition
    if Nend in ["Violated", "NotFound"] or Nend < minN:
        if Nend < minN:
            return "Short"
        return Nend

    # If inflation lasts for more than 70 efolds, reposition ICS s.t. we avoid exponentially large momenta
    if Nend > 70.0:
        
        # Search background from bottom to top
        for row in np.flipud(back):

            # When the background step is 70 efolds away from the "end" of inflation
            if Nend - row[0] > 70:
                
                # Redefine initial conditions to this point & rerun background compution
                initial = row[1:]
                Nend = Nend - row[0]
                
                # For canonical end, this should be straightforward
                if canonical is True:
                    t = np.linspace(0.0, Nend, int((Nend - 0) * 3))
                    back_adj = PyT.backEvolve(t, initial, pvals, tols, False)
                
                # For non-canonical
                else:
                    
                    # *Small* chance that extension may mess up, returning ValueError
                    try:
                        back_adj, Neps = PyS.ExtendedBackEvolve(initial, pvals, PyT)
                    except ValueError:
                        back_adj = None
                
                # Break out of search window
                break
    else:
        back_adj = back

    # If rescaling procedure fails: Revert to initial background compution
    if type(back_adj) == np.ndarray:
        pass
    else:
        back = back_adj

    # Attempt computation of momenta at horizon crossing, this will be used in all calculations of observables
    kExit = PyS.kexitN(Nend - Nexit, back, pvals, PyT)
    
    # Asses success by data type of momenta result
    if type(kExit) not in [float, np.float, np.float32, np.float64]:
        return "kexitn_fail"

    # We now test for an adiabatic limit: We begin by assuming this is True, then test for violating conditions
    adiabatic = True

    if cfg.accept_criteria['TestAdiabaticLimit'] is True:
        
        print "-- Searching for adiabatic limit: %05d" % modelnumber
        
        # Define number of efolds from end of inflation which should be adiabatic TODO: Wrap this into config. file
        Nadiabatic = 1.
        
        # Find first efoldign to perform eigenvalue test
        Nadiabatic_start = np.where(back.T[0] >= Nend - Nadiabatic)[0][0]
        
        # Compute mass-matrix evolution from this point
        Mij_end = PyS.MijEvolve(back[last5:], pvals, PyT)

        # For each mass-matrix evolution step
        for item in Mij_end:
            
            # Determine if lightest mass is tachyonic
            tachyon = item[1] < 0
            
            # Determine number of hubble scale heavier masses
            rest_large = sum([eig >= 1 for eig in item[2:]]) == nF - 1

            # If there is no tachyon or the remaining masses aren't large for the portion of evolution
            if not (tachyon and rest_large):
                
                # Change adiabatic property to False, and break out of search window
                adiabatic = False
                
                print "-- Adiabatic limit not found: %05d" % modelnumber
                break

    # Record all sample data
    realization(modelnumber, fvals, vvals, pvals,
                back, adiabatic, Nend, kExit, smp_path, rerun_model)

    print "-- Generated sample: %05d" % modelnumber

    # Return model number, acts as flag to confirm successful sample
    return modelnumber


def DemandSample(modelnumber, sample_stats_dir):
    """ Repeat initialization until successful sample is found """
    ii = -1
    
    # If successful sample generation, the model number is returned
    while ii != modelnumber:
        ii = Initialize(modelnumber)
        
        # If str type, ii is a key for stats records
        if type(ii) == str:
            key = ii
            
        # If model number, sum number of iterations
        elif ii == modelnumber:
            key = "end"
        else:
            raise KeyError, ii
        
        # Record statistics on background efforts
        record_stats.log_stats(modelnumber, key, sample_stats_dir)
        

def DemandSample_rerun(modelnumber):
    """ Repeat initialization until successful sample is found """
    ii = -1
    n_trials = 1
    n_kexitn_fail = 0
    n_feoi_fail = 0
    
    while ii != modelnumber:
        ii = Initialize(modelnumber, rerun_model=True)
        n_trials += 1
        if ii == "kexitn_fail":
            n_kexitn_fail += 1
            ii = -1
        elif ii = "feoi_fail":
            n_feoi_fail += 1
        else:
            pass


def Mij(modelnumber):
    # Begin by loading pickled sample
    sampler_path = os.path.join(smp_path, "{}.sample".format(modelnumber))
    assert os.path.exists(sampler_path), "Unable to locate sample location: {}".format(sampler_path)
    s = open(sampler_path, "rb")
    with s:
        sample = pk.load(s)

    # Cross check model numbers to ensure that we are reading and computing with the correct data
    assert sample.modelnumber == modelnumber, "Pickled model data does not match! {m1} != {m2}".format(
        m1=modelnumber, m2=sample.modelnumber
    )

    # Unload core data for observable calculation
    back = sample.background
    backT = back.T
    Nend = sample.Nend
    assert Nend is not None, "Nend = None should not be a saved sample"
    pvals = np.asarray(sample.parameters)

    # Define locations and savenames for spectral index and running
    Mij_savepath = os.path.join(Mij_path, "{}.Mij".format(modelnumber))

    nF = PyT.nF()

    # Begin timer
    ti = time.clock()

    Nstar = Nend - cfg.sampler['ExitTime'];
    Nevo = backT[0]

    if Nstar in Nevo:
        Nstar_idx = [Nevo.index(Nstar)]
        back_star = back[Nstar_idx]

    else:
        # We consider the evolution of the mass matrix at +/- 1 eFold about the pivot scale
        idx_above = np.where(np.logical_and(Nevo>Nstar, Nevo<Nstar+1.))[0]
        idx_below = np.where(np.logical_and(Nevo < Nstar, Nevo > Nstar-1.))[0]
        Nstar_idx = list(np.concatenate([idx_below, idx_above]))

        # Define efold range to "smooth" over
        Nsmooth = [Nevo[idx] for idx in Nstar_idx]

        # Construct splines for the field evolution w.r.t. N time
        fsplines = [
            UnivariateSpline(
                Nsmooth, [backT[i][idx] for idx in Nstar_idx]) for i in range(1, 1 + nF)
        ]
        
        # Construct splines for the velocity evolution w.r.t. N time
        vsplines = [
            UnivariateSpline(
                Nsmooth, [backT[i][idx] for idx in Nstar_idx]) for i in range(1 + nF, 1 + 2 * nF)
        ]

        # Evaluate the background at N_star = Nend - Nexit
        back_star = np.concatenate([np.array([Nstar]), np.array([s(Nstar) for s in fsplines + vsplines])])
        
    # Compute mass matrix eigenvalues
    Mij = PyS.MijEvolve(np.array([back_star]), pvals, PyT, scale_eigs=True)[0][1:]

    # Write dictionary values for masses
    Mij_eig = {}
    for i in range(len(Mij)):
        Mij_eig['m{}'.format(i)] = Mij[i]

    # End timer
    Mij_eig['time'] = time.clock() - ti

    Mij_file = open(Mij_savepath, "wb")
    with Mij_file:
        pk.dump(Mij_eig, Mij_file)
        

def SpectralIndex(modelnumber):
    # Begin by loading pickled sample
    sampler_path = os.path.join(smp_path, "{}.sample".format(modelnumber))
    assert os.path.exists(sampler_path), "Unable to locate sample location: {}".format(sampler_path)
    s = open(sampler_path, "rb")
    with s:
        sample = pk.load(s)

    # Cross check model numbers to ensure that we are reading and computing with the correct data
    assert sample.modelnumber == modelnumber, "Pickled model data does not match! {m1} != {m2}".format(
        m1=modelnumber, m2=sample.modelnumber
    )

    # Unload core data for observable calculation
    back = sample.background
    Nend = sample.Nend;
    assert Nend is not None, "Nend = None should not be a saved sample"
    pvals = np.asarray(sample.parameters)
    
    kExit = sample.kExit

    # Define savepath for spectral index and running
    twoPt_savepath = os.path.join(path_2pf, "{}.2pf".format(modelnumber))

    # Begin timing calculation
    ti = time.clock()

    try:
        
        # Define range of N times about exit scale of interest
        Nrange = [Nend - Nexit - 2.0 + float(i) for i in [0, 1, 3, 4]]

        # Compute modes that exit at corresponding times
        kExitrange = [PyS.kexitN(N, back, pvals, PyT) for N in Nrange]
        
        Nrange.insert(2, Nend - Nexit)
        kExitrange.insert(2, kExit)

        print "Nexit", Nrange
        print "kExit", kExitrange

        # Check that momenta size is monotonically increasing
        assert all(kExitrange[i]<kExitrange[i+1] for i in range(len(kExitrange)-1))

        # Get initial conditions for momenta
        ICsEvos = [PyS.ICs(subevo, kE, back, pvals, PyT) for kE in kExitrange]

        # Check at initial conditions array(s) are correct length of correct dtype
        for item in ICsEvos:
            assert len(item) == 2, "Initial conditions incorrect length"
        for item in ICsEvos:
            assert type(item[1]) == np.ndarray, "Initial conditions are wrong dtype"

        # Define efoldings at which we want to evaluate spectral index at
        Nevals = [np.array([Nstart[0], Nend]) for Nstart in ICsEvos]

        # Compute the two point function for each mode at evaluation times
        twoPt = [PyT.sigEvolve(NkBack[0], NkBack[1], NkBack[2][1], pvals, tols, False).T for NkBack in zip(
            Nevals, kExitrange, ICsEvos
        )]

        # Log power spectra and Fourier modes
        logPz = [np.log(xx[1][-1]) for xx in twoPt]
        logkE = [np.log(kE) for kE in kExitrange]

        # Compute spectral index and running for pivot scale at end of inflation
        ns = UnivariateSpline(
            logkE, logPz, k=4, s=1e-15).derivative()(np.log(kExitrange[2])) + 4.
        alpha = UnivariateSpline(
            logkE, logPz, k=4, s=1e-15).derivative(2)(np.log(kExitrange[2])) + 0.
        
        print ns, alpha
        


    except (AssertionError, AttributeError) as e:
        
        print "-- ns Failure"

        ns = None
        alpha = None

    # End calculation timer
    twoPt_dict = {'ns': ns, 'alpha': alpha, 'time': time.clock() - ti}

    twoPt_file = open(twoPt_savepath, "wb")

    with twoPt_file:
        pk.dump(twoPt_dict, twoPt_file)


# fNL calculation
def fNL(modelnumber, which3pf):
    # print "-- Begin fNL calculation: Model {m}".format(m = modelnumber)

    # Begin by loading pickled sample
    sampler_path = os.path.join(smp_path, "{}.sample".format(modelnumber))
    assert os.path.exists(sampler_path), "Unable to locate sample location: {}".format(sampler_path)
    s = open(sampler_path, "rb")
    with s:
        sample = pk.load(s)

    # Cross check model numbers to ensure that we are reading and computing with the correct data
    assert sample.modelnumber == modelnumber, "Pickled model data does not match! {m1} != {m2}".format(
        m1=modelnumber, m2=sample.modelnumber
    )

    # Unload core data for observable calculation
    back = sample.background
    Nend = sample.Nend
    assert Nend is not None, "Nend = None should not be a saved sample"
    kExit = sample.kExit
    pvals = np.asarray(sample.parameters)

    # Unpack dictionary data
    name = which3pf['config_name']
    alpha = which3pf['alpha']
    beta = which3pf['beta']

    # Define location and savepath for fNL computation
    fNL_savepath = os.path.join(path_3pf, "{num}.{ext}".format(num=modelnumber, ext=name))

    # Begin timing calculation
    ti = time.clock()

    # Apply parameterization to get triangle momenta
    k1 = kExit / 2. - beta * kExit / 2.;
    k2 = kExit * (1. + alpha + beta) / 4.;
    k3 = kExit * (1. - alpha + beta) / 4.

    # Find ICs for largest scale as this will cross the Horizon first
    kmin = np.min([k1, k2, k3])

    try:

        ICs = PyS.ICs(subevo, kmin, back, pvals, PyT)
        assert len(ICs) == 2, "Initial conditions incorrect length: model {m}, config {c}".format(
            m=modelnumber, c=name
        )

        assert type(ICs[1]) == np.ndarray, "Initial conditions incorrect dtype: model {m}, config {c}".format(
            m=modelnumber, c=name
        )

        # Based on our configuration, we then calulate the 3pf
        threePt = PyT.alphaEvolve(np.array([ICs[0], Nend]), k1, k2, k3, ICs[1], pvals, tols, True).T

        # Get Bispectrum of zeta, and power spectra for modes at end of inflation
        Pz1, Pz2, Pz3, Bz = [threePt[i][-1] for i in range(1, 5)]

        # Compute array of fNL values
        fNL = (5. / 6.) * Bz / (Pz1 * Pz2 + Pz2 * Pz3 + Pz1 * Pz3)

    # Except TypeError, typically propagates from bad BackExitMinus, i.e. float rather than ndarray
    except (AssertionError,  AttributeError) as e:
        fNL = None

    fNL_dict = {'{}'.format(name): fNL, 'time': time.clock() - ti}

    fNL_file = open(fNL_savepath, "wb")

    with fNL_file:
        pk.dump(fNL_dict, fNL_file)
        

def computations(mn_calc):
    modelnumber, calculation = mn_calc
    which3pf = cfg.which3pf

    with warnings.catch_warnings(record=True) as w:

        if calculation == "Mij":
            print "-- Start Mij: model %06d" % modelnumber
            Mij(modelnumber)
            print "--   End Mij: model %06d" % modelnumber

        if calculation == "2pf":
            print "-- Start 2pf: model %06d" % modelnumber
            SpectralIndex(modelnumber)
            print "--   End 2pf: model %06d" % modelnumber

        for config in which3pf:
            if config['config_name'] == calculation:
                print "-- Start fNL: model %06d, config {}".format(config['config_name']) % modelnumber
                fNL(modelnumber, config)
                print "--   End fNL: model %06d, config {}".format(config['config_name']) % modelnumber

        if len(w) > 0:
            print "-- {t} TASK FAILURE, MODEL {m}".format(t=calculation, m=modelnumber)