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

# Retrieve paths from environment variables
envKeys = [
    'PyTS_pathPyT', 'PyTS_pathRoot', 'PyTS_path2pf', 'PYTS_path3pf', 'PyTS_pathMasses', 'PyTS_pathSamples',
    'PyTS_pathStats', 'PyTS_pathClasses', 'PyTS_pathMethods', 'PyTS_pathLocalData'
]

pathPyT, pathRoot, path2pf, path3pf, pathMasses, \
pathSamples, pathStats, pathClasses, pathMethods, pathLocalData = [os.environ[eK] for eK in envKeys]


# get transport data dictionary
transportPath = os.path.join(pathLocalData, "transport.localdata")
transportFile = open(transportPath, "rb")
with transportFile: transportDict = pk.load(transportFile)


# configure internal pytransport paths
sys.path.append(pathPyT)
import PyTransSetup
PyTransSetup.pathSet()


# Import module and get other required data for computing observables
PyT = importlib.import_module(transportDict['transportModule'])
import PyTransScripts as PyS

transportKeys = ["exitN", "subN", "intTols", "adiabaticN", "minN"]
Nexit, subevo, tols, adiabaticN, minN = [transportDict[tK] for tK in transportKeys]


# Get fNL data
fNLPath = os.path.join(pathLocalData, "fNL.localdata")
fNLFile = open(fNLPath, "rb")
with fNLFile: fNLDict = pk.load(fNLFile)


# Import sample generator & other useful class structures
from generator import genSample
from realization import realization as new_sample
import record_stats



def Initialize(modelnumber, rerun_model=False):
    
    tmax_bg = int(os.environ['tmax_bg'])
    
    if rerun_model is False:
        # Generate new sample with parameter & field priori
        n, fvals, vvals, pvals = genSample(modelnumber)

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

    # If SR is in vvals, compute initial velocity value via slowroll equation
    if "SR" in vvals:
        
        # Construct np float array
        vvals_ = np.zeros(np.shape(vvals), dtype=float)
        
        # Compute derivatices of the potential
        dV = PyT.dV(fvals, pvals)

        # Compute velocities via slowroll equation
        vvalsSR = -dV / np.sqrt(3 * V)

        # For each velocity def.
        for jj in range(len(vvals)):
            
            #  If def. is SR, assign slow-roll velocity to float array
            if vvals[jj] == "SR":
                vvals_[jj] = vvalsSR[jj]
            else:
                # Otherwise force float-type to other array elements
                vvals_[jj] = float(vvals[jj])
                
        # Redefine vvals
        vvals = vvals_

    # Concatenate field and velocity values into initial conditions array
    initial = np.concatenate((fvals, vvals))

    """ We have the following flags in place:
    
    Timeout flags:
    -10 : findEndOfInflation
    -11 : backEvolve
    -12 : sigEvolve
    -13 : alphaEvolve
    
    Integration error flags:
    -20 : findEndOfInflation
    -21 : beckEvolve
    -22 : sigEvolve
    -23 : alphaEvolve
    
    Misc. PyT:
    -30 : Nend < Nmin
    -31 : k not found
    -32 : No ICsBM (2pf)
    -33 : No ICsBM (3pf)
    -34 : Model violation
    -35 : Eternal inflation
    
    """
    
    acceptIfAnyPath = os.path.join(pathLocalData, "acceptIfAny.localdata")
    rejectIfAnyPath = os.path.join(pathLocalData, "rejectIfAny.localdata")
    acceptIfAllPath = os.path.join(pathLocalData, "acceptIfAll.localdata")
    rejectIfAllPath = os.path.join(pathLocalData, "rejectIfAll.localdata")
    
    # Define canonical end to inflation if there are no field conditions to search through
    canonical = not (os.path.exists(acceptIfAnyPath) or os.path.exists(rejectIfAnyPath)
                     or os.path.exists(acceptIfAllPath) or os.path.exists(rejectIfAllPath))

    # If we choose to end inflation with eps = 1
    if canonical is True:
        
        breakFlag = True
        
        # Compute end of inflation
        Nend = PyT.findEndOfInflation(initial, pvals, tols, 0.0, 10000, tmax_bg, True)

        # Integration fails / eternal
        if type(Nend) is tuple: return Nend
        
        # Attempt computation of background
        back = PyT.backEvolve(np.linspace(0, Nend, 1000), initial, pvals, tols, False, tmax_bg, True)
        
        # If not numpy array, background computation failed
        if type(back) is tuple: return back

        # Get number of rows in bg evo
        rowCount = len(back)

    else:
        # Load model constriants
        
        if os.path.exists(acceptIfAnyPath):
            f = open(acceptIfAnyPath, "rb")
            with f:
                acceptIfAny = pk.load(f)
                acceptIfAnyIdx = acceptIfAny['idx']
                acceptIfAnyMin = acceptIfAny['min']
                acceptIfAnyMax = acceptIfAny['max']
                
        else: acceptIfAny = None
        
        if os.path.exists(rejectIfAnyPath):
            f = open(rejectIfAnyPath, "rb")
            with f:
                rejectIfAny = pk.load(f)
                rejectIfAnyIdx = rejectIfAny['idx']
                rejectIfAnyMin = rejectIfAny['min']
                rejectIfAnyMax = rejectIfAny['max']
                
        else: rejectIfAny = None
        
        if os.path.exists(acceptIfAllPath):
            f = open(acceptIfAllPath, "rb")
            with f:
                acceptIfAll = pk.load(f)
                
        else: acceptIfAll = None

        if os.path.exists(rejectIfAllPath):
            f = open(rejectIfAllPath, "rb")
            with f:
                rejectIfAll = pk.load(f)
                
        else: rejectIfAll = None
        
        # Start background row counter
        rowCount = 0
        
        # We flag to break out of loop cycle if there is a condition that ends inflation
        breakFlag = False if (acceptIfAny is not None or acceptIfAll is not None) else True
        
        # If we are looking for a specific exit condition, then we compute extended BG evo
        if breakFlag is False:
            
            # Integrate over extra-large window until integrator falls over / timeout
            back = PyT.backEvolve(np.linspace(0, 2 * 1e4, 1e4), initial, pvals, tols, False, tmax_bg, True)

            # Try and extend by a few more iterations with reduced-tols
            extTols = tols *1e-1
            tbgi = time.clock()
            n_bg_iter = 2
            for n in range(n_bg_iter):
                EXT = PyT.backEvolve(
                    np.linspace(back[-1][0], back[-1][0] + 2 * 1e4, 1e4),
                    back[-1][1:], pvals, extTols, False, tmax_bg, True)
                if len(EXT) > 1:
                    back = np.vstack((back, EXT[1:]))
                else: break
            tbgf = time.clock()
            
            # If too short and time taken exceeds max timers x n_iter, timeout flag
            if back[-1][0] < minN and tbgf-tbgi > n_bg_iter*(tmax_bg)-0.5:
                # timeout
                return (-11, back[-1][0])
            
            # If too short and time taken below max timer x n_iter, integrator flag
            elif back[-1][0] < minN and tbgf-tbgi < n_bg_iter*(tmax_bg)-0.5:
                # int
                return (-21, back[-1][0])
            
            # Simply too short
            elif back[-1][0] < minN:
                # minN
                return (-30, back[-1][0])
            
            else:
                pass
                
                

            # Simply change this to return if int? All flags *should" be handled
            if type(back) is tuple: return back
            
            
        else:
            
            # Integrate over extra-large window until eps > 1
            back = PyT.backEvolve(np.linspace(0, 2*1e4, 1e6), initial, pvals, tols, True, tmax_bg, True)
            
            if type(back) is tuple: return back
            
            # If inflation is too short, may as well return flag now
            if back[-1][0] < minN: return
        
        for row in back:
            
            """ Check for violations """

            if rejectIfAny is not None:
                if np.any(np.logical_and(row[rejectIfAnyIdx] > rejectIfAnyMin, row[rejectIfAnyIdx] < rejectIfAnyMax)):
                    # bad
                    return -34
            
            if rejectIfAll is not None:
                for d in rejectIfAll:
                    idx = rejectIfAll[d]['idx']
                    min = rejectIfAll[d]['min']
                    max = rejectIfAll[d]['max']
                    if np.all(np.logical_and(row[idx] > min, row[idx] < max)):
                        # bad
                        return -34
                
            if acceptIfAny is not None:
                if np.any(np.logical_and(row[acceptIfAnyIdx] > acceptIfAnyMin, row[acceptIfAnyIdx] < acceptIfAnyMax)):
                    # good
                    breakFlag = True
                    break
                
            if acceptIfAll is not None:
                for d in acceptIfAll:
                    idx = acceptIfAll[d]['idx']
                    min = acceptIfAll[d]['min']
                    max = acceptIfAll[d]['max']
                    if np.all(np.logical_and(row[idx] > min, row[idx] < max)):
                        # good
                        breakFlag = True
                        break
                
                if breakFlag: break
                    
            Nend = row[0]
            rowCount += 1
    
    if Nend < minN: return -30

    if breakFlag is not True: return -34
    
    back = back[:rowCount]
    
    # If inflation lasts for more than 70 efolds, reposition ICS s.t. we avoid exponentially large momenta
    if Nend > 70.0:
        
        # Search background from bottom to top
        for row in np.flipud(back):

            # When the background step is 70 efolds away from the "end" of inflation
            if Nend - row[0] > 70:
                
                # Redefine initial conditions to this point & rerun background compution
                initial = row[1:]
                Nend = Nend - row[0]
                
                t = np.linspace(0.0, Nend, int((Nend - 0) * 3))
                back_adj = PyT.backEvolve(t, initial, pvals, tols, False, tmax_bg, True)
                
                if type(back_adj) != np.ndarray:
                    back_adj = back
                
                if back_adj[-1][0] != back[-1][0]:
                    back_adj = back
                
                # Break out of search window
                break
    else:
        back_adj = back

    # Attempt computation of momenta at horizon crossing, this will be used in all calculations of observables
    kExit = PyS.kexitN(Nend - Nexit, back_adj, pvals, PyT)
    
    # Asses success by data type of momenta result
    if type(kExit) not in [float, np.float, np.float32, np.float64] or np.isnan(kExit) or np.isinf(kExit):
        return -31

    # We now test for an adiabatic limit: We begin by assuming this is True, then test for violating conditions
    adiabatic = True
    
    # Define number of efolds from end of inflation which should be adiabatic TODO: Wrap this into config. file
    
    # Find first efoldign to perform eigenvalue test
    try:
        adiabaticN_start = np.where(back.T[0] >= Nend - adiabaticN)[0][0]
    except IndexError:
        return -21
    
    # Compute mass-matrix evolution from this point
    Mij_end = PyS.evolveMasses(back[adiabaticN_start:], pvals, PyT, scale_eigs=False, hess_approx=False, covariant=False)

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
            
            break

    # Record all sample data
    new_sample(modelnumber, fvals, vvals, pvals, back_adj, adiabatic,
               Nend, kExit, pathSamples, rerun_model)

    print "-- Generated sample: %05d" % modelnumber
    
    assert type(modelnumber) == int, modelnumber

    # Return model number, acts as flag to confirm successful sample
    return modelnumber


def DemandSample(modelnumber):
    """ Repeat initialization until successful sample is found """

    # Get directory for sample stats log
    sample_path = os.path.join(pathStats, "bg", "{}.bg".format(modelnumber))

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
    flagDict = {str(f): 0 for f in allFlags}

    # Start timer
    tstart = time.clock()

    ii = -1

    # If successful sample generation, the model number is returned
    while ii != modelnumber:
        
    
        ii = Initialize(modelnumber)

        print flagDict, ii
    
        # If fail flag received: log statistic
        if type(ii) is tuple:
            flagKey = str(ii[0])
            flagDict[flagKey] += 1
            
        elif type(ii) is int and ii in allFlags:
            flagKey = str(ii)
            flagDict[flagKey] += 1
    
        # If model number, sum number of iterations
        elif ii == modelnumber:
            tfinal = time.clock() - tstart
        
            flagDict['time'] = tfinal
        
            f = open(sample_path, "wb")
        
            with f:
                pk.dump(flagDict, f)
    
        else:
            raise KeyError, ii
        

def DemandSample_rerun(modelnumber):
    """ Repeat initialization until successful sample is found : using rerun samples"""

    # Get directory for sample stats log
    sample_path = os.path.join(pathStats, "bg", "{}.bg".format(modelnumber))

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
    flagDict = {str(f): 0 for f in allFlags}

    # Start timer
    tstart = time.clock()

    ii = -1

    # If successful sample generation, the model number is returned
    while ii != modelnumber:
    
        ii = Initialize(modelnumber, rerun_model=True)
    
        # If fail flag received: log statistic
        if type(ii) is tuple:
            flagKey = str(ii[0])
            flagDict[flagKey] += 1
    
        elif type(ii) is int and ii in allFlags:
            flagKey = str(ii)
            flagDict[flagKey] += 1
    
        # If model number, sum number of iterations
        elif ii == modelnumber:
            tfinal = time.clock() - tstart
        
            flagDict['time'] = tfinal
        
            f = open(sample_path, "wb")
        
            with f:
                pk.dump(flagDict, f)
    
        else:
            raise KeyError, ii


def Mij(modelnumber):
    
    # Unload sample data
    path = os.path.join(pathSamples, "{}.sample".format(modelnumber))
    assert os.path.exists(path), "Unable to locate sample location: {}".format(path)
    s = open(path, "rb")
    with s:
        sample = pk.load(s)

    # Cross-check model number identifier
    assert sample.modelnumber == modelnumber, "Pickled model data does not match! {m1} != {m2}".format(
        m1=modelnumber, m2=sample.modelnumber
    )

    # Unpack background data
    back = sample.background
    Nend = sample.Nend
    assert Nend is not None, "Nend = None -> Sample construction failure!"
    pvals = np.asarray(sample.parameters)

    # Transport background for col evolution rather than row evolution
    backT = back.T
    
    # Build save path for mass-matrix eigenvalues
    Mij_savepath = os.path.join(pathMasses, "{}.masses".format(modelnumber))

    # Get number of fields
    nF = PyT.nF()

    # Begin timer
    ti = time.clock()

    # Define Nstar, i.e. exit time of interest
    Nstar = Nend - Nexit
    
    # Unpack efold evolution
    Nevo = backT[0]

    # If Nstar falls exactly in the background, simply load this background step
    if Nstar in Nevo:
        Nstar_idx = [Nevo.index(Nstar)]
        back_star = back[Nstar_idx]

    # Otherwise we build a spline about Nstar
    else:
        
        # Get background index values +/- 1 efold around the prescribed exit time
        Nfit = 0
        dN = 0.5
        while Nfit < 3:
            dN += 0.5
            idx_above = np.where(np.logical_and(Nevo>Nstar, Nevo<Nstar+dN))[0]
            idx_below = np.where(np.logical_and(Nevo < Nstar, Nevo > Nstar-dN))[0]
            
            Nstar_idx = list(np.concatenate([idx_below, idx_above]))
            Nfit = len(Nstar_idx)

        # Define efold range to "smooth" over
        Nsmooth = [Nevo[idx] for idx in Nstar_idx]
        
        # We check for repeated efold numbers as this will throw the interpolation scheme
        repLog = []
        for ii in range(len(Nsmooth)-1):
            if Nsmooth[ii] == Nsmooth[ii+1]:
                repLog.append(ii+1)
        for index in sorted(repLog, reverse=True):
            del Nstar_idx[index]
            
        
        Nsmooth = [Nevo[idx] for idx in Nstar_idx]
        
        
        # BY default, we choose a cubic spline, but if this is not possible, we reduce the order
        if len(Nstar_idx) <= 3:
            k = len(Nstar_idx) - 1
        else:
            k = 3


        # Construct splines for the field evolution over Nsmooth
        fsplines = [
            UnivariateSpline(
                Nsmooth, [backT[i][idx] for idx in Nstar_idx], k=k) for i in range(1, 1 + nF)
        ]
        
        # Construct splines for the velocity evolution over Nsmooth
        vsplines = [
            UnivariateSpline(
                Nsmooth, [backT[i][idx] for idx in Nstar_idx], k=k) for i in range(1 + nF, 1 + 2 * nF)
        ]

        # Evaluate the background at N_star = Nend - Nexit
        back_star = np.concatenate([np.array([Nstar]), np.array([s(Nstar) for s in fsplines + vsplines])])
        
    # Compute mass matrix eigenvalues at the exit time
    Mij = PyS.evolveMasses(
        np.array([back_star]), pvals, PyT, scale_eigs=False, hess_approx=False, covariant=False)[0][1:]

    # Build dictionary of masses
    Mij_eig = {}
    for i in range(len(Mij)):
        Mij_eig['m{}'.format(i)] = Mij[i]

    # End timer
    Mij_eig['time'] = time.clock() - ti

    # Write binary file
    Mij_file = open(Mij_savepath, "wb")
    with Mij_file:
        pk.dump(Mij_eig, Mij_file)
        
    return 0
        

def SpectralIndex(modelnumber):
    
    # Get timeout
    tmax_2pf = int(os.environ['tmax_2pf'])
    
    # Unload sample data
    path = os.path.join(pathSamples, "{}.sample".format(modelnumber))
    assert os.path.exists(path), "Unable to locate sample location: {}".format(path)
    s = open(path, "rb")
    with s:
        sample = pk.load(s)

    # Cross-check model number identifier
    assert sample.modelnumber == modelnumber, "Pickled model data does not match! {m1} != {m2}".format(
        m1=modelnumber, m2=sample.modelnumber
    )

    # Unpack background information
    back = sample.background
    Nend = sample.Nend
    
    assert Nend is not None, "Nend = None should not be a saved sample"
    pvals = np.asarray(sample.parameters)
    
    # Get exit time momenta
    kExit = sample.kExit

    # Build path for two-point function result
    twoPt_savepath = os.path.join(path2pf, "{}.2pf".format(modelnumber))

    # Begin timer
    ti = time.clock()

    # Compute two-point function

    # Define an efolding range 2-efolds above and below the exit time
    Nrange = [Nend - Nexit - 2.0 + float(i) for i in [0, 1, 2, 3, 4]]

    # Compute the corresponding momenta at horizon exit times
    kExitrange = [PyS.kexitN(N, back, pvals, PyT) for N in Nrange]
    
    for k in kExitrange:
        if np.isnan(k) or np.isinf(k):
            # Build null result file, dump and return error data
            twoPt_dict = {'ns': None, 'alpha': None, 'time': None}
    
            # Write binary file
            twoPt_file = open(twoPt_savepath, "wb")
            with twoPt_file:
                pk.dump(twoPt_dict, twoPt_file)
    
            return {"mn": modelnumber, "flag": -31, "ext": "2pf"}

    # Check momenta are monotonically increasing
    assert all(kExitrange[i]<kExitrange[i+1] for i in range(len(kExitrange)-1)), kExitrange

    # Get initial condition for each momenta
    ICsEvos = [PyS.ICsBM(subevo, kE, back, pvals, PyT) for kE in kExitrange]

    # Check initial conditions array(s)
    for item in ICsEvos:
        if item == (np.nan, np.nan):
            
            # Build null result file, dump and return error data
            twoPt_dict = {'ns': None, 'alpha': None, 'time': None}
    
            # Write binary file
            twoPt_file = open(twoPt_savepath, "wb")
            with twoPt_file:
                pk.dump(twoPt_dict, twoPt_file)
                
            return {"mn": modelnumber, "flag": -32, "ext": "2pf"}

    # Prescribe Nstart and Nend as evaluation times for 2pf run
    Nevals = [np.array([Nstart[0], Nend]) for Nstart in ICsEvos]
    
    # Call sigEvolve to compute power spectra
    twoPt = [PyT.sigEvolve(NkBack[0], NkBack[1], NkBack[2][1], pvals, tols, False, tmax_2pf, True) for NkBack in zip(
        Nevals, kExitrange, ICsEvos
    )]
    
    # If flag return for sigEvolve, return flag data
    for item in twoPt:
        if type(item) is tuple:
            
            # Build null result file, dump and return error data
            twoPt_dict = {'ns': None, 'alpha': None, 'time': None}
    
            # Write binary file
            twoPt_file = open(twoPt_savepath, "wb")
            with twoPt_file:
                pk.dump(twoPt_dict, twoPt_file)
                
            return {"mn": modelnumber, "flag": item[0], "ext": "2pf"}
        
    # Otherwise tranpose data
    twoPt = [item.T for item in twoPt]

    # Log power spectra and momenta
    logPz = [np.log(xx[1][-1]) for xx in twoPt]
    logkE = [np.log(kE) for kE in kExitrange]

    # Compute spectral index and running via spline derivatives of logged values
    ns = UnivariateSpline(
        logkE, logPz, k=4, s=1e-15).derivative()(np.log(kExitrange[2])) + 4.
    
    alpha = UnivariateSpline(
        logkE, logPz, k=4, s=1e-15).derivative(2)(np.log(kExitrange[2])) + 0.

    # End calculation timer and build results dictionary
    twoPt_dict = {'ns': ns, 'alpha': alpha, 'time': time.clock() - ti}

    # Write binary file
    twoPt_file = open(twoPt_savepath, "wb")
    with twoPt_file:
        pk.dump(twoPt_dict, twoPt_file)

    return 0

def fNL(modelnumber, configName):
    
    tmax_3pf = int(os.environ['tmax_3pf'])
    
    # Gather configuration data from configName
    for d in fNLDict:
        if d['name'] == configName:
            name = configName
            alpha = d['alpha']
            beta  = d['beta']
            break

    # Unload sample data
    path = os.path.join(pathSamples, "{}.sample".format(modelnumber))
    assert os.path.exists(path), "Unable to locate sample location: {}".format(path)
    s = open(path, "rb")
    with s:
        sample = pk.load(s)

    # Cross-check model number identifier
    assert sample.modelnumber == modelnumber, "Pickled model data does not match! {m1} != {m2}".format(
        m1=modelnumber, m2=sample.modelnumber
    )

    # Unpack background data
    back = sample.background
    Nend = sample.Nend
    assert Nend is not None, "Nend = None should not be a saved sample"
    kExit = sample.kExit
    pvals = np.asarray(sample.parameters)

    # Build savee path
    fNL_savepath = os.path.join(path3pf, "{num}.{ext}".format(num=modelnumber, ext=name))

    # Begin timer
    ti = time.clock()

    # Build Fourier triangle via Fergusson Shellard convention
    k1 = kExit / 2. - beta * kExit / 2.
    k2 = kExit * (1. + alpha + beta) / 4.
    k3 = kExit * (1. - alpha + beta) / 4.

    # Find smallest momentum: Build ICs from largest scale mode
    kmin = np.min([k1, k2, k3])
    
    # Compute three-point function

    # Compute initial conditions
    ICs = PyS.ICsBM(subevo, kmin, back, pvals, PyT)
    
    if ICs == (np.nan, np.nan):
        # Build null results dict, dump and return error data
        fNL_dict = {'{}'.format(name): None, 'time': None}
    
        # Save binary file
        fNL_file = open(fNL_savepath, "wb")
        with fNL_file:
            pk.dump(fNL_dict, fNL_file)
    
        return {"mn": modelnumber, "flag": -33, "ext": name}

    # Compute 3pt function
    threePt = PyT.alphaEvolve(np.array([ICs[0], Nend]), k1, k2, k3, ICs[1], pvals, tols, True, tmax_3pf, True)

    # If flag has been returned when computing 3pf, return flag data
    if type(threePt) is tuple:
        
        # Build null results dict, dump and return error data
        fNL_dict = {'{}'.format(name): None, 'time': None}
    
        # Save binary file
        fNL_file = open(fNL_savepath, "wb")
        with fNL_file:
            pk.dump(fNL_dict, fNL_file)
            
        return {"mn": modelnumber, "flag": threePt[0], "ext": name}
    
    # Transpose data into rows
    threePt = threePt.T

    # Get power spectra and bispectrum for triangule
    Pz1, Pz2, Pz3, Bz = [threePt[i][-1] for i in range(1, 5)]

    # Compute fNL
    fNL = (5. / 6.) * Bz / (Pz1 * Pz2 + Pz2 * Pz3 + Pz1 * Pz3)

    # Build result dictionary
    fNL_dict = {'{}'.format(name): fNL, 'time': time.clock() - ti}

    # Save binary file
    fNL_file = open(fNL_savepath, "wb")
    with fNL_file:
        pk.dump(fNL_dict, fNL_file)
        
    return 0
        

def computations(mn_calc):
    """ Simple handler for computations called from sampler """
    
    # Unpack model number and calculation type
    modelnumber, calculation = mn_calc

    # Prevent computation fails from hanging via catch_warnings
    with warnings.catch_warnings(record=True) as w:

        if calculation == "masses":
            print "-- Start MAb: model %06d" % modelnumber
            r = Mij(modelnumber)
            print "--   End MAb: model %06d" % modelnumber

        elif calculation == "2pf":
            print "-- Start 2pf: model %06d" % modelnumber
            r = SpectralIndex(modelnumber)
            print "--   End 2pf: model %06d" % modelnumber

        elif calculation not in ["masses", "2pf"]:
            print "-- Start fNL: model %06d, config {}".format(calculation) % modelnumber
            r = fNL(modelnumber, calculation)
            print "--   End fNL: model %06d, config {}".format(calculation) % modelnumber

        else:
            raise ValueError, "Undefined calculation: {}".format(calculation)
        
        if r != 0:
            
            # To do: Add error stats for mass-matrix. Though this should be stable.
            
            if calculation == "masses":
                pass
            else:
                
                # Determine subdirectory for error records
                if calculation == "2pf": subdir = calculation
                else: subdir = "3pf"
                
                # Get / remove extension key from dictionary
                pkEXT = r.pop('ext')
                
                # Dump dictionary object
                pkPath = os.path.join(pathStats, subdir, "{}.{}".format(modelnumber, pkEXT))
                
                pkFile = open(pkPath, "wb")
                with pkFile: pk.dump(r, pkFile)
                