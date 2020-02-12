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


def buildICPs(modelnumber, rerun_model=False):
    """ Builds initial conditions & parameters, or retrieves them if in rerun mode"""
    
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
            
            return np.concatenate(fvals, vvals)
    
    # If SR is in vvals, compute initial velocity value via slowroll equation
    if "SR" in vvals:
        
        # Compute intial velocities from approximation
        vvalsSR = PyT.dotfieldsSR(fvals, pvals)
        
        # Populate empty array with relevant values
        vvals_ = np.zeros(PyT.nF())
        
        # For each velocity def.
        for jj in range(PyT.nF()):
            
            #  If def. is SR, assign slow-roll velocity to float array
            if vvals[jj] == "SR":
                vvals_[jj] = vvalsSR[jj]
            else:
                # Otherwise force float-type to other array elements
                vvals_[jj] = float(vvals[jj])
        
        # Redefine vvals
        vvals = vvals_
        
    # Concatenate field and velocity values into initial conditions array
    return np.concatenate((fvals, vvals)), pvals


def Initialize(modelnumber, rerun_model=False):
    
    tmax_bg = int(os.environ['tmax_bg'])
    
    initial, pvals = buildICPs(modelnumber, rerun_model)

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
    acceptIfAllPath = os.path.join(pathLocalData, "acceptIfAll.localdata")
    
    rejectIfAnyPath = os.path.join(pathLocalData, "rejectIfAny.localdata")
    rejectIfAllPath = os.path.join(pathLocalData, "rejectIfAll.localdata")
    
    # We define non-canonical conditions to end inflation, and those which violate inflation
    endAnyNC = os.path.exists(acceptIfAnyPath) # NC end
    endAllNC = os.path.exists(acceptIfAllPath)
    
    exitAnyNC = os.path.exists(rejectIfAnyPath) # NC violate
    exitAllNC = os.path.exists(rejectIfAllPath)
    
    canonicalEnd = not (endAnyNC or endAllNC)
    canonicalExit = not (exitAnyNC or exitAllNC)

    # if canonical end or not, suitable to start with cheap epsilon search
    Nepsilon = PyT.findEndOfInflation(initial, pvals, tols, 0.0, 10000, tmax_bg, True)

    # First we check models that end with SR violation
    if canonicalEnd:
        
        # Check for integration / time out error
        if type(Nepsilon) is tuple:
            return Nepsilon
        
        # Check for enough inflation
        elif Nepsilon < minN:
            return (-30, Nepsilon)
        
        # Looks good, go compute background
        else:
            back = PyT.backEvolve(
                np.linspace(0, Nepsilon, int(3*Nepsilon)), initial, pvals, tols, True, tmax_bg, True)
            
        Nend = Nepsilon

    # For models that don't end with SR violation
    else:
    
        # Check for integration / time out error
        if type(Nepsilon) is tuple:
            
            # Instead, draw back farthest efold reached by 10% and try and compute fid. background
            Nepsilon = 0.9*Nepsilon[1]
            
        # Try and obtain fiducial background up until epsilon = 1
        backFid = PyT.backEvolve(
            np.linspace(0, Nepsilon, 3*(int(Nepsilon)+1)), initial, pvals, tols, 0, tmax_bg, True)
        
        """ We will now attempt to extend the background by up to 100 efolds passed the epsilon definition
            Note that getting further than this is highly unlikely, since integration becomes hard
            Note that we increase desired tols by an order of magnitude to try and achieve this """

        backExt = PyT.backEvolve(
            np.linspace(backFid[-1][0], backFid[-1][0]+100, 300), backFid[-1][1:], pvals, tols*1e-1, 0, tmax_bg, True)
        
        # Combine arrays, omiting the first row of the extension since this will match the last of fid.
        back = np.vstack((backFid, backExt[1:]))

        
        # We now will search through the background and look for Nend condition
        endIndex = len(back)
        foundEnd = False
        
        # Load conditions
        if endAnyNC:
            with open(acceptIfAnyPath, "rb") as f: acceptIfAny = pk.load(f)
            acceptIfAnyIdx = acceptIfAny['idx']
            acceptIfAnyMin = acceptIfAny['min']
            acceptIfAnyMax = acceptIfAny['max']

        if endAllNC:
            with open(acceptIfAllPath, "rb") as f: acceptIfAll = pk.load(f)
            acceptIfAllIdx = acceptIfAll['idx']
            acceptIfAllMin = acceptIfAll['min']
            acceptIfAllMax = acceptIfAll['max']
        
        for idx in range(len(back)):
            
            row = back[idx]
            
            if endAnyNC and np.any(
                    np.logical_and(row[acceptIfAnyIdx] > acceptIfAnyMin, row[acceptIfAnyIdx] < acceptIfAnyMax)):
                
                endIndex = idx
                Nend = row[0]
                foundEnd = True
                break

            if endAllNC and np.all(
                    np.logical_and(row[acceptIfAllIdx] > acceptIfAllMin, row[acceptIfAllIdx] < acceptIfAllMax)):
                endIndex = idx
                Nend = row[0]
                foundEnd = True
                break
    
        back = back[:endIndex]
        
    if canonicalExit is False:
        
        if exitAnyNC:
            with open(rejectIfAnyPath, "rb") as f: rejectIfAnyPath = pk.load(f)
            rejectIfAnyIdx = rejectIfAnyPath['idx']
            rejectIfAnyMin = rejectIfAnyPath['min']
            rejectIfAnyMax = rejectIfAnyPath['max']

        if exitAllNC:
            with open(rejectIfAllPath, "rb") as f: rejectIfAllPath = pk.load(f)
            rejectIfAllIdx = rejectIfAllPath['idx']
            rejectIfAllMin = rejectIfAllPath['min']
            rejectIfAllMax = rejectIfAllPath['max']
            
        for row in back:
    
            if exitAnyNC and np.any(
                    np.logical_and(row[rejectIfAnyIdx] > rejectIfAnyMin, row[rejectIfAnyIdx] < rejectIfAnyMax)):
                Nend = row[0]
                
                return -34, Nend
                
            if exitAllNC and np.all(
                    np.logical_and(row[rejectIfAllIdx] > rejectIfAllMin, row[rejectIfAllIdx] < rejectIfAllMax)):
                Nend = row[0]
                
                return -34, Nend
            
    # Now perform the final check that the background
    if canonicalEnd:
        if Nend < minN:
            return -30, Nend
    else:
        
        Nend = back[-1][0]
        
        if Nend > minN and not foundEnd:
            return -35, Nend
        
        if Nend < minN and foundEnd:
            return -30, Nend
        
        if Nend < minN and not foundEnd:
            return -21, Nend

    # Track *actual* initial conditions prior to repositioning
    initial_0 = initial
    
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
        rest_large = sum([eig >= 1 for eig in item[2:]]) == PyT.nF() - 1

        # If there is no tachyon or the remaining masses aren't large for the portion of evolution
        if not (tachyon and rest_large):
            
            # Change adiabatic property to False, and break out of search window
            adiabatic = False
            
            break

    # Record all sample data
    new_sample(modelnumber, initial_0[:PyT.nF()], initial_0[PyT.nF():], pvals, back_adj, adiabatic,
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

    # We will compute the 2pf based on momenta that horizon exit 1.2-efolds about k pivot
    Npivot = Nend - Nexit
    Nscatter = 0.6
    Nsteps = 2
    DeltaN = [Npivot - Nscatter + ii*Nscatter/float(Nsteps) for ii in range(Nsteps*2 + 1)]

    kVals = np.array([])
    for NN in DeltaN:
        
        k = PyS.kexitN(NN, back, pvals, PyT)
        
        if np.isinf(k) or np.isnan(k):
            
            # Build null result file, dump and return error data
            twoPt_dict = {'ns': None, 'alpha': None, 'time': None}
    
            # Write binary file
            twoPt_file = open(twoPt_savepath, "wb")
            with twoPt_file:
                pk.dump(twoPt_dict, twoPt_file)
    
            return {"mn": modelnumber, "flag": -31, "ext": "2pf"}
        
        kVals = np.append(kVals, k)

    # Check momenta are monotonically increasing
    assert all(kVals[i]<kVals[i+1] for i in range(len(kVals)-1)), kVals

    # Get pivot scale
    kPivot = kVals[Nsteps]
    
    # Build power spectra values
    pZetaVals = np.array([])
    for k in kVals:
        
        # Get ICs from massless condition
        Nstart, ICs = PyS.ICsBM(subevo, k, back, pvals, PyT)
        
        if type(ICs) != np.ndarray and np.isnan(ICs):
            
            # Build null result file, dump and return error data
            twoPt_dict = {'ns': None, 'alpha': None, 'time': None}
    
            # Write binary file
            twoPt_file = open(twoPt_savepath, "wb")
            with twoPt_file:
                pk.dump(twoPt_dict, twoPt_file)
                
            return {"mn": modelnumber, "flag": -32, "ext": "2pf"}

        twoPf = PyT.sigEvolve(np.array([Nstart, Nend]), k, ICs, pvals, tols, False, tmax_2pf, True)
        
        if type(twoPf) is tuple:
            
            # Build null result file, dump and return error data
            twoPt_dict = {'ns': None, 'alpha': None, 'time': None}
    
            # Write binary file
            twoPt_file = open(twoPt_savepath, "wb")
            with twoPt_file:
                pk.dump(twoPt_dict, twoPt_file)
                
            return {"mn": modelnumber, "flag": twoPf[0], "ext": "2pf"}
        
        pZetaVals = np.append(pZetaVals, twoPf.T[1][-1])
        
    twoPtSpline = UnivariateSpline(
        np.log(kVals/kPivot), np.log(pZetaVals), k=4
    )
    
    nsSpline = twoPtSpline.derivative()
    alphaSpline = nsSpline.derivative()
    
    ns_ = lambda k: nsSpline(np.log(k/kPivot)) + 4.0
    alpha_ = lambda k: alphaSpline(np.log(k/kPivot))
    
    ns = ns_(kPivot)
    alpha = alpha_(kPivot)

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
    Nstart, ICs = PyS.ICsBM(subevo, kmin, back, pvals, PyT)
    
    if type(ICs) != np.ndarray and np.isnan(ICs):
        # Build null results dict, dump and return error data
        fNL_dict = {'{}'.format(name): None, 'time': None}
    
        # Save binary file
        fNL_file = open(fNL_savepath, "wb")
        with fNL_file:
            pk.dump(fNL_dict, fNL_file)
    
        return {"mn": modelnumber, "flag": -33, "ext": name}

    # Compute 3pt function
    threePt = PyT.alphaEvolve(np.array([Nstart, Nend]), k1, k2, k3, ICs, pvals, tols, True, tmax_3pf, True)

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
                

