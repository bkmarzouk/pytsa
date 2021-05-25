# # Import numerical tools
# import numpy as np
# from scipy.interpolate import UnivariateSpline
#
# # Import system tools
# import os
# import sys
# import importlib
# import pickle as pk
# import time
# import warnings
#
# # Retrieve paths from environment variables
# envKeys = [
#     'PyTS_pathPyT', 'PyTS_pathRoot', 'PyTS_pathSamples',
#     'PyTS_pathStats', 'PyTS_pathClasses', 'PyTS_pathMethods', 'PyTS_pathLocalData'
# ]
#
# pathPyT, pathRoot, \
# pathSamples, pathStats, pathClasses, pathMethods, pathLocalData = [os.environ[eK] for eK in envKeys]
#
# # get transport data dictionary
# transportPath = os.path.join(pathLocalData, "transport.localdata")
# transportFile = open(transportPath, "rb")
# with transportFile: transportDict = pk.load(transportFile)
#
# # configure internal pytransport paths
# sys.path.append(pathPyT)
# import PyTransSetup
#
# PyTransSetup.pathSet()
#
# # Import module and get other required data for computing observables
# PyT = importlib.import_module(transportDict['transportModule'])
# import PyTransScripts as PyS
#
# transportKeys = ["exitN", "subN", "intTols", "adiabaticN", "minN"]
# Nexit, subevo, tols, adiabaticN, minN = [transportDict[tK] for tK in transportKeys]
#
# # Get fNL data
# fNLPath = os.path.join(pathLocalData, "fNL.localdata")
# fNLFile = open(fNLPath, "rb")
# with fNLFile: fNLDict = pk.load(fNLFile)
#
# # Import sample generator & other useful class structures
# from generator import genSample
# from realization import realization as new_sample
#
# sys.path.append(pathMethods)
# from writer import update_sample

import numpy as np

from PyTransport.PyTransSetup import pathSet

pathSet()

import importlib


def compute_background(data):
    PyT = importlib.import_module(data.PyT)
    model_number = data.model_number
    catalogue_path = data.path
    ics_params = np.zeros((2 * PyT.nF() + PyT.nP()), dtype=np.float64)
    ics_params[::] = np.load(catalogue_path, mmap_mode="r")[model_number]

    ics = ics_params[:2 * PyT.nF()]
    params = ics_params[2 * PyT.nF():]

    tols = data.tols
    N_min = data.N_min
    N_sub = data.N_sub
    N_eps = PyT.findEndOfInflation(ics, params, tols, True)

    N_evo = np.linspace(0, N_eps, int(1.25 * N_eps))
    bg_eps = PyT.backEvolve(N_evo, ics, params, tols, True, -1, True)

    print(N_eps, bg_eps[-1][0])

#
# def buildICPs(modelnumber, rerun_model=False):
#     """ Builds initial conditions & parameters, or retrieves them if in rerun mode"""
#
#     if rerun_model is False:
#         # Generate new sample with parameter & field priori
#         n, fvals, vvals, pvals = genSample(modelnumber)
#
#     else:
#         # Load an existing set of parameters & model priori
#         model_path = os.path.join(saveloc, "samples", "{}.sample".format(modelnumber))
#         model_file = open(model_path, "rb")
#
#         with model_file as mf:
#             model = pk.load(mf)
#             fvals = model.fields
#             vvals = model.velocities
#             pvals = model.parameters
#
#             return np.concatenate(fvals, vvals), pvals
#
#     # If SR is in vvals, compute initial velocity value via slowroll equation
#     if "SR" in vvals:
#
#         # Compute intial velocities from approximation
#         vvalsSR = PyT.dotfieldsSR(fvals, pvals)
#
#         # Populate empty array with relevant values
#         vvals_ = np.zeros(PyT.nF())
#
#         # For each velocity def.
#         for jj in range(PyT.nF()):
#
#             #  If def. is SR, assign slow-roll velocity to float array
#             if vvals[jj] == "SR":
#                 vvals_[jj] = vvalsSR[jj]
#             else:
#                 # Otherwise force float-type to other array elements
#                 vvals_[jj] = float(vvals[jj])
#
#         # Redefine vvals
#         vvals = vvals_
#
#     # Concatenate field and velocity values into initial conditions array
#     return np.concatenate((fvals, vvals)), pvals
#
#
# def Initialize(modelnumber, rerun_model=False, initial=None, pvals=None, returnSample=False):
#     tmax_bg = int(os.environ['tmax_bg'])
#
#     if initial is None and pvals is None:
#         initial, pvals = buildICPs(modelnumber, rerun_model)
#
#     """ We have the following flags in place:
#
#     Timeout flags:
#     -10 : findEndOfInflation
#     -11 : backEvolve
#     -12 : sigEvolve
#     -13 : alphaEvolve
#
#     Integration error flags:
#     -20 : findEndOfInflation
#     -21 : beckEvolve
#     -22 : sigEvolve
#     -23 : alphaEvolve
#
#     Misc. PyT:
#     -30 : Nend < Nmin
#     -31 : k not found
#     -32 : No ICsBM (2pf)
#     -33 : No ICsBM (3pf)
#     -34 : Model violation
#     -35 : Eternal inflation
#
#     """
#
#     acceptIfAnyPath = os.path.join(pathLocalData, "acceptIfAny.localdata")
#     acceptIfAllPath = os.path.join(pathLocalData, "acceptIfAll.localdata")
#
#     rejectIfAnyPath = os.path.join(pathLocalData, "rejectIfAny.localdata")
#     rejectIfAllPath = os.path.join(pathLocalData, "rejectIfAll.localdata")
#
#     # We define non-canonical conditions to end inflation, and those which violate inflation
#     endAnyNC = os.path.exists(acceptIfAnyPath)  # NC end
#     endAllNC = os.path.exists(acceptIfAllPath)
#
#     exitAnyNC = os.path.exists(rejectIfAnyPath)  # NC violate
#     exitAllNC = os.path.exists(rejectIfAllPath)
#
#     canonicalEnd = not (endAnyNC or endAllNC)
#     canonicalExit = not (exitAnyNC or exitAllNC)
#
#     # if canonical end or not, suitable to start with cheap epsilon search
#     Nepsilon = PyT.findEndOfInflation(initial, pvals, tols, 0.0, 10000, tmax_bg, True)
#
#     # First we check models that end with SR violation
#     if canonicalEnd:
#
#         # Check for integration / time out error
#         if type(Nepsilon) is tuple:
#             return Nepsilon
#
#         # Check for enough inflation
#         elif Nepsilon < minN:
#             return (-30, Nepsilon)
#
#         # Looks good, go compute background
#         else:
#             back = PyT.backEvolve(
#                 np.linspace(0, Nepsilon, int(3 * Nepsilon)), initial, pvals, tols, True, tmax_bg, True)
#
#         Nend = Nepsilon
#
#     # For models that don't end with SR violation
#     else:
#
#         # Check for integration / time out error
#         if type(Nepsilon) is tuple:
#             # Instead, draw back farthest efold reached by 10% and try and compute fid. background
#             Nepsilon = 0.9 * Nepsilon[1]
#
#         # Try and obtain fiducial background up until epsilon = 1
#         backFid = PyT.backEvolve(
#             np.linspace(0, Nepsilon, 3 * (int(Nepsilon) + 1)), initial, pvals, tols, 0, tmax_bg, True)
#
#         """ We will now attempt to extend the background by up to 100 efolds passed the epsilon definition
#             Note that getting further than this is highly unlikely, since integration becomes hard
#             Note that we increase desired tols by an order of magnitude to try and achieve this """
#
#         backExt = PyT.backEvolve(
#             np.linspace(backFid[-1][0], backFid[-1][0] + 100, 300), backFid[-1][1:], pvals, tols * 1e-1, 0, tmax_bg,
#             True)
#
#         # Combine arrays, omiting the first row of the extension since this will match the last of fid.
#         back = np.vstack((backFid, backExt[1:]))
#
#         # We now will search through the background and look for Nend condition
#         endIndex = len(back)
#         foundEnd = False
#
#         # Load conditions
#         if endAnyNC:
#             with open(acceptIfAnyPath, "rb") as f: acceptIfAny = pk.load(f)
#             acceptIfAnyIdx = acceptIfAny['idx']
#             acceptIfAnyMin = acceptIfAny['min']
#             acceptIfAnyMax = acceptIfAny['max']
#
#         if endAllNC:
#             with open(acceptIfAllPath, "rb") as f: acceptIfAll = pk.load(f)
#             acceptIfAllIdx = acceptIfAll['idx']
#             acceptIfAllMin = acceptIfAll['min']
#             acceptIfAllMax = acceptIfAll['max']
#
#         for idx in range(len(back)):
#
#             row = back[idx]
#
#             if endAnyNC and np.any(
#                     np.logical_and(row[acceptIfAnyIdx] > acceptIfAnyMin, row[acceptIfAnyIdx] < acceptIfAnyMax)):
#                 endIndex = idx
#                 Nend = row[0]
#                 foundEnd = True
#                 break
#
#             if endAllNC and np.all(
#                     np.logical_and(row[acceptIfAllIdx] > acceptIfAllMin, row[acceptIfAllIdx] < acceptIfAllMax)):
#                 endIndex = idx
#                 Nend = row[0]
#                 foundEnd = True
#                 break
#
#         endIndex = min(endIndex + 1, len(back))
#
#         back = back[:endIndex]
#
#     if canonicalExit is False:
#
#         if exitAnyNC:
#             with open(rejectIfAnyPath, "rb") as f: rejectIfAnyPath = pk.load(f)
#             rejectIfAnyIdx = rejectIfAnyPath['idx']
#             rejectIfAnyMin = rejectIfAnyPath['min']
#             rejectIfAnyMax = rejectIfAnyPath['max']
#
#         if exitAllNC:
#             with open(rejectIfAllPath, "rb") as f: rejectIfAllPath = pk.load(f)
#             rejectIfAllIdx = rejectIfAllPath['idx']
#             rejectIfAllMin = rejectIfAllPath['min']
#             rejectIfAllMax = rejectIfAllPath['max']
#
#         for row in back:
#
#             if exitAnyNC and np.any(
#                     np.logical_and(row[rejectIfAnyIdx] > rejectIfAnyMin, row[rejectIfAnyIdx] < rejectIfAnyMax)):
#                 Nend = row[0]
#
#                 return -34, Nend
#
#             if exitAllNC and np.all(
#                     np.logical_and(row[rejectIfAllIdx] > rejectIfAllMin, row[rejectIfAllIdx] < rejectIfAllMax)):
#                 Nend = row[0]
#
#                 return -34, Nend
#
#     # Now perform the final check that the background
#     if canonicalEnd:
#         if Nend < minN:
#             return -30, Nend
#     else:
#
#         Nend = back[-1][0]
#
#         if Nend > minN and not foundEnd:
#             return -35, Nend
#
#         if Nend < minN and foundEnd:
#             return -30, Nend
#
#         if Nend < minN and not foundEnd:
#             return -21, Nend
#
#     # Track *actual* initial conditions prior to repositioning
#     initial_0 = initial
#
#     # Rescale background to start at Nend - minN to avoid exp. large k
#     back_adj = PyS.rescaleBack(back, Nr=minN + 10)
#
#     # Attempt computation of momenta at horizon crossing, the value used will be in the ballpark of requirements
#     # for 2pf and 3pf tasks. Hence, we discard the realization at this point if it's prospects look bad.
#
#     Npiv = PyS.matchKExitN(back_adj, pvals, PyT, 0.002)
#
#     if np.isnan(Npiv):
#         return -31
#
#     kExit = PyS.kexitN(Nend - Npiv, back_adj, pvals, PyT)
#
#     # Asses success by data type of momenta result
#     if type(kExit) not in [float, np.float, np.float32, np.float64] or np.isnan(kExit) or np.isinf(kExit):
#         return -31
#
#     # We now test for an adiabatic limit: We begin by assuming this is True, then test for violating conditions
#     adiabatic = True
#
#     # Define number of efolds from end of inflation which should be adiabatic TODO: Wrap this into config. file
#
#     # Find first efoldign to perform eigenvalue test
#     try:
#         adiabaticN_start = np.where(back.T[0] >= back_adj[-1][0] - adiabaticN)[0][0]
#     except IndexError:
#         return -21
#
#     # Compute mass-matrix evolution from this point
#     Mij_end = PyS.evolveMasses(back[adiabaticN_start:], pvals, PyT, scale_eigs=False, hess_approx=False,
#                                covariant=False)
#
#     # For each mass-matrix evolution step
#     for item in Mij_end:
#
#         # Determine if lightest mass is tachyonic
#         tachyon = item[1] < 0
#
#         # Determine number of hubble scale heavier masses
#         rest_large = sum([eig >= 1 for eig in item[2:]]) == PyT.nF() - 1
#
#         # If there is no tachyon or the remaining masses aren't large for the portion of evolution
#         if not (tachyon and rest_large):
#             # Change adiabatic property to False, and break out of search window
#             adiabatic = False
#
#             break
#
#     # Record all sample data
#     new_sample(modelnumber, initial_0[:PyT.nF()], initial_0[PyT.nF():], pvals, back_adj, adiabatic,
#                Nend, kExit, pathSamples, rerun_model)
#
#     print
#     "-- Generated sample: %05d" % modelnumber
#
#     assert type(modelnumber) == int, modelnumber
#
#     # Return model number, acts as flag to confirm successful sample
#     return modelnumber
#
#
# def updateStats(modelnumber, flagKey=None, timeEnd=None):
#     """ Rewrites stats object during demand sample routine """
#
#     # Get directory for sample stats log
#     sample_path = os.path.join(pathStats, "bg", "{}.bg".format(modelnumber))
#
#     if os.path.exists(sample_path):
#         with open(sample_path, "r") as f:
#             flagDict = pk.load(f)
#         pathExists = True
#
#     else:
#         pathExists = False
#
#         # Flag defs.
#         timeoutFlags = [
#             -10, -11, -12, -13  # end of inflation, background, 2pf, 3pf
#         ]
#
#         integratorFlags = [
#             -20, -21, -22, -23  # end of inflation, background, 2pf, 3pf
#         ]
#
#         samplerFlags = [
#             -30, -31, -32, -33, -34, -35, -36  # N < Nmin, k not found, no ICs 2pf, no ICs 3pf, model violation, eternal
#         ]
#
#         allFlags = timeoutFlags + integratorFlags + samplerFlags
#         flagDict = {str(f): 0 for f in allFlags}
#
#     if flagKey in flagDict:
#         flagDict[flagKey] += 1
#     elif flagKey is None and timeEnd is not None:
#         flagDict['time'] = timeEnd
#     else:
#         raise KeyError, flagKey
#
#     if pathExists:
#         os.remove(sample_path)
#
#     with open(sample_path, "w") as f:
#         pk.dump(flagDict, f)
#
#
# def DemandSample(modelnumber, rerun=False):
#     """ Repeat initialization until successful sample is found """
#
#     # # Get directory for sample stats log
#     # sample_path = os.path.join(pathStats, "bg", "{}.bg".format(modelnumber))
#     #
#     # Flag defs.
#     timeoutFlags = [
#         -10, -11, -12, -13  # end of inflation, background, 2pf, 3pf
#     ]
#
#     integratorFlags = [
#         -20, -21, -22, -23  # end of inflation, background, 2pf, 3pf
#     ]
#
#     samplerFlags = [
#         -30, -31, -32, -33, -34, -35, -36
#         # N < Nmin, k not found, no ICs 2pf, no ICs 3pf, model violation, eternal, mij eig fail
#     ]
#
#     allFlags = timeoutFlags + integratorFlags + samplerFlags
#     # flagDict = {str(f): 0 for f in allFlags}
#
#     # Start timer
#     tstart = time.clock()
#
#     ii = -1
#
#     # If successful sample generation, the model number is returned
#     while ii != modelnumber:
#
#         ii = Initialize(modelnumber, rerun_model=rerun)
#
#         # If fail flag received: log statistic
#         if type(ii) is tuple:
#             flagKey = str(ii[0])
#             timeEnd = None
#             # flagDict[flagKey] += 1
#
#         elif type(ii) is int and ii in allFlags:
#             flagKey = str(ii)
#             timeEnd = None
#             # flagDict[flagKey] += 1
#
#         # If model number, sum number of iterations
#         elif ii == modelnumber:
#             flagKey = None
#             timeEnd = time.clock() - tstart
#
#         else:
#             raise ValueError, modelnumber
#
#         updateStats(modelnumber, flagKey=flagKey, timeEnd=timeEnd)
#
#
# def DemandSample_rerun(modelnumber):
#     """
#     Repeat initialization until successful sample is found : using rerun samples
#     """
#     DemandSample(modelnumber, rerun=True)
#
#
# def Mij(modelnumber):
#     """
#     Computes the eigenvalues of the mass-square-matrix at horizon crossing
#     """
#
#     # Unload sample data
#     path = os.path.join(pathSamples, "{}.sample".format(modelnumber))
#     assert os.path.exists(path), "Unable to locate sample location: {}".format(path)
#     s = open(path, "rb")
#     with s:
#         sample = pk.load(s)
#
#     # Cross-check model number identifier
#     assert sample.modelnumber == modelnumber, "Pickled model data does not match! {m1} != {m2}".format(
#         m1=modelnumber, m2=sample.modelnumber
#     )
#
#     # Unpack background data
#     back = sample.background
#     Nend = sample.Nend
#     assert Nend is not None, "Nend = None -> Sample construction failure!"
#     pvals = np.asarray(sample.parameters)
#
#     # Transport background for col evolution rather than row evolution
#     backT = back.T
#
#     # Get number of fields
#     nF = PyT.nF()
#
#     # Begin timer
#     ti = time.clock()
#
#     # Define Nstar, i.e. exit time of interest
#     Nstar = Nend - Nexit
#
#     # Unpack efold evolution
#     Nevo = backT[0]
#
#     # If Nstar falls exactly in the background, simply load this background step
#     if Nstar in Nevo:
#         Nstar_idx = [Nevo.index(Nstar)]
#         back_star = back[Nstar_idx]
#
#     # Otherwise we build a spline about Nstar
#     else:
#
#         # Get background index values +/- 1 efold around the prescribed exit time
#         Nfit = 0
#         dN = 0.5
#         while Nfit < 3:
#             dN += 0.5
#             idx_above = np.where(np.logical_and(Nevo > Nstar, Nevo < Nstar + dN))[0]
#             idx_below = np.where(np.logical_and(Nevo < Nstar, Nevo > Nstar - dN))[0]
#
#             Nstar_idx = list(np.concatenate([idx_below, idx_above]))
#             Nfit = len(Nstar_idx)
#
#         # Define efold range to "smooth" over
#         Nsmooth = [Nevo[idx] for idx in Nstar_idx]
#
#         # We check for repeated efold numbers as this will throw the interpolation scheme
#         repLog = []
#         for ii in range(len(Nsmooth) - 1):
#             if Nsmooth[ii] == Nsmooth[ii + 1]:
#                 repLog.append(ii + 1)
#         for index in sorted(repLog, reverse=True):
#             del Nstar_idx[index]
#
#         Nsmooth = [Nevo[idx] for idx in Nstar_idx]
#
#         # BY default, we choose a cubic spline, but if this is not possible, we reduce the order
#         if len(Nstar_idx) <= 3:
#             k = len(Nstar_idx) - 1
#         else:
#             k = 3
#
#         # Construct splines for the field evolution over Nsmooth
#         fsplines = [
#             UnivariateSpline(
#                 Nsmooth, [backT[i][idx] for idx in Nstar_idx], k=k) for i in range(1, 1 + nF)
#         ]
#
#         # Construct splines for the velocity evolution over Nsmooth
#         vsplines = [
#             UnivariateSpline(
#                 Nsmooth, [backT[i][idx] for idx in Nstar_idx], k=k) for i in range(1 + nF, 1 + 2 * nF)
#         ]
#
#         # Evaluate the background at N_star = Nend - Nexit
#         back_star = np.concatenate([np.array([Nstar]), np.array([s(Nstar) for s in fsplines + vsplines])])
#
#     try:
#         # Compute mass matrix eigenvalues at the exit time
#         Mij = PyS.evolveMasses(
#             np.array([back_star]), pvals, PyT, scale_eigs=False, hess_approx=False, covariant=False)[0][1:]
#
#     except np.linalg.LinAlgError:
#
#         # Build dictionary of masses
#         Mij_eig = {}
#         for i in range(PyT.nF()):
#             Mij_eig['m{}'.format(i)] = None
#
#         Mij_eig['T_masses'] = None
#
#         sample.update_observables(Mij_eig)
#
#         return {"mn": modelnumber, "flag": -36, "ext": "mij"}
#
#     # Build dictionary of masses
#     Mij_eig = {}
#     for i in range(len(Mij)):
#         Mij_eig['m{}'.format(i)] = Mij[i]
#
#     # End timer
#     Mij_eig['T_masses'] = time.clock() - ti
#
#     # Update sample observables
#     sample.update_observables(Mij_eig)
#
#     return 0
#
#
# def SpectralIndex(modelnumber):
#     # Get timeout
#     tmax_2pf = int(os.environ['tmax_2pf'])
#
#     # Unload sample data
#     path = os.path.join(pathSamples, "{}.sample".format(modelnumber))
#     assert os.path.exists(path), "Unable to locate sample location: {}".format(path)
#     s = open(path, "rb")
#     with s:
#         sample = pk.load(s)
#
#     # Cross-check model number identifier
#     assert sample.modelnumber == modelnumber, "Pickled model data does not match! {m1} != {m2}".format(
#         m1=modelnumber, m2=sample.modelnumber
#     )
#
#     # def spectralIndex(back, pvals, Nexit, tols, subevo, MTE, kPivot=None, returnRunning=True, errorReturn=False, tmax=None):
#     ti = time.clock()
#     ns_alpha = PyS.spectralIndex(
#         sample.background, sample.parameters, Nexit, tols, subevo, PyT, errorReturn=True, tmax=tmax_2pf)
#
#     # Check for failed ns calculation
#     if isinstance(ns_alpha, tuple):
#         if ns_alpha[0] is ValueError:
#
#             errKey = ns_alpha[1]
#
#             # Define error dictionary to return
#             if errKey == "k":
#                 retDict = {"mn": modelnumber, "flag": -31, "ext": "2pf"}
#             elif errKey == "ics":
#                 retDict = {"mn": modelnumber, "flag": -32, "ext": "2pf"}
#             elif errKey == "2pf":
#                 retDict = {"mn": modelnumber, "flag": ns_alpha[2][0], "ext": "2pf"}
#             else:
#                 raise KeyError, "Unknown error: {}".format(errKey)
#
#             # Update sample with null result
#             sample.update_observables({'ns': None, 'alpha': None, 'T_2pf': None})
#
#             return retDict
#
#     # Unpack 2pt observables
#     ns, alpha = ns_alpha
#
#     # End calculation timer and build results dictionaryobs_
#     twoPt_dict = {'ns': ns, 'alpha': alpha, 'T_2pf': time.clock() - ti}
#
#     # Update sample observables
#     sample.update_observables(twoPt_dict)
#
#     return 0
#
#
# def fNL(modelnumber, configName):
#     tmax_3pf = int(os.environ['tmax_3pf'])
#
#     name = None
#     # Gather configuration data from configName
#     for d in fNLDict:
#         if d['name'] == configName:
#             name = configName
#             alpha = d['alpha']
#             beta = d['beta']
#             break
#
#     # Check that fNL definition has been found
#     assert name is not None, "Unable to locate fNL configuration in localdata: {}".format(name)
#
#     # Unload sample data
#     path = os.path.join(pathSamples, "{}.sample".format(modelnumber))
#     assert os.path.exists(path), "Unable to locate sample location: {}".format(path)
#     s = open(path, "rb")
#     with s:
#         sample = pk.load(s)
#
#     # Cross-check model number identifier
#     assert sample.modelnumber == modelnumber, "Pickled model data does not match! {m1} != {m2}".format(
#         m1=modelnumber, m2=sample.modelnumber
#     )
#     # Begin timer
#     ti = time.clock()
#
#     # Compute fNL
#     fNL = PyS.fNL(
#         sample.background, sample.parameters, Nexit, tols, subevo, PyT, alpha=alpha, beta=beta, stdConfig=name,
#         errorReturn=True, tmax=tmax_3pf)
#
#     # Check for failed fNL calculation calculation
#     if isinstance(fNL, tuple):
#         if fNL[0] is ValueError:
#
#             errKey = fNL[1]
#
#             # Define error dictionary to return
#             if errKey == "k":
#                 retDict = {"mn": modelnumber, "flag": -31, "ext": name}
#             elif errKey == "ics":
#                 retDict = {"mn": modelnumber, "flag": -33, "ext": name}
#             elif errKey == "3pf":
#                 retDict = {"mn": modelnumber, "flag": fNL[2][0], "ext": name}
#             else:
#                 raise KeyError, "Unknown error: {}".format(errKey)
#
#             # Update sample with null result
#             sample.update_observables({'{}'.format(name): None, 'T_{}'.format(configName): None})
#
#             return retDict
#
#     # Build result dictionary
#     threePt_Dict = {'{}'.format(name): fNL, 'T_{}'.format(configName): time.clock() - ti}
#
#     # Update sample observables
#     sample.update_observables(threePt_Dict)
#
#     return 0
#
#
# def computations(mn_calc):
#     """ Simple handler for computations called from sampler """
#
#     # Unpack model number and calculation type
#     modelnumber, calculation = mn_calc
#
#     # Prevent computation fails from hanging via catch_warnings
#     with warnings.catch_warnings(record=True):
#
#         if calculation == "masses":
#             print
#             "-- Start MAb: model %06d" % modelnumber
#             r = Mij(modelnumber)
#             status = "" if r == 0 else ", failed"
#             print
#             "--   End MAb: model %06d" % modelnumber + status
#
#         elif calculation == "2pf":
#             print
#             "-- Start 2pf: model %06d" % modelnumber
#             r = SpectralIndex(modelnumber)
#             status = "" if r == 0 else ", failed"
#             print
#             "--   End 2pf: model %06d" % modelnumber + status
#
#         elif calculation not in ["masses", "2pf"]:
#             print
#             "-- Start fNL: model %06d, config {}".format(calculation) % modelnumber
#             r = fNL(modelnumber, calculation)
#             status = "" if r == 0 else ", failed"
#             print
#             "--   End fNL: model %06d, config {}".format(calculation) % modelnumber + status
#
#         else:
#             raise ValueError, "Undefined calculation: {}".format(calculation)
#
#         if r != 0:
#
#             # Determine subdirectory for error records
#             # For now, we will store masses fails in 2pf dir to avoid restructuring
#             if calculation == "masses":
#                 subdir = "mij"
#             elif calculation == "2pf":
#                 subdir = "2pf"
#             else:
#                 subdir = "3pf"
#
#             # Get / remove extension key from dictionary, used to discriminate 3pf varieties
#             pkEXT = r.pop('ext')
#
#             # Dump dictionary object
#             pkPath = os.path.join(pathStats, subdir, "{}.{}".format(modelnumber, pkEXT))
#
#             pkFile = open(pkPath, "wb")
#             with pkFile:
#                 pk.dump(r, pkFile)
