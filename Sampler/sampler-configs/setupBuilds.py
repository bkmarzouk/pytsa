import os, sys, numpy as np, importlib

pytpath = (
    os.path.abspath(
        os.path.join(
            os.getcwd(), "../..", "PyTransport", "PyTransport"
        )
    )
)
assert os.path.exists(pytpath), "Cannot locate PyTransport installation: {}".format(pytpath)
sys.path.append(pytpath)

import PyTransSetup
PyTransSetup.pathSet()


class icpsCfgTemplate:
    
    def __init__(self, PyT):
        
        self.transportModule = PyT
        self.nF = PyT.nF()
        self.nP = PyT.nP()
        
        
    def __setValues__(self, attr_, command, fpNums, *requiredModules):
        
        # Locate number of expected value assignments for a given attribute
        if attr_ == "fields":
            if not hasattr(self, "fields"): self.fields = {}
            nFP = self.nF
        elif attr_ == "dotfields":
            if not hasattr(self, "dotfields"): self.dotfields = {}
            nFP = self.nF
        elif attr_ == "parameters":
            if not hasattr(self, "parameters"): self.parameters = {}
            nFP = self.nP
        else:
            raise AttributeError, "Attempting to set invalid attribute: {}".format(attr_)
        
        # Check command will make sense by interpreter
        assert type(command) in [
            str, float, int], "Write functional command in '' markers, or present a constant numeric value"
        
        # Try and iterate over attribute indices
        try:
            iter(fpNums)
            
            # If iterable
            for fN in fpNums:
                
                # Check sensible index number
                assert type(fN) is int, "Field/Parameter number must be integer: {}".format(fN)
                assert fN < nFP and fN >= 0, "Invalid field/parameter number: {}".format(fN)
                
                # Assign dictionary key value for index
                getattr(self, attr_)[str(fpNums)] = command
            
        # If not iterable
        except TypeError:
            
            # Check sensible idnex number
            assert type(fpNums) is int, "Field/Parameter number must be integer: {}".format(fpNums)
            assert fpNums < nFP and fpNums >= -1, "Invalid field/parameter number: {}".format(fpNums)

            # If -1 is passed as the index, assign all index keys the same command
            if fpNums == -1:
                for ii in range(nFP):
                    getattr(self, attr_)[str(ii)] = command
            
            # Otherwise assign the passed index the passed command
            else:
                getattr(self, attr_)[str(fpNums)] = command
                
        # If any required modules are passed
        for rM in requiredModules:
            
            # Check that the module name is given as a string
            assert type(rM) is str, "Pass any required modules in string format"
            
            # Initialize list for requiredModule attribute
            if not hasattr(self, "requiredModules"): self.requiredModules = []
            
            # Test that module passed is importable, if so, add to list
            try:
                importlib.import_module(rM)
                self.requiredModules.append(rM)
            except ImportError:
                raise ImportError, "Was unable to import: {}".format(rM)
            
            
    def setInitialFieldValues(self, fpNums, command, *requiredModules):
        
        self.__setValues__("fields", command, fpNums, *requiredModules)


    def setInitialFieldVelocities(self, fpNums, command, *requiredModules):
    
        self.__setValues__("dotfields", command, fpNums, *requiredModules)


    def setParameterValues(self, fpNums, command, *requiredModules):
    
        self.__setValues__("parameters", command, fpNums, *requiredModules)


class bispectrumCfgTemplate:
    
    def __init__(self):
        self.configurations = None
    
    
    def addBispectrumConfiguration(self, name, latex, alpha, beta):
        assert type(name) is str, "configuration name must be string"
        assert type(latex) is str, "latex definition must be string"
        assert type(alpha) is float, "alpha parameter must be float"
        assert type(beta) is float, "beta parameter must be float"
        
        if not hasattr(self, "configurations"): self.configurations = []
        
        self.configurations.append(
            {
                "name" : name,
                "latex": latex,
                "alpha": alpha,
                "beta" : beta
            }
        )


class PyTransportSampler(icpsCfgTemplate, bispectrumCfgTemplate):
    
    def __init__(self, PyT):
        
        self.transportModule = PyT
        self.nF = PyT.nF()
        self.nP = PyT.nP()
        
        self.name = None
        self.saveLocation = None
        self.efoldsBeforeExit = None
        self.subHorizonEvolution = None
        self.integratorTols = None
        self.adiabaticN = None
        
        self.goodExit = None
        self.badExit = None
        
        self.recordMasses = True
    
    
    def setCoreParams(self, name, saveLocation="default", efoldsBeforeExit=55., subHorizonEvolution=6.,
                      integratorTols=[1e-8, 1e-8], adiabaticN=1.):
        
        assert type(name) is str, "sampler name must be string to build directories"
        assert saveLocation == "default" or os.path.exists(saveLocation), "Path not found: {}".format(saveLocation)
        assert type(efoldsBeforeExit) in [float, int], "Specify N, s.t. Nend - N defines the pivot scale"
        assert type(subHorizonEvolution) in [float, int], "Specify the amount of subhorizon evolution"
        assert type(integratorTols) is list, "Specify a list of tols. for the integrator"
        assert len(integratorTols) == 2 and all([type(t) is float for t in integratorTols]), "set 2 floats"
        
        self.name = name
        self.saveLocation = saveLocation
        self.efoldsBeforeExit = efoldsBeforeExit
        self.subHorizonEvolution = subHorizonEvolution
        self.integratorTols = integratorTols
        self.adiabaticN = adiabaticN
    
    
    def addEndRegion(self, fieldNumber, minFieldValue, maxFieldValue):
        
        assert type(fieldNumber) is int, "field number must be an integer corresponding to the installation"
        assert type(minFieldValue) in [int, float], "minFieldValue must be an integer or float"
        assert type(maxFieldValue) in [int, float], "maxFieldValue must be an integer or float"
        
        if self.goodExit is None: self.goodExit = []
        
        self.goodExit.append(
            {
                "fieldNumber"  : fieldNumber,
                "minFieldValue": minFieldValue,
                "maxFieldValue": maxFieldValue
            }
        )
    
    
    def addViolatedRegion(self, fieldNumber, minFieldValue, maxFieldValue):
        
        assert type(fieldNumber) is int, "field number must be an integer corresponding to the installation"
        assert type(minFieldValue) in [int, float], "minFieldValue must be an integer or float"
        assert type(maxFieldValue) in [int, float], "maxFieldValue must be an integer or float"
        
        if self.badExit is None: self.badExit = []
        
        self.badExit.append(
            {
                "fieldNumber"  : fieldNumber,
                "minFieldValue": minFieldValue,
                "maxFieldValue": maxFieldValue
            }
        )
        
    
    def buildSampler(self):
    
    