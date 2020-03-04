import os, sys, numpy as np, importlib, shutil, pickle as pk

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
            str, float, np.float, np.float32, np.float64, int], "Write functional command in '' markers, or present a constant numeric value"
        
        # Try and iterate over attribute indices
        try:
            iter(fpNums)
            
            # If iterable
            for fN in fpNums:
                
                # Check sensible index number
                assert type(fN) is int, "Field/Parameter number must be integer: {}".format(fN)
                assert fN < nFP and fN >= 0, "Invalid field/parameter number: {}".format(fN)
                
                # Assign dictionary key value for index
                getattr(self, attr_)[str(fN)] = command
            
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
            
        # Finally check that the command can be evaluated
        if command == "SR" and attr_ == "dotfields":
            pass
        
        else:
            try:
                
                # Execute import statements
                if hasattr(self, "requiredModules"):
                    for m in self.requiredModules:
                        exec("import {}".format(m))
                
                # Attempt evaluation of result
                constTypes = [int, float, np.float, np.float32, np.float64]
                
                result = command if type(command) in constTypes else eval(command)
                
                # Check result data type
                assert type(result) in constTypes, "Invalid result: {} -> {}".format(command, result)
                
            except:
                raise ValueError, "Could not evaluate command: {}".format(command)
            
            
    def setInitialFieldValues(self, fpNums, command, *requiredModules):
        
        self.__setValues__("fields", command, fpNums, *requiredModules)


    def setInitialFieldVelocities(self, fpNums, command, *requiredModules):
    
        self.__setValues__("dotfields", command, fpNums, *requiredModules)


    def setParameterValues(self, fpNums, command, *requiredModules):
    
        self.__setValues__("parameters", command, fpNums, *requiredModules)
    
    
    def __setLaTeX__(self, attr_, num, LaTeX):
        
        if not hasattr(self, "latex"): self.latex = {}
        
        # Locate number of expected value assignments for a given attribute
        if attr_ == "fields":
            prefix = "f"
            assert type(num) is int and num in range(0, self.nF), "Invalid field index: {}".format(self.nF)
        elif attr_ == "dotfields":
            prefix = "v"
            assert type(num) is int and num in range(0, self.nF), "Invalid field index: {}".format(self.nF)
        elif attr_ == "parameters":
            prefix = "p"
            assert type(num) is int and num in range(0, self.nP), "Invalid param index: {}".format(self.nP)
        else:
            raise AttributeError, "Attempting to set invalid attribute: {}".format(attr_)

        key = prefix + str(num)
        
        assert key not in self.latex, "LaTeX definition already prescribed: {} -> {}".format(key, self.latex[key])
        
        assert type(LaTeX) is str, "LaTeX definition must be string: {}".format(LaTeX)

        self.latex[key] = LaTeX
        

    def recordFieldValue(self, fieldIndex, latex):
        self.__setLaTeX__("fields", fieldIndex, latex)
        
    
    def recordDotFieldValue(self, fieldIndex, latex):
        self.__setLaTeX__("dotfields", fieldIndex, latex)
        
        
    def recordParameterValue(self, parameterIndex, latex):
        self.__setLaTeX__("parameters", parameterIndex, latex)
        
        
    def checkicps(self):
        
        field_range = []
        dotfield_range = []
        parameter_range = []
        
        """TODO: Iterable definitions not working!! """
        for key in self.fields: field_range.append(int(key))
        for key in self.dotfields: dotfield_range.append(int(key))
        for key in self.parameters: parameter_range.append(int(key))
        
        field_range.sort()
        dotfield_range.sort()
        parameter_range.sort()
        
        assert field_range == range(
            self.nF), "Invalid field indices for initial conditions: {} != {}".format(field_range, range(self.nF))
        assert dotfield_range == range(
            self.nF), "Invalid dotfield indices for initial conditions: {} != {}".format(dotfield_range, range(self.nF))
        assert parameter_range == range(
            self.nP), "Invalid parameter indices for model: {} != {}".format(parameter_range, range(self.nP))
            

class bispectrumCfgTemplate:
    
    def __init__(self):
        self.fNLConfigs = None
    
    
    def addBispectrumConfiguration(self, name, latex, alpha, beta):
        assert type(name) is str, "configuration name must be string"
        assert type(latex) is str, "latex definition must be string"
        assert type(alpha) is float, "alpha parameter must be float"
        assert type(beta) is float, "beta parameter must be float"
        
        if not hasattr(self, "fNLConfigs"): self.fNLConfigs = []
        
        self.fNLConfigs.append(
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
        self.minN = None
    
        self.rejectIfAny = None
        self.rejectIfAll = None
        
        self.acceptIfAny = None
        self.acceptIfAll = None
    
        self.groupNames = None
        
        self.recordMasses = True
    
    
    def setCoreParams(self, name, saveLocation="default", efoldsBeforeExit=55., subHorizonEvolution=6.,
                      integratorTols=np.array([1e-8, 1e-8]), adiabaticN=1., minN=60.):
        
        assert type(name) is str, "sampler name must be string to build directories"
        assert saveLocation == "default" or os.path.exists(saveLocation), "Path not found: {}".format(saveLocation)
        assert type(efoldsBeforeExit) in [float, int], "Specify N, s.t. Nend - N defines the pivot scale"
        assert type(subHorizonEvolution) in [float, int], "Specify the amount of subhorizon evolution"
        assert type(integratorTols) is np.ndarray, "Integrator tols must be numpy array"
        assert np.size(integratorTols) == 2, "Set abs. and rel. error (2 args) for tols."
        
        self.name = name
        self.saveLocation = saveLocation
        self.efoldsBeforeExit = efoldsBeforeExit
        self.subHorizonEvolution = subHorizonEvolution
        self.integratorTols = integratorTols
        self.adiabaticN = adiabaticN
        self.minN = minN
        
        assert efoldsBeforeExit < minN, "Horizon exit time greater than duration of inflation:" \
                                        "{} > {}".format(efoldsBeforeExit, minN)
    

    def addCriticalValueGroup(self, groupName, nature):
        
        assert type(groupName) is str, "groupName should be string type: {}".format(groupName)
        assert type(nature) is str and nature in [
            "acceptSample", "rejectSample"], "nature should belong to {}".format(str(["acceptSample",
                                                                                          "rejectSample"]))
    
        if self.groupNames is None: self.groupNames = []
        
        assert groupName not in self.groupNames, "groupName already use: {}".format(groupName)

        self.groupNames.append(groupName)

        if nature == "acceptSample":
    
            if self.acceptIfAll is None: self.acceptIfAll = {}
    
            self.acceptIfAll[groupName] = {'idx': np.array([], dtype=int), 'min': np.array([]), 'max': np.array([])}

        if nature == "rejectSample":
    
            if self.rejectIfAll is None: self.rejectIfAll = {}
    
            self.rejectIfAll[groupName] = {'idx': np.array([], dtype=int), 'min': np.array([]), 'max': np.array([])}
            
        else: raise KeyError, nature
        
    
    def __addCriticalValue__(self, fieldNumber, minValue, maxValue, attr_, nature=None, groupName=None):
        
        assert type(fieldNumber) is int and fieldNumber < self.nF and fieldNumber >= 0, fieldNumber
        assert type(minValue) in [int, float], minValue
        assert type(maxValue) in [int, float], maxValue
        
        if attr_ == "fields":
            idx = 1 + fieldNumber
        elif attr_ == "dotFields":
            idx = 1 + fieldNumber + self.nF
        else:
            raise KeyError, attr_
        
        if groupName is None:
            
            if nature == "acceptSample":

                if self.acceptIfAny is None:
                    self.acceptIfAny = {}
                    self.acceptIfAny['idx'] = np.array([idx], dtype=int)
                    self.acceptIfAny['min'] = np.array([minValue])
                    self.acceptIfAny['max'] = np.array([maxValue])
                else:
                    self.acceptIfAny['idx'] = np.concatenate((self.acceptIfAny['idx'], np.array([idx], dtype=int)))
                    self.acceptIfAny['min'] = np.concatenate((self.acceptIfAny['min'], np.array([minValue])))
                    self.acceptIfAny['max'] = np.concatenate((self.acceptIfAny['max'], np.array([maxValue])))
                
            elif nature == "rejectSample":

                if self.rejectIfAny is None:
                    self.rejectIfAny = {}
                    self.rejectIfAny['idx'] = np.array([idx], dtype=int)
                    self.rejectIfAny['min'] = np.array([minValue])
                    self.rejectIfAny['max'] = np.array([maxValue])
                else:
                    self.rejectIfAny['idx'] = np.concatenate((self.rejectIfAny['idx'], np.array([idx], dtype=int)))
                    self.rejectIfAny['min'] = np.concatenate((self.rejectIfAny['min'], np.array([minValue])))
                    self.rejectIfAny['max'] = np.concatenate((self.rejectIfAny['max'], np.array([maxValue])))
                    
            else: raise KeyError, nature
            
        else:
            
            assert groupName in self.groupNames, groupName

            if nature == "acceptSample":
                self.acceptIfAll[groupName]['idx'] = np.concatenate(
                    (self.acceptIfAll[groupName]['idx'], np.array([idx], dtype=int)))
                self.acceptIfAll[groupName]['min'] = np.concatenate(
                    (self.acceptIfAll[groupName]['min'], np.array([minValue])))
                self.acceptIfAll[groupName]['max'] = np.concatenate(
                    (self.acceptIfAll[groupName]['max'], np.array([maxValue])))

            elif nature == "rejectSample":
                self.rejectIfAll[groupName]['idx'] = np.concatenate(
                    (self.rejectIfAll[groupName]['idx'], np.array([idx], dtype=int)))
                self.rejectIfAll[groupName]['min'] = np.concatenate(
                    (self.rejectIfAll[groupName]['min'], np.array([minValue])))
                self.rejectIfAll[groupName]['max'] = np.concatenate(
                    (self.rejectIfAll[groupName]['max'], np.array([maxValue])))

            else: raise KeyError, nature


    def addAcceptSampleFieldValue(self, fieldNumber, minValue, maxValue, groupName=None):
        self.__addCriticalValue__(fieldNumber, minValue, maxValue, "fields", "acceptSample", groupName)
    
    
    def addRejectSampleFieldValue(self, fieldNumber, minValue, maxValue, groupName=None):
        self.__addCriticalValue__(fieldNumber, minValue, maxValue, "fields", "rejectSample", groupName)
    
    
    def addAcceptSampleDotFieldValue(self, fieldNumber, minValue, maxValue, groupName=None):
        self.__addCriticalValue__(fieldNumber, minValue, maxValue, "dotFields", "acceptSample", groupName)
    
    
    def addRejectSampleDotFieldValue(self, fieldNumber, minValue, maxValue, groupName=None):
        self.__addCriticalValue__(fieldNumber, minValue, maxValue, "dotFields", "rejectSample", groupName)
        

    def buildSampler(self, update=False):
        
        self.checkicps()
        if hasattr(self, "requiredModules"):
            self.requiredModules = list(set(self.requiredModules))
        else:
            self.requiredModules = []
    
        if self.saveLocation != "default":
            assert os.path.exists(self.saveLocation), "Savelocation does not exist: {}".format(self.saveLocation)
        else:
            self.saveLocation = os.path.abspath(os.path.join(os.getcwd(), "../sampler-builds"))
        
        # define core paths for directory structure
        self.pytpath = pytpath
        self.root = os.path.join(self.saveLocation, self.name)
        self.stats_dir = os.path.join(self.root, "stats")
        self.stats_dir_bg = os.path.join(self.stats_dir, "bg")
        self.stats_dir_2pf = os.path.join(self.stats_dir, "2pf")
        self.stats_dir_3pf = os.path.join(self.stats_dir, "3pf")
        self.samples_dir = os.path.join(self.root, "samples")
        self.outputs_dir = os.path.join(self.root, "outputs")
        self.localdata_dir = os.path.join(self.root, ".localdata")
        
        # Add additional locations of classes and methods
        self.classes_dir = os.path.abspath(os.path.join(os.getcwd(), "../classes"))
        self.methods_dir = os.path.abspath(os.path.join(os.getcwd(), "../sampler-methods"))
        
        # Build paths as required, subject to whether the model file is being updated
        for item in [self.root, self.stats_dir, self.stats_dir_bg,
                     self.stats_dir_2pf, self.stats_dir_3pf, self.samples_dir, self.outputs_dir, self.localdata_dir]:
            
            if update is False:
                assert not os.path.exists(item), "Directory already exists! {}".format(item)
                os.makedirs(item)
                
            else:
                if not os.path.exists(item):
                    os.makedirs(item)

        # Add numpy to required modules list s.t. we can build arrays for sample specimens
        if "numpy" not in self.requiredModules: self.requiredModules.append("numpy")

        # Setup import statements to write to file
        import_lines = ["import {}\n".format(m) for m in self.requiredModules]
        
        # Start function definition
        func_lines = ["\ndef genSample(n):\n"]
        
        # Start fields, dotfields and parameter arrays
        fields_lines = ["\tfvals = numpy.array([\n"]
        dotfields_lines = ["\tvvals = numpy.array([\n"]
        parameters_lines = ["\tpvals = numpy.array([\n"]
        
        # For every field index
        for ii in range(self.nF):
            
            # Load the execution method for the initial field and velocity value(s)
            f_command = self.fields[str(ii)]
            v_command = self.dotfields[str(ii)]
            
            f_command = str(f_command) if type(f_command) != str else f_command
            v_command = str(v_command) if type(v_command) != str else v_command
            
            # Add to lines list
            fval_line = "\t\t" + f_command
            vval_line = "\t\t" + v_command if v_command != "SR" else "\t\t" + "'" + v_command + "'"
        
            # Add end of line, delimiter by comma if not at end
            if ii < self.nF - 1:
                fval_line += ",\n"
                vval_line += ",\n"
            else:
                fval_line += "\n"
                vval_line += "\n"
            
            fields_lines.append(fval_line)
            dotfields_lines.append(vval_line)
        
        # Close initial condition definitions
        fields_lines.append("\t])\n\n")
        dotfields_lines.append("\t])\n\n")
        
        # Repeat process for parameters
        for ii in range(self.nP):
            
            p_command = self.parameters[str(ii)]
            
            p_command = str(p_command) if type(p_command) != str else p_command
            
            pval_line = "\t\t" + p_command
            
            if ii < self.nP - 1:
                pval_line += ",\n"
            else:
                pval_line += "\n"
            
            parameters_lines.append(pval_line)
        
        parameters_lines.append("\t])\n\n")
        
        # Combine lines to write function
        func_lines.append(fields_lines)
        func_lines.append(dotfields_lines)
        func_lines.append(parameters_lines)
        func_lines.append("\treturn n, fvals, vvals, pvals")
        
        # Define path to generator function
        genSamplePath = os.path.join(self.localdata_dir, "generator.py")
        
        # Build generator subject to update criteria
        if update is False:
            assert not os.path.exists(genSamplePath), "Generator already exists: {}".format(genSamplePath)
        else:
            if os.path.exists(genSamplePath): os.remove(genSamplePath)
        
        # Write to file
        f = open(genSamplePath, "w")
        with f:
            for line in import_lines + func_lines: f.writelines(line)
        
        # Define paths to copy sampler routine to local model directory
        runSamplerPath_keep = os.path.abspath(os.path.join(os.getcwd(), "../sampler-methods/run_sampler.py"))
        runSamplerPath_copy = os.path.join(self.root, "run_sampler.py")
        
        if update is False:
            assert not os.path.exists(runSamplerPath_copy), "Object already exists: {}".format(runSamplerPath_copy)
        else:
            if os.path.exists(runSamplerPath_copy):
                os.remove(runSamplerPath_copy)

        shutil.copyfile(runSamplerPath_keep, runSamplerPath_copy)
        
        # Define paths to (hidden) local data
        bispectraObjPath   = os.path.join(self.localdata_dir, "fNL.localdata")
        environmentObjPath = os.path.join(self.localdata_dir, "env.localdata")
        transObjPath = os.path.join(self.localdata_dir, "transport.localdata")
        latexObjPath = os.path.join(self.localdata_dir, "latex.localdata")
        
        # Define dictionary of key transport data
        transDict = {
            "transportModule" : self.transportModule.__name__,
            "exitN"           : self.efoldsBeforeExit,
            "subN"            : self.subHorizonEvolution,
            "intTols"         : self.integratorTols,
            "adiabaticN"      : self.adiabaticN,
            "minN"            : self.minN
        }
        
        # Define dictionary of key environment varaibles
        envDict = {
            'PyTS_pathPyT'       : self.pytpath,
            'PyTS_pathRoot'      : self.root,
            'PyTS_pathSamples'   : self.samples_dir,
            'PyTS_pathStats'     : self.stats_dir,
            'PyTS_pathLocalData' : self.localdata_dir,
            'PyTS_pathClasses'   : self.classes_dir,
            'PyTS_pathMethods'   : self.methods_dir
        }
        
        
        # Zip paths / objects and build. Add field space conditions if these are defined
        localPaths = [bispectraObjPath, environmentObjPath, transObjPath, latexObjPath]
        localDicts = [self.fNLConfigs, envDict, transDict, self.latex]
        
        if self.acceptIfAny is not None:
            localPaths.append(os.path.join(self.localdata_dir, "acceptIfAny.localdata"))
            localDicts.append(self.acceptIfAny)
            
        if self.rejectIfAny is not None:
            localPaths.append(os.path.join(self.localdata_dir, "rejectIfAny.localdata"))
            localDicts.append(self.rejectIfAny)
            
        if self.acceptIfAll is not None:
            localPaths.append(os.path.join(self.localdata_dir, "acceptIfAll.localdata"))
            localDicts.append(self.acceptIfAll)
            
        if self.rejectIfAll is not None:
            localPaths.append(os.path.join(self.localdata_dir, "rejectIfAll.localdata"))
            localDicts.append(self.rejectIfAll)
            
            
        localFiles = zip(localPaths, localDicts)
        
        for item in localFiles:
            lF, o = item
            
            if update is False:
                assert not os.path.exists(lF), "Local file object already exists: {}".format(lF)
            else:
                if os.path.exists(lF):
                    os.remove(lF)
                    
            f = open(lF, "wb")
            with f: pk.dump(o, f)
        