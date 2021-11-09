import os
import pickle as pk

import numpy as np

rootpath = os.environ['PyTS_pathRoot']
localPath = os.environ['PyTS_pathLocalData']


class realization:
    
    def __init__(self, modelnumber, fields, velocities, parameters,
                 background, adiabatic, Nend, kExit, savelocation, overwrite=False):
        """ Object tracks sample data; used as a ledger for ICs / params / results """
        
        # Core data (can rerun sample with this)
        self.modelnumber = modelnumber  # Model number
        self.fields = fields  # Field positions
        self.velocities = velocities  # Field velocities
        self.parameters = parameters  # Model parameters
        
        # Background data (prevents recomputing)
        self.background = background  # Background evolution
        self.Nend = Nend  # Record end of inflation
        self.kExit = kExit  # record momenta at horizon exit
        
        # Results / write info.
        self.observables = {}  # Will be populated with results
        self.adiabatic = adiabatic  # Adiabitic limit found
        self.savepath = os.path.join(savelocation, "{}.sample".format(modelnumber))
        
        self.reject = False  # default is keep
        
        # Save object
        self.save_init(overwrite=overwrite)
    
    def check_reject(self):
        
        reject = False
        
        for key in self.observables:
            if self.observables[key] is None:
                reject = True
                
        self.reject = reject
        
        self.save_init(overwrite=True)
    
    def update_observables(self, obsDict):
        
        for key in obsDict:
            
            # TODO: Possible write-safe issue here (comment out)... Though testing suggests there shouldn't be an issue.
            # Turns out it's now simpler to allow rewrites, now we have mitigate directory structure defs. for obs.
            
            # assert key not in self.observables, "Attempting to rewrite observable data! key={}".format(key)
            
            self.observables[key] = obsDict[key]
            
            if obsDict[key] is None: self.reject = True
        
        self.save_init(overwrite=True)
    
    
    def line_dict(self):
        """"""
        
        """ We initialize a dictionary for writing lines to results file"""
        line = {'weight': 1.0, 'like': 1.0}
        
        """ Load definitions of sample data we want to record in results """
        latexFile = open(os.path.join(localPath, "latex.localdata"), "rb")
        
        with latexFile as f:
            lDefs = pk.load(f)
        
        for key in lDefs:
            x, num = key[0], int(key[1:])
            if x == "p":
                line[key] = self.parameters[num]
            elif x == "f":
                line[key] = self.fields[num]
            elif x == "v":
                line[key] = self.velocities[num]
            else:
                raise KeyError, "Unrecognized key: {}".format(key)
        
        """ Add results for observables """
        for o in self.observables: line[o] = self.observables[o]
        
        return line
    
    
    def save_init(self, overwrite=False):
        """ Saves / initializes binary file for sample object. Deletes existing object if overwrite is True """
        
        """ If overwrite is True, we remove the existing model data and rewrite the file """
        if overwrite is True:
            assert os.path.exists(self.savepath), "Sample does not exist: {}".format(self.savepath)
            os.remove(self.savepath)
            f = open(self.savepath, "wb")
            with f:
                pk.dump(self, f)
        
        else:
            """ Avoid dumping to same binary file"""
            assert not os.path.exists(self.savepath), "Sample already exists: {}".format(self.savepath)
            f = open(self.savepath, "wb")
            with f:
                pk.dump(self, f)
    
    
    def load_sample(self):
        """"""
        
        """ We load the binary file corresponding to the model realization """
        f = open(self.savepath, "rb")
        with f:
            sample = pk.load(f)
        
        return sample
