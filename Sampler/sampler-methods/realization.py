import numpy as np
import os
import pickle as pk

rootpath = os.environ['PyTS_pathRoot']
localPath = os.environ['PyTS_pathLocalData']

class realization:

    def __init__(self, modelnumber, fields, velocities, parameters,
                 background, adiabatic, Nend, kExit, savelocation, overwrite=False):
        """"""

        """ Record core model data for reproducibility """
        self.modelnumber = modelnumber # Model number
        self.fields      = fields      # Field positions
        self.velocities  = velocities  # Field velocities
        self.parameters  = parameters  # Model parameters
        
        """ Record background data which is recycled in observable tasks"""
        self.background  = background  # Background evolution
        self.Nend        = Nend        # Record end of inflation
        self.kExit       = kExit       # record momenta at horizon exit

        """ Record results and classification data"""
        self.observables = {}          # Will be populated with results
        self.adiabatic   = adiabatic   # Adiabitic limit found
        self.savepath    = os.path.join(savelocation, "{}.sample".format(modelnumber))
    
        """ Save sample upon initialization"""
        self.save_init(overwrite=overwrite)

    def update_observables(self, obsdict, name_timer):
        """"""
        
        """ We update the model's observable dictionary """
        for key in obsdict:
            if key != 'time':
                assert key not in self.observables, "Realization already has observable data! {}".format(key)
                self.observables[key] = obsdict[key]
            
            else:
                time_key = "{}_t".format(name_timer)
                self.observables[time_key] = obsdict[key]

    def line_dict(self):
        """"""
        
        """ TODO: Add parameter localmodel file such that we can recover latex data """
        
        """ We initialize a dictionary for writing lines to results file"""
        line = {'weight': 1.0, 'like': 1.0}
        
        """ Load definitions of sample data we want to record in results """
        latexFile = open(os.path.join(localPath, "latex.localdata"), "rb")
        
        with latexFile as f: lDefs = pk.load(f)
        pNums = []
        pVals = []
        pDefs = []
        
        fNums = []
        fVals = []
        fDefs = []
        
        vNums = []
        vVals = []
        vDefs = []
        
        for key in lDefs:
            x, num = key[0], int(key[1:])
            if x=="p":
                line[key] = self.parameters[num]
                # pNums.append(num)
                # pVals.append(self.parameters[num])
                # pDefs.append(lDefs[key])
            elif x=="f":
                line[key] = self.fields[num]
                # fNums.append(num)
                # fVals.append(self.fields[num])
                # fDefs.append(lDefs[key])
            elif x =="v":
                line[key] = self.velocities[num]
                # vNums.append(num)
                # vVals.append(self.velocities[num])
                # vDefs.append(lDefs[key])
            else:
                raise KeyError, "Unrecognized key: {}".format(key)

        """ Add results for observables """
        for o in self.observables:
            line[o] = self.observables[o]

        return line

    def save_init(self, overwrite=False):
        """"""
        
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