import numpy as np
import os
import pickle as pk

rootpath = os.environ['PyTSamplerRoot']; import config as cfg

class realization:
    """ Initialise model realisaiton object:
        model number is used as an identifier for unifying observables
        rootdir """

    def __init__(self, modelnumber, parameters, background, adiabatic, Nend, kExit, savelocation):

        self.modelnumber = modelnumber # Track identifying model number
        self.parameters  = parameters  # Track which model parameters
        self.background  = background  # Track background evolution
        self.Nend        = Nend        # Record end of inflation
        self.kExit       = kExit       # Record scale that exits the horizon at Nexit

        self.observables = {}          # Will be populated with results

        self.adiabatic   = adiabatic   # Update with adiabatic limit from BG evolution

        self.savepath    = os.path.join(savelocation, "{}.sample".format(modelnumber))


        self.save_init()

    def update_observables(self, obsdict, name_timer):

        for key in obsdict:
            if key != 'time':
                assert key not in self.observables, "Realization already has observable data! {}".format(key)
                self.observables[key] = obsdict[key]
            else:
                time_key = "{}_t".format(name_timer)
                self.observables[time_key] = obsdict[key]

    def line_dict(self):
        line = {'weight': 1.0, 'like': 1.0}

        pinfo = cfg.parameter_values

        for p in pinfo:
            if p['LaTeX'] is not None:
                pn = p['ParameterNumber']    # Get parameter number of interest
                pv = self.parameters[pn]     # Get parameter value
                line['p{}'.format(pn)] = pv  # Add parameter to dictionary

        # Add all observable dictionary items to line
        for o in self.observables: line[o] = self.observables[o]

        return line

    def save_init(self):
        assert not os.path.exists(self.savepath), "Sample already exists: {}".format(self.savepath)
        f = open(self.savepath, "wb")
        with f: pk.dump(self, f)