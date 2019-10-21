import numpy as np
import os
import pickle as pk

class realization:
    """ Initialise model realisaiton object:
        model number is used as an identifier for unifying observables
        rootdir """

    def __init__(self, modelnumber, background, Nend, kExit, savelocation):

        self.modelnumber = modelnumber # Track identifying model number
        self.background  = background  # Track background evolution
        self.Nend        = Nend        # Record end of inflation
        self.kExit       = kExit       # Record scale that exits the horizon at Nexit

        self.observables = [] # Will be populated with results
        self.parameters  = [] # Will be populated with parameters

        self.adiabatic   = None # Update with adiabatic limit from BG evolution

        self.savepath    = os.path.join(savelocation, "{}.sample".format(modelnumber)) # Construct path for save

        self.save_init() # save sample object

    def Nend(self):
        return self.background.T[0][-1]

    def update_observables(self, obs_name, value, accept):
        assert type(obs_name) == str, "Must pass string definition for observable dictionary"

        observable = {}
        observable['name']   = obs_name
        observable['value']  = value
        observable['accept'] = accept
        assert observable not in self.observables, "Observable '{o}' already exists for sample: {m}".format(
            o = obs_name, m = self.modelnumber
        )
        self.observables.append(observable)

    def update_parameters(self, pname, pvalue):
        assert type(pname) == str, "Must pass string definition for parameter name"

        parameter = {}
        parameter['name'] = pname
        parameter['value'] = pvalue
        self.parameters.append(parameter)

    def line_dict(self):
        line = {}

        for item in self.observables:
            line[item['name']] = item['value']

        for item in self.parameters:
            line[item['name']] = item['value']

        return line

    def save_init(self):
        assert not os.path.exists(self.savepath), "Sample already exists: {}".format(self.savepath)
        f = open(self.savepath, "wb")
        with f: pk.dump(self, f)

