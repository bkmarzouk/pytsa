# find tree for build directory (one level down)
import os
import sys
tree = os.path.abspath(os.path.join(__file__, "../.."))
sys.path.append(tree)

# import parsers and ast for literal evaluation of dictionaries in configuration files
import ast
import parsers
m_parser = parsers.get_moduleparser()
b_parser = parsers.get_bispectraparser()

# Create list for holding latex definition
latex = []

# Find options in parameter values
p_options = m_parser.options("ParameterValues")

# for each item
for item in p_options:

    # Literally evaluate to obtain dictionaries for parameters
    p = ast.literal_eval(m_parser.get("ParameterValues", item))

    # If the id is for all of the fields, demand None for latex definition to avoid ambiguity
    if p['id'] == 'ALL':

        # Get number of parameters
        np = m_parser.getint("Parameters", "nparams")
        for i in range(np):
            latex.append(None)

        # break out of loop
        break

    # If there is a list definition for paramters
    elif type(p['id']) == list:

        # For each item in the range, append None, otherwise there will be ambiguity in the definition
        for i in range(p['id'][0], p['id'][1]):
            latex.append(None)

    # Otherwise, add the latex defintion
    else:
        latex.append(p['latex'])


# We will need to understand calculation types
flavs = ["ns", "alpha"]
for item in b_parser.options("Shapes"):
    flavs.append(item)

# Define results class, to concatenate observables
class result:

    # Initialise with a number of observables, *obs
    def __init__(self, obs):

        import sampleclass
        import observableclass
        import numpy as np

        # We setup rejection counters from 0
        rcounter = 0

        # Record model numbers from obs arga
        modelnums = []
        ttotal = 0

        # We will generate a dictionary object to write our getdist file with
        rdict = {}

        # Iterate through observables
        for item in obs:

            modelnums.append(item.sample.mn)

            # Check that observables are the corrwct type
            assert item.__class__ == observableclass.observable, "Results can only be generated" \
                                                                 "from combining 'observable' type objects."

            # Add model number to pool
            modelnums.append(item.sample.mn)

            # check that model numbers are consistent
            assert item.sample.mn == modelnums[0],"Model number must be shared between" \
                                                  "observables for complete result."

            # check flavours and update values
            rdict[item.flavour] = item.result
            if item.flavour != "alpha":
                rdict["{}t".format(item.flavour)] = item.time
                ttotal += item.time
            rcounter += int(item.tde)

            # if nan result detected, add to counter
            try:
                if np.isnan(item.result) or type(item.result) == type(None):
                    rcounter += 1
            except TypeError:
                rcoutner += 1

        rdict['model'] = modelnums[0]
        rdict['ttotal'] = ttotal

        if rcounter > 0:
            self.rdict = None
        else:
            self.rdict = rdict

        # for consistent model number, each obs has common sample origin
        self.sample = obs[0].sample

    # Generate dictionary for writing results
    def dictionary(self):

        # Call dictionary that has been constructed by combining observables
        dict = self.rdict

        # If combined results form a None type dictionary, immediately return None
        if dict is None:
            return None

        else:
            dict['weight'] = 1.0
            dict['like'] = 1.0

            # For each model parameter
            for i in range(len(self.sample.pvals)):

                # if there exists a latex definition, i.e. one we want to write
                if latex[i] != None:

                    # Add arbitrary parameter name & assign entry value
                    dict["p{}".format(i)] = self.sample.pvals[i]

            return dict