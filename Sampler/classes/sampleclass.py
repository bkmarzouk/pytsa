# Import core moduels for file and directory handling
import os
import pickle
import sys

# find tree for build directory (one level down)
tree = os.path.abspath(os.path.join(__file__, "../.."))
sys.path.append(tree)

# we define the sample class object now
from generator import gensample as gs #  testgen -> samplegenerator in final product

# We now define the sample template
class sample:

    # initialise sample with gensample function
    def __init__(self, mn, saveloc):

        # Generate values for sample and assign core attributes
        mn, fvals, vvals, pvals, latex = gs(mn)
        self.mn = mn
        self.fvals = fvals
        self.vvals = vvals
        self.pvals = pvals
        self.latex = latex
        self.saveloc = saveloc

        # Assign 'None' to attributes that will be determined by the end of inflation routin
        self.Nstart        = None
        self.Nend          = None
        self.NB            = None
        self.back          = None
        self.backExitMinus = None
        self.k             = None

    def Ndata(self):
        return [self.Nstart, self.Nend, self.NB]

    def kdata(self):
        return [self.back, self.backExitMinus, self.k]

    # define dictionary function to write initial ensemble information
    def dictionary(self):

        # Initialise column tuple with model identifier
        cols = ('model',)
        vals = (self.mn,)

        # Add arbitrary headers for initial field values, velocities and parameters & assign values
        for i in range(len(self.fvals)):
            cols+=('f{}'.format(i),)
            vals+=(self.fvals[i],)

        for i in range(len(self.vvals)):
            cols+=('v{}'.format(i),)
            vals+=(self.vvals[i],)

        for i in range(len(self.pvals)):
            vals+=(self.pvals[i],)
            cols+=('p{}'.format(i),)

        # return dictionary of data
        return dict(zip(cols, vals))

    # # we can then compute the background evolution of the sample, updating previous None types
    # def compute_evolution(self):
    #     #self.Nstart, self.Nend, self.NB, self.back, self.backExitMinus, self.k =
    #     oscripts.modelsetup(self)

    # finally we write a function to save the sample object
    def save(self):
        savename = os.path.join(self.saveloc, "m{}.pkl".format(self.mn))

        f = file(savename, 'wb')
        pickle.dump(self, f)
        f.close()
