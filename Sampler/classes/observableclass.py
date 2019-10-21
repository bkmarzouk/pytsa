# Import os and sys to navigate one level lower in directories: add to system path
import os
import sys
tree = os.path.abspath(os.path.join(__file__, "../.."))
sys.path.append(tree)

# import parsers and get
import parsers
b_parser = parsers.get_bispectraparser()

# Define observables class
class observable:

    # Initialise obj with its original Sample data, flavour, tiemdependance, result and compute time
    def __init__(self, Sample, Flavour, TimeDependent, Result, Time):

        import numpy as np
        import sampleclass

        # Define possible flavours for observables
        Flavours = ["ns", "alpha"]
        shapes = b_parser.options("Shapes")
        for item in shapes:
            Flavours.append(item)

        # Assert correct data types for sample
        assert Sample.__class__ == sampleclass.sample, "Observable must be constructed with object 'sample'."

        # Assert that the flavour is understandable
        assert Flavour in Flavours, "'Flavour' of calculation unrecognised."

        # Assert that the time dependence is a bool type
        assert type(TimeDependent) == bool, "'TimeDependent' must be bool type."

        # Ensure standard float type, rather than numpy float
        if type(Result) == np.float64:
            Result = float(Result)

        # We allow passage of None type results as they correspond to rejected samples and are written accordingly.
        # Hence assert results are float or None type
        assert type(Result) == float or type(Result) == type(None), "'Result' must be float type or 'None'."
        assert type(Time) == float, "'Time' must be float type."


        # If the observable satisfies all this criteria, it is a datatype we can understand. Make attributes of observable callable.
        self.sample = Sample
        self.flavour = Flavour
        self.tde = TimeDependent
        self.result = Result
        self.time = Time

    # Function to save class object via pickle
    def save(self, savedir):

        # Import pickle module
        import pickle
        import os
        #import sampleclass #

        # Construct savename from modelnumber and some input savedir
        savename = os.path.join( savedir, "m{}.pkl".format(self.sample.mn) )

        # File handler
        f = file(savename, "wb")

        # Dump with pickle
        pickle.dump(self, f)

        # close file
        f.close()
