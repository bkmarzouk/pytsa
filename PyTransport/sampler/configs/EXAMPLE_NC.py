""" Example sampler setup file : Double-quadratic field inflation (non-canonical fieldspace conditions) """


# Ensure that the setupBuilds module is imported
from . setup_sampler import SamplerMethods

# Import the _PyTransport module you want to perform sampling with
import PyTrans2Quad as PyT

# Initialize the model file using the _PyTransport module
sampler_setup = SamplerMethods(PyT)

# Set core params
sampler_setup.set("nc_exmaple")

# Set priors on initial conditions and parameters
sampler_setup.setInitialFieldValues(-1, "numpy.random.uniform(5, 20)", "numpy")
sampler_setup.setInitialFieldVelocities(-1, "SR")
sampler_setup.setParameterValues(-1, "numpy.exp(numpy.random.uniform(-6, -1))", "numpy")

# Set configurations for the reduced-Bispectrum
sampler_setup.addBispectrumConfiguration("eq", "$f_{NL}^{eq}", 1. / 3., 1. / 3.)

# # Once the parameters have been set, we simply execute buildSampler

sampler_setup.recordFieldValue(0, "\chi_{0}")
sampler_setup.recordDotFieldValue(0, "\dot\chi_{0}")
sampler_setup.recordParameterValue(0, "M_\chi")

sampler_setup.recordFieldValue(1, "\phi_{0}")
sampler_setup.recordDotFieldValue(1, "\dot\phi_{0}")
sampler_setup.recordParameterValue(1, "M_\phi")

""" We add non-canonical conditions to the fieldspace - in particular, if the 0th field travels below -10 Mpl, reject.
    We assume a hybrid-exit triggered by reaching 10 Mp l"""

sampler_setup.addAcceptSampleFieldValue(0, -np.inf, 0)

sampler_setup.addCriticalValueGroup("test", "rejectSample")

sampler_setup.addRejectSampleFieldValue(0, -np.inf, -5, groupName="test")
sampler_setup.addRejectSampleFieldValue(1, -np.inf, -5, groupName="test")

sampler_setup.minN = 30

sampler_setup.buildSampler(update=True)