""" Example sampler setup file : Double-quadratic field inflation (non-canonical fieldspace conditions) """


# Ensure that the setupBuilds module is imported
from setupBuilds import *

# Import the _PyTransport module you want to perform sampling with
import PyTrans2Quad as PyT

# Initialize the model file using the _PyTransport module
model = PyTransportSampler(PyT)

# Set core params
model.setCoreParams("nc_exmaple")

# Set priors on initial conditions and parameters
model.setInitialFieldValues(-1, "numpy.random.uniform(5, 20)", "numpy")
model.setInitialFieldVelocities(-1, "SR")
model.setParameterValues(-1, "numpy.exp(numpy.random.uniform(-6, -1))", "numpy")

# Set configurations for the reduced-Bispectrum
model.addBispectrumConfiguration("eq", "$f_{NL}^{eq}", 1./3., 1./3.)

# # Once the parameters have been set, we simply execute buildSampler

model.recordFieldValue(0, "\chi_{0}")
model.recordDotFieldValue(0, "\dot\chi_{0}")
model.recordParameterValue(0, "M_\chi")

model.recordFieldValue(1, "\phi_{0}")
model.recordDotFieldValue(1, "\dot\phi_{0}")
model.recordParameterValue(1, "M_\phi")

""" We add non-canonical conditions to the fieldspace - in particular, if the 0th field travels below -10 Mpl, reject.
    We assume a hybrid-exit triggered by reaching 10 Mp l"""

model.addAcceptSampleFieldValue(0, -np.inf, 0)

model.addCriticalValueGroup("test", "rejectSample")

model.addRejectSampleFieldValue(0, -np.inf, -5, groupName="test")
model.addRejectSampleFieldValue(1, -np.inf, -5, groupName="test")

model.minN = 30

model.buildSampler(update=True)