""" Example sampler setup file : Double-quadratic field inflation """

# Ensure that the setupBuilds module is imported
from setupBuilds import *

# Import the _PyTransport module you want to perform sampling with
import PyTrans2Quad as PyT

# Initialize the model file using the _PyTransport module
model = PyTransportSampler(PyT)

# Set core params
model.setCoreParams("EXAMPLE_dquad")

# Set priors on initial conditions and parameters
model.setInitialFieldValues(-1, "numpy.random.uniform(-20, 20)", "numpy")
model.setInitialFieldVelocities(-1, "SR")
model.setParameterValues(-1, "numpy.exp(numpy.random.uniform(-6, -1))", "numpy")

# Set configurations for the reduced-Bispectrum
model.addBispectrumConfiguration("eq", "$f_{NL}^{eq}", None, None)
model.addBispectrumConfiguration("fo", "$f_{NL}^{fo}", None, None)
model.addBispectrumConfiguration("sq", "$f_{NL}^{sq}", None, None)

# Once the parameters have been set, we simply execute buildSampler
model.recordFieldValue(0, "\chi_{0}")
model.recordDotFieldValue(0, "\dot{\chi_{0}}")
model.recordParameterValue(0, "M_\chi")

model.recordFieldValue(1, "\phi_{0}")
model.recordDotFieldValue(1, "\dot{\phi_{0}}")
model.recordParameterValue(1, "M_\phi")

model.buildSampler()