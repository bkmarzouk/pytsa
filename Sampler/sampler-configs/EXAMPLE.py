""" Example sampler setup file : Double-quadratic field inflation """

# Ensure that the setupBuilds module is imported
from setupBuilds import *

# Import the PyTransport module you want to perform sampling with
import PyTrans2Quad as PyT

# Initialize the model file using the PyTransport module
model = PyTransportSampler(PyT)

# Set core params
model.setCoreParams("dquad_exmaple")

# Set priors on initial conditions and parameters
model.setInitialFieldValues(-1, "numpy.random.uniform(-20, 20)", "numpy")
model.setInitialFieldVelocities(-1, "SR")
model.setParameterValues(-1, "numpy.exp(numpy.random.uniform(-6, -1))", "numpy")

# Set configurations for the reduced-Bispectrum
model.addBispectrumConfiguration("eq", "$f_{NL}^{eq}", 1./3., 1./3.)

# Once the parameters have been set, we simply execute buildSampler
model.buildSampler(update=True)