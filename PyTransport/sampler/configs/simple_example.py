from pytransport.sampler.configs.setup_sampler import Setup

# Import PyTransSetup module to set paths to compiled cpp module files
from pytransport import pytrans_setup


# Now we can import the compiled module
import dquad_euclidean as PyT

# We will make use of scipy's statistical libraries to sample parameter space
import scipy.stats as stats

sampler_setup = Setup(PyT)

sampler_setup.set_field(0, 1, method=stats.uniform(-20, 20))
sampler_setup.set_dot_field(0, 1, method="sr")
sampler_setup.set_param(0, 1, method=stats.loguniform(1e-5, 1e-2))

sampler_setup.build_sampler()
