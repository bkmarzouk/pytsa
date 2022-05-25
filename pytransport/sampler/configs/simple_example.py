from pytransport.sampler.setup_sampler import SamplerMethods

# Import PyTransSetup module to set paths to compiled cpp module files


# Now we can import the compiled module
import dquad_euclidean as PyT

# We will make use of scipy's statistical libraries to sample parameter space
import scipy.stats as stats

sampler_setup = SamplerMethods(PyT)

sampler_setup.set_field(0, 1, method=stats.uniform(-20, 20))
sampler_setup.set_dot_field(0, 1, method="sr")
sampler_setup.set_param(0, 1, method=stats.loguniform(1e-5, 1e-2))

sampler_setup.build_sampler()
