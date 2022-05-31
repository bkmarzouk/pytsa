from pytsa.models import quad
from pytsa.sampler import setup_sampler
from scipy.stats import uniform, expon

# Example sampler setup file for double quadratic model.

# Initialize sampler object with the model (imported at the tope) and a name (for identification)
setup = setup_sampler.SamplerMethods(quad, "quad_example")

# Setup analysis hyper parameters
setup.set_analysis_params(tols=[1e-5, 1e-5])

# Setup field initial conditions; here we sample on the domain [-20, 20]
setup.set_field(0, method=uniform(-20, 40))

# Setup dot field initial conditions; here we just use slow roll values
setup.set_dot_field(0, method="sr")

# Setup parameter sampling; here we sample each mass parameter from an exponential distribution
setup.set_param(0, method=expon(scale=0.01))

# Finally, build the sampler
setup.build_sampler()
