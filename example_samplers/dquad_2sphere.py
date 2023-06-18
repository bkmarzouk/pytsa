from pytsa.models import dquad_2sphere
from pytsa.sampler import setup_sampler
from scipy.stats import uniform, chi2, expon

# Example sampler setup file for double quadratic model.

# Initialize sampler object with the model (imported at the tope)
# and a name (for identification)
setup = setup_sampler.SamplerMethods(dquad_2sphere, "dquad_2sphere_example")

# Setup analysis hyper parameters
setup.set_analysis_params(tols=[1e-5, 1e-5])

# Setup field initial conditions; here we sample on the domain [-20, 20]
setup.set_field(0, 1, method=uniform(-20, 40))

# Setup dot field initial conditions; here we just use slow roll values
setup.set_dot_field(0, 1, method="sr")

# Setup parameter sampling; here we sample each mass parameter from an
# exponential distribution
setup.set_param(0, 1, method=expon(scale=0.01))

# We also sample the metric parameter, R. Here we sample from the chi2
# distribution.
# We chose this dist. since it is positive definite and the metric depends
# on R^2, so we don't over sample.
setup.set_param(2, method=chi2(3))

# Finally, build the sampler
setup.build_sampler()
