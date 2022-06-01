# PyTransport Sampler (pytsa)
A general sampling tool for the prediction of observables in inflationary models, 
based on the code PyTransport. The code is designed fundamentally to compute the distribution(s) of observables for inflationary models, in particular, by employing random-sampling to circumvent the dimensionality of parameter and priori spaces.

To install the code as a package (recommended) run the following command from the root directory:
```bash
python -m install -e .
```

## Installing inflationary models
Installing inflationary models relies upon symbolic definitions of the potential and field-space metric if appropriate.
The core symbolic engine relies on ```SymPy```, and expressions and parameters should be defined using this package.

Skeleton from setting up models:
```python
from pytsa import pytrans_setup as setup
import sympy as sym

# Define number of fields and parameters for model
nf = <number of fields>
np = <number of params>

# Build symbolic array for fields, each indexed with f[<index>] like a list.
f = sym.symarray('f', nf)

# Build symbolic array for number of parameters, each indexed with p[<index>] like a list. 
p = sym.symarrat('p', np)

# Build potential from combination of fields and params.
pot = <some function of f, p>

# Define field space metric: Note index down covariant expression G_{ab}
met = sym.Matrix([[<metric components>], ..., [<more metric components>]])

# Translate model into c++ source code. Note that we pass the metric G=G
setup.potential(
    pot,  # symbolic definition of the potential
    nf,   # number of fields for the model
    nf,   # number of params for the model
    G=G,  # symbolic definition for the field space metric, optional
    simplify_fmet=True,  # simplifies expressions relating to field space metric
    simplify_pot=True,   # simplifies expressions relating to the potenital
    simplify_covd=True,  # simplifies covariant derivative expressions
    silent=False         # indicates whether to suppress output when building
)

# Run compiling step

model_name = <string corresponding to model name>
setup.compile_module(
    <model_name>,     # A string that defines the model name, this will be used to import
    <non-canonical>,  # True of False, indicating whether we have used a non-canonical field space metric
)
```

Once the script is defined, simply run ``python <installation_script.py>``, 
and an importable verions of the model can be imported within a script with 
```python
from pytsa.models import <model_name>
````

A selection of complete example installation files is provided in the ```./example_models``` directory,
which can be run as-is so long as ``pytsa`` has been installed.

## Building a sampling routine

Once a PyTransport model has been installed, we can build samplers to explore them.

In general, this requires writing a setup script which defines the priors on statistical distributions
that will be sampled from. Statistical distribution can be defined easily using the ``scipy.stats`` modules,
though custom distributions can be defined further so long as they retain OO methods.

Skeleton for constructing sampler routines:
```python
from pytsa.models import <model_name> as model
from pytsa.sampler import setup_sampler

from scipy.stats import <distributions>

# Initialize sampler object
name = <name_for_sampling_routine>
setup = setup_sampler.SamplerMethods(
    model,
    name,
    cache_loc=<location to build sampler>  # defaults to ./samplers
)

# Setup analysis hyper parameters, if not set will use defaults at build
setup.set_analysis_params(
    N_sub_evo=<efolds of subhorizon evolution>,                  # defaults to 6
    N_adiabatic=<last efolds to test adiabatic limit over>,      # defaults to 1
    N_min=<mimumum number of efolds of inflation>,               # defaults to 60
    tmax_bg=<max seconds to spend computing background>,         # defaults to 60
    tmax_2pf=<max seconds to spend computing 2point function>,   # defaults to 300
    tmax_3pf=<max seconds to spend computing 3point function>,   # defaults to 600
    tols=<list of abs. + rel. tols used by intergator >,         # defaults to [1e-5, 1e-5]
    step_density=<number of steps per efold used by integrator>  # defaults to 20
)

# Set priors on initial field values
setup.set_field(
    *<field indices>,                  # specify which fields the definition corresponds to
    method=<statistical distribution>  # define the prior, e.g. scipy.stats.uniform(-20, 40)
)

# Set priors on initial field velocty values
setup.set_dot_field(
    *<field_indices>,                  # same as above
    method=<statistical_distribution>  # same as above, except, can pass "sr" to use slow-roll equation
)

# Set priors on parameter values (analogous to field definitions)
setup.set_param(
    *<param indices>,
    method=<statistical_distribution>  
)

# build sampler
setup.build_sampler()
```

For the field and paramter prior definitions set in the above, note that we can call e.g. 
`setup.set_field` as many times as required to assign all definitions.

Note: If there is a field value / velocity or parameter that is to be set constant, simply pass the
constant value in place of `<statistical_distribution>` kwarg.

There are complete examples of the sampler setup procedure in ``./example_samplers``.