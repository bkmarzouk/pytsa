# PyTransport Sampler (pytsa)
**A general sampling tool for the exploration of inflationary observables.**

This code is designed to derive general predictions for inflationary cosmology with random sampling techniques. 
The goal is not to optimize a model's compatibility with observation, but to make 
predictions from a robust exploration of initial conditions and parameter spaces.
A recent application of this software can be found in [this paper](https://arxiv.org/abs/2105.03637).

The code serves as an extension to 
[PyTransport](https://github.com/jronayne/PyTransport) (see accompanying [paper](https://arxiv.org/abs/1609.00381))
which was developed by David Mulryne and John Ronayne.

**Contents**
- Updated python 3 version of PyTransport, with internal routines for tensor algebra.
- A priori and Latin hypercube sampling methods.
- MPI ready: run sampling jobs locally or on computing clusters.

To install the code as a package run the following command from the root directory:
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
sampler_name = <name_for_sampling_routine>
setup = setup_sampler.SamplerMethods(
    model,
    sampler_name,
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

Once the setup script has been writting, simply run ``python <sampler_setup.py>`` and the 
core defintions and a `run.py` file will be constructed and placed in the cache location.
By default, this will be `./samplers/<sammpler_name>`.

If a misdefinition is detected or if the definitions are insufficient, an error will be thrwon at build time.
Moreover, if an example already exists with the defined name at the cache location
with differing parameters an error will again be raised.

### Running sampling routines

Once the core sampler has been built, we are ready to go. Executing the sampler
is performed via the command line, and this is where `mpi` definitions can further be assigned.

The command structure is as follows:
```bash
<mpi_args> python run.py <routine_name> <other_flags>
```
where `<mpi_args>` may be e.g. `mpiexec -n 32` for a 32-core process. The routine
name is used to help distinguish different execution which may differ 
in incompatible ways based on flags (described below).

**Required**
- `--n_samples <int>` number of samples to compute.
- `--apriori` or `--latin` runs the sampler in apriori or Latin-hypercube sampling mode

**Optional**
- `--ns` indicates computation of the 2pt correlation function: spectral index, running, and scalar amplitude.
- `--eq` indicates computation of equilateral fNL.
- `--fo` indicates computation of folded fNL.
- `--sq` indicates computation of moderately-squeezed fNL.
- `--alpha <*vals>` and `--beta <*vals>`, indicate additional fNL calculations. 
Values should correspond to Fergusson-Shellard parameterization.
- `--verbose` prints finer-grained information whilst sampling.
- `--entropy` define initial system entropy for random number generators. Assures reproducibility.

If only **required** information is given, then the sampler will only report results for
slow-roll parameters (at horizon exit and end of inflation) 
as well as the mass-matrix eigenvalues at the same time-steps.

#### Reading the output
If a sampling job has successfully complete there will be some files at your disposal.
These will be located in the `<sampler_name>/<routine_name>` directory:
- `errors_<routine_name>.txt` containing an overview of rejected samples and overall efficiency.
- `getdist_<routine_name>.txt` and `getdist_<routine_name>.paramnames` containing
samples with complete data. These can be read directly into `GetDist` analysis software.
- `pandas_<routine_name>.df`, a `pandas` data-frame containing all sample data, any failed components of the data 
collection will be represented with `numpy.nan` values. Appropriate post-processing in this file has been
left for the user. In order to load the data, simply call e.g. 
`data = pandas.read_pickle("/path/to/pandas_<routine_name>.df")`. To check what's there you can use `data.columns`,
then index whatever you want with `data['<column_name>]`.