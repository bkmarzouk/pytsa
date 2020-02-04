# PyT-Sample
A general sampling tool for the prediction of observables in inflationary models, based on the code PyTransport. The code is designed fundamentally to compute the distribution(s) of observables for inflationary models, in particular, by employing random-sampling to circumvent the dimensionality of parameter and priori spaces.

## How to use the code
We begin by following the standard process of building a PyTransport installation file. In particular, we define a potential, and fieldspace metric (optional), and execute the compile step. During construction, model parameters are optionally defined -- PyT-Sample will enable these parameters to be fixed, or sampled over. Examples of files can be found in the ``Install'' directory.

After the model has been installed, a sampler parameter file needs to be configured. These should be written in the directory ``sampler-configs'', where there exists an object-oriented module that will be used to construct the sampler, along with all of it's supporting files.

## Initializing a Sampler Object
In any sampler configuration file, we must include the following import statements

```python
from setupBuilds import *
import <PyTransportModule> as PyT
```

Following this, a PyTransportSampler instance can be created. This is initialized with the PyTransport module, i.e. the previously imported `PyT`:

```python
model = PyTransportSampler(PyT)
```

after which, the sampler can be configured as we wish.

#### Assigning core Sampler properties

The core attributes for the sampler should now be assigned in the following way:

```python
model.setCoreParams(<sampler name>, **kwargs)
```

The only required field is the sampler name, this is integral, as it is required to build the sampler's directory structure. There are a number of other keyword args that can be assigned at this point if you wish:
- `savelocation, dtype=str` : This defines the location where the sampler (and hence it's associated output files will be written). By default, all samplers are constructed in the directory *sampler-builds*. You may wish to change this if for example, you are operating on a remote-server and want to write outputs to scratch.
- `efoldsBeforeExit, dtype=float or int` : This defines the number of eFolds before the end of inflation we wish to track modes that exit the horizon. I.e. this defines the pivot-time for compuations. By default, this is *55*.
- `subHorizonEvolution, dtype=float or int` : This define the number of efolds of sub-horizon evolution perturbations should be evolved from. Loosely speaking, the more sub-horizon evolution that is tracked, the more accurate the calculation is - the expense of course being compute time. We refer the reader to the PyTransport docs. for further information. By default, this value is *6*.
- `integratorTols, dtype=np.ndarray` : This defines the absolute, and relative tols. for the integration stepper. By default, the value is *np.array([1e-8, 1e-])*, which usually performs well.
- `adiabaticN, dtype=float or int` : This defines the number of eFolds for which we test an adiabatic limit (see below for further information). By default, the value is *1*, such that we test adiabaicity over the final eFold of inflation.
- `minN, dtype=float or int` : This defines the minimum number of eFolds inflation must last for. By defualt this is *60*.

In principle, after the core parameters have been defined, the sampler can be constructed, though additional configuration will yield much more interesting data. The sampler is constructed via

```python
model.buildSampler(**kwarg)
```

A directory `/<save location>/<sampler name>/` will now be constructed, containing a `run_sampler.py` file. Note that the constructor attempts to cross-check variable assignments, such that potential errors will hopefully be raised prior to runtime. There is a single key word argument `update` which by default is `False`. If this set to `True`, any changes to the sampler parameters will be updated.

### Configuring the Sampler

Generally speaking, we probably want to explore the space of parameters and initial conditions associated with the model. We now describe how to define sampling space on such data fields. We describe first the **basic**, then **advanced** configuration options.

#### Basic: Defining intial conditions and model parameters

Initial conditions and parameters are all defined with the same basic syntax

```python
model.setFieldValues(index, command, *modules)
model.setDotFieldValues(index, command, *modules)
model.setParameterValues(index, command, *modules)
```

The first argument, `index` corresponds to the index value of field *or* parameter associated with the PyTransport installation file (PyT). The second argument, `command` informs the sampler what to do. Let's say that we wish the zeroth parameter to always have a fixed value of 0.1, we would define this in the following way:

```python
model.setParameterValues(0, 0.01)
```

More interestingly is the case where we may wish to draw values from a distrribution. To illustrate. let's suppose that we want to draw values from the a uniform distribution for the field with index 0. This could be achieved quite easily using the ``numpy`` module in the following way:

```python
model.setFieldValues(0, "numpy.random.uniform(-20, 20)", "numpy")
```

Often we may want to prescribe repeated value definitions. Rather than writing similar statements over-and-over, in place of a single index value we can pass a list (or any iterable) containing multiple indices, or `-1`. In the latter case, all available indices will be assigned the same command.

#### Basic: Recording initial conditions and parameter values

By default, PyT-Sample will take the relevant definitions and perform calculations: ultimately saving outputs to a result file. In the case of sampling, we might want to keep track of some additional information. To do this, the following optional commands are available:

```python
model.recordFieldValue(index, LaTeX)
model.recordDotFieldValue(index, LaTeX)
model.recordParameterValue(index, LaTeX)
```

The `index` argument is analagous to that previously discussed. Notabaly however, it does not currently permit multiple-value assignments (this will be updated soon). The `LaTeX` field is a string that can be formatted in *LaTeX* style. This is useful for example if we wish to export our results to `GetDist` for analysis (see later discussion).

#### Basic: Defining reduced-Bispectrum configurations

PyTransport supports 3-point function calculations, which are wrapped into `fNL` calculations by PyT-Sample. During sampler configuration, it is possible to define which fNL calculations will be available at runtime. It is possible to add bispectrum configurations with the following commands

```python
model.addBispectrumConfiguration(name, latex, alpha, beta)
```

The `name` field is simply used to distinguish between configurations. `latex` enables *LaTeX* formatting in analysis files. Finally, `alpha` and `beta` are simply used to parameterize the relative momenta about the pivot-scale:

```python
k1 = kExit / 2. - beta * kExit / 2.
k2 = kExit * (1. + alpha + beta) / 4.
k3 = kExit * (1. - alpha + beta) / 4.
```

For example, we may wish to define the following canonical configurations

```python
model.addBispectrumConfiguration("equilateral", "f_{NL}^{eq}", 1./3., 1./3.)
model.addBispectrumConfiguration("squeezed", "f_{NL}^{sq}", 0.9, 0.01)
model.addBispectrumConfiguration("folded", "f_{NL}^{fo}", -0.5, 0.5)
```

#### Advanced: Non-trivial sample acceptance and rejection

Unless otherwise specified, samples will be accepted if they satisfy the basic condition of supporting `minN` eFolds before the end of inflation, which is canonically defined to be when *\epsilon >= 1*. Some models of inflation may end or be violated if the fieldspace attains some critical value. A motivating example being D-Brane inflation, which ends with brane collision, and is violated in the case of brane ejection.

To add definitions of this kind, the following methods are available

```python
model.addAcceptSampleFieldValue(index, minValue, maxValue, groupName=None)
model.addRejectSampleFieldValue(index, minValue, maxValue, groupName=None)
model.addAcceptSampleDotFieldValue(index, minValue, maxValue, groupName=None)
model.addRejectSampleDotFieldValue(index, minValue, maxValue, groupName=None)
```

which essentially define regions of field/velocity-space that the background data is compared against. `index` is the field index, `minValue` is minimum value of the region, `maxValue` is maximum value of the region and `groupName` enables conditions to be coupled (see below description).

As an example, let's suppose that inflation can only end if the zeroth field value reaches a value of 0.02, and is rejected if it reaches 1.0. This would be implemented as follows:

```python
model.addAcceptSampleFieldValue(0, -np.inf, 0.02)
model.addRejectSampleFieldValue(0, 1.0, np.inf)
```

note that we make use of the `numpy` definition of infinity to dictate that the fieldspace bounds are half-open domains. These definitions mean that the sampler performs the following sequence of tasks: 1. Compute the background equations of motion. 2. If `f_0 >= 1` at any time step, reject the sample 3. If `f_0` <= 0.02 and N >= minN, accept the sample.

#### Advanced: Coupling sampling criteria

In the above example, there were single criteria that determined whether a background trajectory would be accepted or rejected. We now describe how we define simultaneous conditions. To illustrate, let's suppose we wish to permit the zeroth and first field to attain negative fieldspace values, but neither at the same time. Firstly, we would have to define a critical value group. The general syntax for doing so is

```python
model.addCriticalValueGroup(groupName, nature)
```

The `groupName` argument simply serves as an identifier, which can be called whatever you like. The "nature" argument indicates whether the group will imply sample rejection, or acceptence. Hence nature must be passed as "rejectSample" or "acceptSample". Hence we could initialize a group of conditions for our example via

```python
model.addCriticalGroupValue("negative", "rejectSample")
```

which would enable us to then link the fieldspace constraints described

```python
model.addRejectSampleFieldValue(0, -np.inf, 0.0, groupName="negative")
model.addRejectSampleFieldValue(1, -np.inf, 0.0, groupName="negative")
```
