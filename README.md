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

#### Assigning optional Sampler parameters

Generally speaking, we probably want to explore the space of parameters and initial conditions associated with the model. We now describe how to define sampling space on such data fields. We describe first the **basic**, then **advanced** configuration options.

**Basic Definitions**

Initial conditions and parameters are all defined with the same basic syntax

```python
model.setFieldValues(index, command, *modules)
model.setDotFieldValues(index, command, *modules)
model.setParameterValues(index, command, *modules)
```

The first argument, `index` corresponds to the index value of field *or* parameter associated with the PyTransport installation file (PyT). The second argument, `command` informs the sampler what to do. E.g. if we want to assign a constant value, simply enter the value. More interestingly is the case where we may wish to draw values from a distrribution. As an example, let's suppose that we want to draw values from the a uniform distribution for the field with index 0. This could be achieved quite easily using the ``numpy`` module in the following way

```python
model.setFieldValues(0, "numpy.random.uniform(-20, 20)", "numpy")
```

```python
model.set
