# PyT-Sample
A general sampling tool for the prediction of observables in inflationary models, based on the code PyTransport. The code is designed fundamentally to compute the distribution(s) of observables for inflationary models, in particular, by employing random-sampling to circumvent the dimensionality of parameter/priori-space.

## How to use the code
We begin by following the standard process of building a PyTransport installation file. In particular, we define a potential, and fieldspace metric (optional), and execute the compile step. During construction, model parameters are optionally defined -- PyT-Sample will enable these parameters to be fixed, or sampled over. Examples of files can be found in the ``Install'' directory.

After the model has been installed, a sampler parameter file needs to be configured. These should be written in the directory ``sampler-configs'', where there exists an object-oriented module that will be used to construct the sampler, along with all of it's supporting files.

## Sampler construction
In any sampler configuration file, we need to import two key modules:

1.  `from setupBuilds import *`
2.  `import <PyTransportModule> as PyT`

After this, we need to generate a PyTransportSampler instance. This should be initialized with the PyTransport module previously imported. You are free to name the object anything, in this readme and in the supplied example we adopt use of "model".

3.  `model = PyTransportSampler(PyT)`
