# PyTransport Sampler (pytsa)
**A general sampling tool for the exploration of inflationary observables.**

This code is designed to derive general predictions for inflationary cosmology with random sampling techniques. 
The goal is not to optimize a model's compatibility with observation, but to make 
predictions from a robust exploration of initial conditions and parameter spaces.
A recent application of this software can be found in [this paper](https://arxiv.org/abs/2105.03637).

The code serves as an extension to 
[PyTransport](https://github.com/jronayne/PyTransport) (see accompanying [paper](https://arxiv.org/abs/1609.00381))
which was developed by David Mulryne and John Ronayne. Some parts of the current technology were developed with David Seery.

**Contents**
- Updated python 3 version of PyTransport, with internal routines for tensor algebra.
- A priori and Latin hypercube sampling methods.
- MPI ready: run sampling jobs locally or on computing clusters.

To install the code as a package run the following command from the root directory:
```bash
python -m pip install -e .
```

The setup should have installed an example model file `dquad_2sphere.py`. To test whether the submodule was installed,
you can run some unittests via
```bash
pytest tests/test_install.py
```

Demo files for single field & multiple field models can be found in the ```/notebooks``` directory. These files explain how to use the code: Installation, Implementing compiled modules, and sampling over models.
