from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration


def configuration(parent_package='', top_path=''):
    return Configuration('', parent_package, top_path)


setup(
    name='pytsa',
    version='0.0.1',
    packages=['pytsa', 'pytsa.pyt', 'pytsa.cppt', 'pytsa.cppt.NC', 'pytsa.cppt.stepper',
              'pytsa.models', 'pytsa.sampler'],
    data_files=[],
    url='https://github.com/bkmarzouk/pytsa',
    author='Kareem Marzouk',
    author_email='bkmarzouk@gmail.com',
    description='This package is an extension to PyTransport, a tool for the compuation of inflationary correlation '
                'functions, which interfaces a random sampling technology to explore the correlations between model '
                'parameter spaces, initial conditions and observables.',
    install_requires=[],
    requires=['numpy', 'sympy', 'scipy', 'dill', 'mpi4py'],
    long_description="This is PyTransport Sampler (pytsa) - A package that builds upon existing the existing "
                     "PyTransport technology to sample inflationary model spaces.",
    configuration=configuration)
