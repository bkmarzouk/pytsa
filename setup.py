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
    url='<todo>',
    author='Kareem Marzouk',
    author_email='<todo>>',
    description='<todo>',
    install_requires=[],
    requires=['numpy', 'sympy', 'scipy', 'mpi4py'],
    long_description="This is PyTransport Sampler (pytsa) - A package that builds upon existing the existing "
                     "PyTransport technology to sample inflationary model spaces.",
    configuration=configuration)
