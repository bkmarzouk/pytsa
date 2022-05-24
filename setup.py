from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration


def configuration(parent_package='', top_path=''):
    return Configuration('', parent_package, top_path)


setup(
    name='PyTransport Sampler',
    version='0.0.1',
    packages=['pytransport', 'pytransport.pyt', 'pytransport.cppt', 'pytransport.cppt.NC', 'pytransport.cppt.stepper',
              'pytransport.models', 'pytransport.sampler'],
    data_files=[],
    url='<todo>',
    author='Kareem Marzouk',
    author_email='<todo>>',
    description='<todo>',
    install_requires=[],
    requires=['numpy', 'sympy', 'scipy', 'mpi4py'],
    long_description="<todo>",
    configuration=configuration)
