import os
from setuptools import setup

setup(
    name='pytsa',
    version='0.1.0',
    author="Kareem Marzouk",
    author_email="bkmarzouk@gmail.com",
    packages=['pytsa'],
    description="This package contains the PyTransport source, interfaced with sampling technology.",
    url="https://github.com/bkmarzouk/pytsa",
    install_requires=[
        "scipy>=1.6.1",
        "numpy>=1.20.1",
        "dill>=0.3.3",
        "mpi4py>=3.0.3",
        "schwimmbad>=0.3.2"
    ],
    python_requires=">=3.8, <4"
)

sym_cache = os.path.join(os.path.dirname(__file__), "pytsa", "sym_cache")
if not os.path.exists(sym_cache):
    os.makedirs(sym_cache)

for sub_dir in ['model', 'fmet', 'pot', 'covd']:
    sub_dir = os.path.join(sym_cache, sub_dir)

    if not os.path.exists(sub_dir):
        os.makedirs(sub_dir)
