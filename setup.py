from setuptools import setup, find_packages

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