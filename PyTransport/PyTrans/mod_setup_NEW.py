from setuptools import setup, Extension
import os
import time
import numpy as np

try:
    n_threads = os.environ['OMP_NUM_THREADS']
except KeyError:
    n_threads = 1

cwd = os.path.dirname(__file__)
template_path = os.path.join(cwd, "PyTrans.cpp")

assert os.path.exists(template_path), template_path

with open(template_path, "r") as f:
    lines = f.readlines()

with open(template_path, "w") as f:
    for line in lines:

        if not line.startswith("// Package recompile"):
            f.write(line)

        if line.startswith("// Package recompile"):
            f.write('// Package recompile attempted at: ' + time.strftime("%c") + '\n')

cppt_dir = os.path.abspath(os.path.join(cwd, "../CppTrans"))
stepper_path = os.path.join(cppt_dir, "stepper/rkf45.cpp")

assert os.path.exists(cppt_dir), cppt_dir
assert os.path.exists(stepper_path), stepper_path

module_extension = Extension(
    'PyTransDQuad',
    sources=[template_path, stepper_path],
    language='c'
)

setup(
    name="PyTransDQuad",
    version=1.0,
    ext_modules=[module_extension],
    include_dirs=[np.get_include(), cppt_dir]
)

