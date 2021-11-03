# This file is part of _PyTransport.

# _PyTransport is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# _PyTransport is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with _PyTransport.  If not, see <http://www.gnu.org/licenses/>.

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

# Do not modify comment at end of following line
mod_name = 'PyTransTEST'  # PYT_MODNAME

# module extension for c++ source contributions
module_extension = Extension(
    mod_name,
    sources=[template_path, stepper_path],
    language='c',
    extra_link_args=["-undefined", "-dynamic_lookup"]
)

setup(
    name=mod_name,
    version=1.0,
    ext_modules=[module_extension],
    include_dirs=[np.get_include(), cppt_dir, os.path.join(cppt_dir, "NC"), os.path.join(cppt_dir, "stepper")]
)
