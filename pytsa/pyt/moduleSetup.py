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
import numpy as np

#  Get path to sources we want to compile and directories that contain header files
cwd = os.path.dirname(__file__)
template_path = os.path.join(cwd, "PyTrans.cpp")
cppt_dir = os.path.abspath(os.path.join(cwd, "../cppt"))
stepper_path = os.path.join(cppt_dir, "stepper/rkf45.cpp")

assert os.path.exists(template_path), template_path
assert os.path.exists(cppt_dir), cppt_dir
assert os.path.exists(stepper_path), stepper_path

# Do not modify comment at end of following line
mod_name = 'dquad_2sphere'  # PYT_MODNAME

# module extension for c++ source contributions
module_extension = Extension(
    mod_name,
    sources=[template_path, stepper_path],
    language='c++',
    extra_compile_args=['-pipe'],
    extra_link_args=["-I{}".format(np.get_include())]
)

setup(
    name=mod_name,
    version="1.0",
    ext_modules=[module_extension],
    include_dirs=[np.get_include(), cppt_dir]
)
