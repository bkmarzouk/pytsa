#This file is part of _PyTransport.

#_PyTransport is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#_PyTransport is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with _PyTransport.  If not, see <http://www.gnu.org/licenses/>.


#module setup script

from distutils.core import setup, Extension
import numpy
import time
import os
dir = os.path.dirname(__file__)
filename = os.path.join(dir, 'PyTrans.cpp')
f = open(filename,"r")
lines = f.readlines()
f.close()
f = open(filename,"w")

for line in lines:
    if not line.startswith("// Package recompile"):
        f.write(line)

    if line.startswith("// Package recompile"):
        f.write('// Package recompile attempted at: '+ time.strftime("%c") +'\n')
f.close()        
filename2 = os.path.join(dir, '../CppTrans/stepper/rkf45.cpp')
dirs = os.path.join(dir, '../CppTrans/')

# don't edit the comment at the end of the setup line below #######################################
setup(name="PyTransdquad_euclidean", version="1.0", ext_modules=[Extension("PyTransdquad_euclidean", [filename, filename2 ])], include_dirs=[numpy.get_include(), dirs], extra_compile_args = ["-std=c++11 -frounding-math -fsignaling-nans -j 6"])#setup
###################################################################################################