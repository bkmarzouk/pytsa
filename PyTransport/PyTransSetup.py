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


# This file contains python scripts used to setup the compiled PyTrans module

import sympy as sym
import sys
import site
import subprocess
import platform
import os
import shutil
import time as t
import pickle as pk
from PyTransport import gravtools_pyt


def directory(NC: bool):
    """
    Configures directory structure for PyTransport

    :param NC: if True methods use non-canonical field metric
    :return: None, configure PyTrans.cpp file
    """
    dir = os.path.dirname(__file__)

    filename = os.path.join(dir, 'PyTrans', 'PyTrans.cpp')

    with open(filename, "r") as f:
        lines = f.readlines()

    with open(filename, "w") as f:

        # Create directory to hold records of curvature objects;
        # We will compute inverse metrics in sympy notation which then can
        # quickly generate christoffel symbols and riemann tensors relevant for mass matrix computations
        curv_dir = os.path.join(dir, 'PyTrans', 'CurvatureRecords')

        if not os.path.exists(curv_dir):
            print("Building curvature directory")
            os.makedirs(curv_dir)

        if NC is False:
            for line in lines:
                if not line.endswith("//evolve\n") and not line.endswith("//moments\n") and not line.endswith(
                        "//model\n") and not line.endswith("//stepper\n"):
                    f.write(line)
                if line.endswith("//evolve\n"):
                    fileT = os.path.join(dir, 'CppTrans', 'evolve.h')
                    f.write('#include' + '"' + fileT + '"' + '//evolve' + '\n')
                if line.endswith("//moments\n"):
                    fileT = os.path.join(dir, 'CppTrans', 'moments.h')
                    f.write('#include' + '"' + fileT + '"' + '//moments' + '\n')
                if line.endswith("//model\n"):
                    fileT = os.path.join(dir, 'CppTrans', 'model.h')
                    f.write('#include' + '"' + fileT + '"' + '//model' + '\n')
                if line.endswith("//stepper\n"):
                    fileT = os.path.join(dir, 'CppTrans', 'stepper', 'rkf45.hpp')
                    f.write('#include' + '"' + fileT + '"' + '//stepper' + '\n')
        else:
            for line in lines:
                if not line.endswith("//evolve\n") and not line.endswith("//moments\n") and not line.endswith(
                        "//model\n") and not line.endswith("//stepper\n"):
                    f.write(line)
                if line.endswith("//evolve\n"):
                    fileT = os.path.join(dir, 'CppTrans', 'NC', 'evolve.h')
                    f.write('#include' + '"' + fileT + '"' + '//evolve' + '\n')
                if line.endswith("//moments\n"):
                    fileT = os.path.join(dir, 'CppTrans', 'NC', 'moments.h')
                    f.write('#include' + '"' + fileT + '"' + '//moments' + '\n')
                if line.endswith("//model\n"):
                    fileT = os.path.join(dir, 'CppTrans', 'NC', 'model.h')
                    f.write('#include' + '"' + fileT + '"' + '//model' + '\n')
                if line.endswith("//stepper\n"):
                    fileT = os.path.join(dir, 'CppTrans', 'stepper', 'rkf45.hpp')
                    f.write('#include' + '"' + fileT + '"' + '//stepper' + '\n')


def pathSet():
    """
    Sets path structure for internal navigation

    :return: None
    """

    dir = os.path.dirname(__file__)
    site.addsitedir(dir)

    p = platform.system()
    if p == 'Windows':
        site.addsitedir(os.path.join(dir, 'PyTrans', 'Python', 'site-packages'))
        site.addsitedir(os.path.join(dir, 'PyTrans', 'Python', 'Lib', 'site-packages'))
        site.addsitedir(os.path.join(dir, 'PyTrans', 'Python' + sys.version[:3].translate(None, '.'), 'site-packages'))
        site.addsitedir(
            os.path.join(dir, 'PyTrans', 'Python' + sys.version[:3].translate(None, '.'), 'Lib', 'site-packages'))
    else:
        site.addsitedir(os.path.join(dir, 'PyTrans', 'lib', 'python', 'site-packages'))
        site.addsitedir(os.path.join(dir, 'PyTrans', 'lib', 'python' + sys.version[:3], 'site-packages'))
        site.addsitedir(os.path.join(dir, 'PyTrans', 'lib', 'site-python'))

    site.addsitedir(os.path.join(dir, 'PyTransScripts'))


def compileName(name, NC=False):
    directory(NC)
    dir = os.path.dirname(__file__)
    location = os.path.join(dir, 'PyTrans')
    filename1 = os.path.join(dir, 'PyTrans', 'moduleSetup.py')

    if NC is True:
        # There should exist a temporary file fom the calling 'potential'
        curv_path = os.path.join(dir, 'PyTrans', 'CurvatureRecords', 'TEMP.curvature')
        assert os.path.exists(curv_path), "Temporary curvature file not found! {}".format(curv_path)

        # Rename temporary curvature file with compilation name!
        new_curv_path = curv_path.replace('TEMP', name)
        os.rename(curv_path, new_curv_path)

    f = open(filename1, "r")
    lines = f.readlines()
    f.close()
    f = open(filename1, "w")
    for line in lines:
        if not line.endswith("#setup\n"):
            f.write(line)
        if line.endswith("#setup\n"):
            f.write(
                # 'setup(name="PyTrans' + name + '", version="1.0", ext_modules=[Extension("PyTrans' + name + '", [filename, filename2 ])], include_dirs=[numpy.get_include(), dirs], extra_compile_args = ["-std=c++11"])#setup\n')
                'setup(name="PyTrans' + name + '", version="1.0", ext_modules=[Extension("PyTrans' + name + '", [filename, filename2 ])], include_dirs=[numpy.get_include(), dirs], extra_compile_args = ["-std=c++11 -frounding-math -fsignaling-nans"])#setup\n')
    f.close()

    filename = os.path.join(dir, 'PyTrans', 'PyTrans.cpp')
    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    f = open(filename, "w")
    for line in lines:
        if not line.endswith("//FuncDef\n") and not line.endswith("//initFunc\n") and not line.endswith("//modDef\n"):
            f.write(line)
        if line.endswith("//FuncDef\n"):
            f.write(
                'static PyMethodDef PyTrans' + name + '_funcs[] = {{"H", (PyCFunction)MT_H,    METH_VARARGS, PyTrans_docs},{"nF", (PyCFunction)MT_fieldNumber,        METH_VARARGS, PyTrans_docs},{"nP", (PyCFunction)MT_paramNumber,        METH_VARARGS, PyTrans_docs},{"V", (PyCFunction)MT_V,            METH_VARARGS, PyTrans_docs},{"dV", (PyCFunction)MT_dV,                METH_VARARGS, PyTrans_docs},  {"ddV", (PyCFunction)MT_ddV,                METH_VARARGS, PyTrans_docs}, {"dotfieldsSR", (PyCFunction)MT_dotfieldsSR,        METH_VARARGS, PyTrans_docs}, {"massMatrix", (PyCFunction)MT_massMatrix,        METH_VARARGS, PyTrans_docs}, {"findEndOfInflation", (PyCFunction)MT_findEndOfInflation,        METH_VARARGS, PyTrans_docs}, {"backEvolve", (PyCFunction)MT_backEvolve,        METH_VARARGS, PyTrans_docs},    {"sigEvolve", (PyCFunction)MT_sigEvolve,        METH_VARARGS, PyTrans_docs},    {"alphaEvolve", (PyCFunction)MT_alphaEvolve,        METH_VARARGS, PyTrans_docs},    {NULL}};//FuncDef\n')
        if line.endswith("//modDef\n"):
            f.write('     //modDef\n')
        if line.endswith("//initFunc\n"):
            f.write(
                'void initPyTrans' + name + '(void)    {        Py_InitModule3("PyTrans' + name + '", PyTrans' + name + '_funcs,                       "Extension module for inflationary statistics");        import_array();   }//initFunc\n')

    f.close()

    my_env = os.environ.copy()
    my_env["PYTHONUSERBASE"] = location
    p = subprocess.Popen(["python", filename1, "install", "--user"], cwd=location, stdin=subprocess.PIPE,
                         stderr=subprocess.PIPE, env=my_env)

    stdout, stderr = p.communicate()

    pathSet()

    shutil.rmtree(os.path.join(location, 'build'), ignore_errors=True)


def compileName3(name, NC=False):
    directory(NC)
    dir = os.path.dirname(__file__)
    location = os.path.join(dir, 'PyTrans')
    filename1 = os.path.join(dir, 'PyTrans', 'moduleSetup.py')

    if NC is True:
        # There should exist a temporary file fom the calling 'potential'
        curv_path = os.path.join(dir, 'PyTrans', 'CurvatureRecords', 'TEMP.curvature')
        assert os.path.exists(curv_path), "Temporary curvature file not found! {}".format(curv_path)

        # Rename temporary curvature file with compilation name!
        new_curv_path = curv_path.replace('TEMP', name)
        os.rename(curv_path, new_curv_path)

    f = open(filename1, "r")
    lines = f.readlines()
    f.close()

    f = open(filename1, "w")
    for line in lines:
        if not line.endswith("#setup\n"):
            f.write(line)
        if line.endswith("#setup\n"):
            f.write(
                # 'setup(name="PyTrans' + name + '", version="1.0", ext_modules=[Extension("PyTrans' + name + '", [filename, filename2 ])], include_dirs=[numpy.get_include(), dirs], extra_compile_args = ["-std=c++11"])#setup\n')
                'setup(name="PyTrans' + name + '", version="1.0", ext_modules=[Extension("PyTrans' + name + '", [filename, filename2 ])], include_dirs=[numpy.get_include(), dirs], extra_compile_args = ["-std=c++11 -frounding-math -fsignaling-nans"])#setup\n')
    f.close()

    filename = os.path.join(dir, 'PyTrans', 'PyTrans.cpp')
    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    f = open(filename, "w")
    for line in lines:
        if not line.endswith("//FuncDef\n") and not line.endswith("//initFunc\n") and not line.endswith("//modDef\n"):
            f.write(line)
        if line.endswith("//FuncDef\n"):
            f.write(
                'static PyMethodDef PyTrans' + name + '_funcs[] = {{"H", (PyCFunction)MT_H,    METH_VARARGS, PyTrans_docs},{"nF", (PyCFunction)MT_fieldNumber,        METH_VARARGS, PyTrans_docs},{"nP", (PyCFunction)MT_paramNumber,        METH_VARARGS, PyTrans_docs},{"V", (PyCFunction)MT_V,            METH_VARARGS, PyTrans_docs},{"dV", (PyCFunction)MT_dV,                METH_VARARGS, PyTrans_docs},  {"ddV", (PyCFunction)MT_ddV,                METH_VARARGS, PyTrans_docs}, {"dotfieldsSR", (PyCFunction)MT_dotfieldsSR,        METH_VARARGS, PyTrans_docs}, {"massMatrix", (PyCFunction)MT_massMatrix,        METH_VARARGS, PyTrans_docs}, {"findEndOfInflation", (PyCFunction)MT_findEndOfInflation,        METH_VARARGS, PyTrans_docs}, {"backEvolve", (PyCFunction)MT_backEvolve,        METH_VARARGS, PyTrans_docs},    {"sigEvolve", (PyCFunction)MT_sigEvolve,        METH_VARARGS, PyTrans_docs},    {"alphaEvolve", (PyCFunction)MT_alphaEvolve,        METH_VARARGS, PyTrans_docs},    {NULL}};//FuncDef\n')

        if line.endswith("//modDef\n"):
            f.write(
                'static struct PyModuleDef PyTransModule = {PyModuleDef_HEAD_INIT, "PyTrans' + name + '", PyTrans_docs, -1, PyTrans' + name + '_funcs}; //modDef\n')

        if line.endswith("//initFunc\n"):
            f.write(
                'PyMODINIT_FUNC PyInit_PyTrans' + name + '(void)    {    PyObject *m = PyModule_Create(&PyTransModule); import_array(); return m;} //initFunc\n')
    f.close()

    my_env = os.environ.copy()
    my_env["PYTHONUSERBASE"] = location
    p = subprocess.Popen(["python", filename1, "install", "--user"], cwd=location, stdin=subprocess.PIPE,
                         stderr=subprocess.PIPE, env=my_env)

    stdout, stderr = p.communicate()

    pathSet()

    shutil.rmtree(os.path.join(location, 'build'), ignore_errors=True)


def deleteModule(name):
    location = os.path.join(dir, 'PyTrans')
    [os.remove(os.path.join(location, f)) for f in os.listdir(location) if f.startswith("PyTrans" + name)]


def tol(rtol, atol):
    dir = os.path.dirname(__file__)
    filename = os.path.join(dir, 'PyTrans', 'PyTrans.cpp')
    f = open(filename, "r")


def potential(V, nF, nP, simpleGeometric=False, simplePotentials=False, G="canonical", silent=True):
    # Define symbols for fields and parameters
    f = sym.symarray('f', nF)
    p = sym.symarray('p', nP)

    # build curvature records directory if not found
    dir = os.path.dirname(__file__)
    curv_dir = os.path.join(dir, 'PyTrans', 'CurvatureRecords')
    if not os.path.exists(curv_dir):
        os.makedirs(curv_dir)

    # Make a temporary name for the curvature file (will be renamed by in compile process)
    curv_tmp = os.path.join(curv_dir, "TEMP.curvature")
    if not silent:
        timer = t.process_time()
        print('[{time}] constructing curvature class instance'.format(time=t.ctime()))

    # Initialize curvature class instance: This will handle all symbolic computations requried
    """ TODO: add kwargs for simplifying curv. and pots. """
    curv_instance = gravtools_pyt.curvatureObject(
        G, f, V, params=p, simpleGeometric=simpleGeometric, simplePotentials=simplePotentials)

    if not silent:
        print('[{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer))

    # Save record of curvature class (this can then be callable from PyTransScripts)
    curv_file = open(curv_tmp, 'wb')
    with curv_file as cf:
        pk.dump(curv_instance, cf)

    if not silent:
        timer = t.process_time()
        print('[{time}] writing curvature expressions'.format(time=t.ctime()))

    # Pass curvature instance to fieldmetric function: Will write PyTrans files
    fieldmetric(G, nF, nP, simple=simpleGeometric, silent=silent, curv_obj=curv_instance)

    # Construct (flattened) symbolic arrays to hold 2st, 2nd and 3rd order derivatives
    vd = sym.symarray('vd', nF)
    vdd = sym.symarray('vdd', nF * nF)
    vddd = sym.symarray('vddd', nF * nF * nF)

    # Load precomputed arrays and remap to flattened configuration
    dV_arr = curv_instance.dV
    ddV_arr = curv_instance.ddV
    dddV_arr = curv_instance.dddV

    # We then flatten these arrays and (hopefully) match indices for the cpp writer
    ddV_arr_flat = ddV_arr.flatten()
    dddV_arr_flat = dddV_arr.T.flatten()  # Transpose of canonical matrix coords. + flatten seems to agree

    for i in range(nF):
        vd[i] = dV_arr[i]
        for j in range(nF):
            vdd[i + j * nF] = ddV_arr_flat[i + j * nF]
            for k in range(nF):
                vddd[i + j * nF + k * nF * nF] = dddV_arr_flat[i + j * nF + k * nF * nF]

    if not silent:
        print('[{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer))

    dir = os.path.dirname(__file__)
    filename1 = os.path.join(dir, 'CppTrans', 'potentialProto.h')
    filename2 = os.path.join(dir, 'CppTrans', 'potential.h')
    f = open(filename1, 'r')
    g = open(filename2, 'w')

    if not silent:
        timer = t.process_time()
        print('[{time}] writing to potential.h'.format(time=t.ctime()))

    for line in f:

        g.write(line)

        if line == "// #Rewrite\n":
            g.write('// Potential file rewriten at' + ' ' + t.strftime("%c") + '\n')

        if line == "// #FP\n":
            g.write('nF=' + str(nF) + ';\n' + 'nP=' + str(nP) + ';\n')

        if line == "// Pot\n":

            # extract common subexpressions from V
            if not silent:
                timer_cse = t.process_time()
                print('[{time}] performing CSE for V'.format(time=t.ctime()))
            decls, new_expr = sym.cse(V, order='none')
            if not silent:
                print('[{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer_cse))

            # emit declarations for CSE variables
            write_cse_decls(decls, g, nF, nP)

            # emit main expression
            emit_expr = sym.printing.cxxcode(new_expr[0])
            rw_expr = rewrite_indices(emit_expr, nF, nP)
            g.write('  sum=' + str(rw_expr) + ';\n')

        if line == "// dPot\n":

            if not silent:
                timer_cse = t.process_time()
                print('[{time}] performing CSE for dV'.format(time=t.ctime()))
            decls, new_exprs = sym.cse(vd, order='none')
            if not silent:
                print('[{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer_cse))

            # emit declarations for CSE variables
            write_cse_decls(decls, g, nF, nP)

            for i in range(nF):
                emit_expr = sym.printing.cxxcode(new_exprs[i])
                rw_expr = rewrite_indices(emit_expr, nF, nP)
                g.write('\n sum[' + str(i) + ']=' + str(rw_expr) + ';\n')

        if line == "// ddPot\n":

            if not silent:
                timer_cse = t.process_time()
                print('[{time}] performing CSE for ddV'.format(time=t.ctime()))
            decls, new_exprs = sym.cse(vdd, order='none')
            if not silent:
                print('[{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer_cse))

            # emit declarations for CSE variables
            write_cse_decls(decls, g, nF, nP)

            for i in range(nF):
                for j in range(nF):
                    emit_expr = sym.printing.cxxcode(new_exprs[i + nF * j])
                    rw_expr = rewrite_indices(emit_expr, nF, nP)
                    g.write('\n sum[' + str(i + nF * j) + ']=' + str(rw_expr) + ';\n')

        if line == "// dddPot\n":

            if not silent:
                timer_cse = t.process_time()
                print('[{time}] performing CSE for dddV'.format(time=t.ctime()))
            decls, new_exprs = sym.cse(vddd, order='none')
            if not silent:
                print('[{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer_cse))

            # emit declarations for CSE variables
            write_cse_decls(decls, g, nF, nP)

            for i in range(nF):
                for j in range(nF):
                    for k in range(nF):
                        emit_expr = sym.printing.cxxcode(new_exprs[i + nF * j + nF * nF * k])
                        rw_expr = rewrite_indices(emit_expr, nF, nP)
                        g.write('\n sum[' + str(i + nF * j + nF * nF * k) + ']=' + str(rw_expr) + ';\n')

    if not silent:
        print('[{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer))

    g.close()
    f.close()


def write_cse_decls(decls, g, nF, nP):
    # emit declarations for common subexpressions
    for rule in decls:
        symb = sym.printing.cxxcode(rule[0])
        expr = sym.printing.cxxcode(rule[1])
        new_expr = rewrite_indices(expr, nF, nP)
        g.write('  auto ' + symb + ' = ' + new_expr + ';\n')


def rewrite_indices(expr, nF, nP):
    new_expr = expr

    for l in range(max(nP, nF)):
        l = max(nP, nF) - 1 - l
        new_expr = new_expr.replace("_" + str(l), "[" + str(l) + "]")

    return new_expr


def fieldmetric(G, nF, nP, simple=False, silent=True, curv_obj=None):
    f = sym.symarray('f', nF)
    p = sym.symarray('p', nP)

    Ga = curv_obj._christoffelSymbols
    Ri = 0  # This isn't actually used...
    Rm = curv_obj._riemannTensor
    covDRm = curv_obj._riemannTensorCovD

    dir = os.path.dirname(__file__)
    filename1 = os.path.join(dir, 'CppTrans', 'fieldmetricProto.h')
    filename2 = os.path.join(dir, 'CppTrans', 'fieldmetric.h')
    e = open(filename1, 'r')
    h = open(filename2, 'w')

    # Build symbolic arrays to hold quantities: These will be converted into cpp expressions later
    G_array = curv_obj.G_array
    Gamma_array = sym.symarray('Gamma', 2 * nF * 2 * nF * 2 * nF)
    R_array = sym.symarray('Riemann', nF * nF * nF * nF)
    gradR_array = sym.symarray('gradRiemann', nF * nF * nF * nF * nF)

    # populate connexion matrix
    for i in range(2 * nF):
        for j in range(2 * nF):
            for k in range(2 * nF):
                if i < nF:
                    ii = -i - 1
                else:
                    ii = i - (nF - 1)
                if j < nF:
                    jj = -j - 1
                else:
                    jj = j - (nF - 1)
                if k < nF:
                    kk = -k - 1
                else:
                    kk = k - (nF - 1)

                if kk < 0 or jj < 0 or ii > 0:
                    Gamma_array[(2 * nF) * (2 * nF) * i + (2 * nF) * j + k] = sym.simplify(0)
                else:
                    Gamma_array[(2 * nF) * (2 * nF) * i + (2 * nF) * j + k] = Ga(
                        ii, jj, kk) if curv_obj is None else curv_obj.getChristoffel(abs(ii) - 1, jj - 1, kk - 1)

                    if curv_obj is None and simple is True:
                        Gamma_array[(2 * nF) * (2 * nF) * i + (2 * nF) * j + k] = sym.simplify(
                            Gamma_array[(2 * nF) * (2 * nF) * i + (2 * nF) * j + k]
                        )

    # populate Riemann matrix
    for i in range(nF):
        for j in range(nF):
            for k in range(nF):
                for l in range(nF):
                    ii = i + 1
                    jj = j + 1
                    kk = k + 1
                    ll = l + 1

                    if curv_obj is None:
                        if simple is True:
                            R_array[(nF) * (nF) * (nF) * i + (nF) * (nF) * j + (nF) * k + l] = sym.simplify(
                                Rm(ii, jj, kk, ll))
                        else:
                            R_array[(nF) * (nF) * (nF) * i + (nF) * (nF) * j + (nF) * k + l] = Rm(ii, jj, kk, ll)
                    else:
                        R_array[(nF) * (nF) * (nF) * i + (nF) * (nF) * j + (nF) * k + l] = Rm[
                            ii - 1, jj - 1, kk - 1, ll - 1]

    # populate covariant-derivative of Riemann matrix
    for i in range(nF):
        for j in range(nF):
            for k in range(nF):
                for l in range(nF):
                    for m in range(nF):
                        ii = i + 1
                        jj = j + 1
                        kk = k + 1
                        ll = l + 1
                        mm = m + 1

                        if curv_obj is None:
                            if simple is True:
                                gradR_array[(nF) * (nF) * (nF) * (nF) * i + (nF) * (nF) * (nF) * j + (nF) * (nF) * k + (
                                    nF) * l + m] = sym.simplify(Rm.covariantD(ii, jj, kk, ll, mm))
                            else:
                                gradR_array[(nF) * (nF) * (nF) * (nF) * i + (nF) * (nF) * (nF) * j + (nF) * (nF) * k + (
                                    nF) * l + m] = Rm.covariantD(ii, jj, kk, ll, mm)

                        else:
                            gradR_array[(nF) * (nF) * (nF) * (nF) * i + (nF) * (nF) * (nF) * j + (nF) * (nF) * k + (
                                nF) * l + m] = covDRm[ii - 1, jj - 1, kk - 1, ll - 1, mm - 1]

    for line in e:
        h.write(line)
        if line == "// #FP\n":
            # h.write('nF='+str(nF)+';\n')
            h.write('nF=' + str(nF) + ';\n' + 'nP=' + str(nP) + ';\n')

        if line == "// metric\n":

            if not silent:
                timer_cse = t.process_time()
                print('[{time}] performing CSE for field metric'.format(time=t.ctime()))
            decls, new_expr = sym.cse(G_array, order='none')
            if not silent:
                print('[{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer_cse))

            # emit declarations for CSE variables
            write_cse_decls(decls, h, nF, nP)

            for i in range(2 * nF):
                for j in range(2 * nF):
                    # emit main expression
                    emit_expr = sym.printing.cxxcode(new_expr[(2 * nF) * i + j])
                    rw_expr = rewrite_indices(emit_expr, nF, nP)
                    h.write('\n FM[' + str((2 * nF) * i + j) + ']=' + str(rw_expr) + ';\n')

        if line == "// Christoffel\n":

            if not silent:
                timer_cse = t.process_time()
                print('    [{time}] performing CSE for Christoffel symbols'.format(time=t.ctime()))
            decls, new_expr = sym.cse(Gamma_array, order='none')
            if not silent:
                print('[{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer_cse))

            # emit declarations for CSE variables
            write_cse_decls(decls, h, nF, nP)

            for i in range(2 * nF):
                for j in range(2 * nF):
                    for k in range(2 * nF):
                        # emit main expression
                        emit_expr = sym.printing.cxxcode(new_expr[(2 * nF) * (2 * nF) * i + (2 * nF) * j + k])
                        rw_expr = rewrite_indices(emit_expr, nF, nP)
                        h.write(
                            '\n CS[' + str((2 * nF) * (2 * nF) * i + (2 * nF) * j + k) + ']=' + str(rw_expr) + ';\n')

        if line == "// Riemann\n":

            if not silent:
                timer_cse = t.process_time()
                print('[{time}] performing CSE for Riemann tensor'.format(time=t.ctime()))

            decls, new_expr = sym.cse(R_array, order='none')
            if not silent:
                print('[{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer_cse))

            # emit declarations for CSE variables
            write_cse_decls(decls, h, nF, nP)

            for i in range(nF):
                for j in range(nF):
                    for k in range(nF):
                        for l in range(nF):
                            # emit main expression
                            emit_expr = sym.printing.cxxcode(
                                new_expr[(nF) * (nF) * (nF) * i + (nF) * (nF) * j + (nF) * k + l])
                            rw_expr = rewrite_indices(emit_expr, nF, nP)
                            h.write(
                                '\n RM[' + str((nF) * (nF) * (nF) * i + (nF) * (nF) * j + (nF) * k + l) + ']=' + str(
                                    rw_expr) + ';\n')

        if line == "// Riemanncd\n":

            if not silent:
                timer_cse = t.process_time()
                print('[{time}] performing CSE for Riemann tensor'.format(time=t.ctime()))
            decls, new_expr = sym.cse(gradR_array, order='none')
            if not silent:
                print('[{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer_cse))

            # emit declarations for CSE variables
            write_cse_decls(decls, h, nF, nP)

            for i in range(nF):
                for j in range(nF):
                    for k in range(nF):
                        for l in range(nF):
                            for m in range(nF):
                                # emit main expression
                                emit_expr = sym.printing.cxxcode(new_expr[
                                                                     (nF) * (nF) * (nF) * (nF) * i + (nF) * (nF) * (
                                                                         nF) * j + (nF) * (nF) * k + (nF) * l + m])
                                rw_expr = rewrite_indices(emit_expr, nF, nP)
                                h.write('\n RMcd[' + str(
                                    (nF) * (nF) * (nF) * (nF) * i + (nF) * (nF) * (nF) * j + (nF) * (nF) * k + (
                                        nF) * l + m) + ']=' + str(rw_expr) + ';\n')

    h.close()
    e.close()

    # return g, Ga, Ri, Rm
