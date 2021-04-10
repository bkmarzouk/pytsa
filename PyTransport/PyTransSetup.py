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
import time

import sympy as sym
import sys
import site
import subprocess
import platform
import os
import shutil
import time as t
from PyTransport import sym_tools

try:
    import multiprocessing

    n_cores = multiprocessing.cpu_count()
    j_avail = True
except ImportError:

    multiprocessing = None
    n_cores = 0
    j_avail = False


def compile_str(name, use_j=False):
    if not use_j:
        return 'setup(name="PyTrans' + name + '", version="1.0", ext_modules=[Extension("PyTrans' + name + '", [filename, filename2 ])], include_dirs=[numpy.get_include(), dirs], extra_compile_args = ["-std=c++11 -frounding-math -fsignaling-nans"])#setup\n'

    j_str = " -j {}".format(int(1.5 * n_cores))

    return 'setup(name="PyTrans' + name + '", version="1.0", ext_modules=[Extension("PyTrans' + name + '", [filename, filename2 ])], include_dirs=[numpy.get_include(), dirs], extra_compile_args = ["-std=c++11 -frounding-math -fsignaling-nans' + j_str + '"])#setup\n'


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

def delta_ctime(a, b):

    h = lambda t: int(t.split(":")[0][-2:])
    m = lambda t: int(t.split(":")[1])
    s = lambda t: int(t.split(":")[2][:2])

    i = h(a) * 60 * 60 + m(a) * 60 + s(a)
    j = h(b) * 60 * 60 + m(b) * 60 + s(b)

    return j - i

def compileName(name, NC=False, use_j=False):

    t_start = t.ctime()
    print('[{time}] compile start'.format(time=t_start))

    directory(NC)
    dir = os.path.dirname(__file__)
    location = os.path.join(dir, 'PyTrans')
    filename1 = os.path.join(dir, 'PyTrans', 'moduleSetup.py')

    if use_j and j_avail:
        use_j = True
    else:
        use_j = False

    # if NC is True:
    #     # There should exist a temporary file fom the calling 'potential'
    #     curv_path = os.path.join(dir, 'PyTrans', 'CurvatureRecords', 'TEMP.curvature')
    #     assert os.path.exists(curv_path), "Temporary curvature file not found! {}".format(curv_path)
    #
    #     # Rename temporary curvature file with compilation name!
    #     new_curv_path = curv_path.replace('TEMP', name)
    #     os.rename(curv_path, new_curv_path)

    # assert 0, compile_str(name, use_j=use_j)

    f = open(filename1, "r")
    lines = f.readlines()
    f.close()

    f = open(filename1, "w")
    for line in lines:
        if not line.endswith("#setup\n"):
            f.write(line)
        if line.endswith("#setup\n"):
            f.write(compile_str(name, use_j=use_j))
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

    t_end = t.ctime()

    print('[{time}] compile complete'.format(time=t_end))

    print("\n-- Compiled source in {} seconds".format(delta_ctime(t_start, t_end)))


def deleteModule(name):
    location = os.path.join(dir, 'PyTrans')
    [os.remove(os.path.join(location, f)) for f in os.listdir(location) if f.startswith("PyTrans" + name)]


def tol(rtol, atol):
    dir = os.path.dirname(__file__)
    filename = os.path.join(dir, 'PyTrans', 'PyTrans.cpp')
    f = open(filename, "r")


def potential(V, nF, nP, simplify_fmet=False, simplify_pot=False, simplify_covd=False, G: sym.Matrix or None = None,
              silent=False, recache=False):
    # Define symbols for fields and parameters
    # f = sym.symarray('f', nF)
    # p = sym.symarray('p', nP)
    #
    # # build curvature records directory if not found
    # dir = os.path.dirname(__file__)
    # curv_dir = os.path.join(dir, 'PyTrans', 'CurvatureRecords')
    # if not os.path.exists(curv_dir):
    #     os.makedirs(curv_dir)

    covd_sym = sym_tools.CovDSym(nF, nP, G, V, simplify_fmet=simplify_fmet, simplify_pot=simplify_pot,
                                 simplify=simplify_covd, recache=recache)

    fieldmetric(nF, nP, G, V, recache=recache, simple_fmet=simplify_fmet, simple_potential=simplify_pot,
                simple_covd=simplify_covd, silent=silent)

    # fmet_sym = covd_sym.fmet_sym
    # pot_sym = covd_sym.pot_sym

    # # Make a temporary name for the curvature file (will be renamed by in compile process)
    # curv_tmp = os.path.join(curv_dir, "TEMP.curvature")
    # if not silent:
    #     timer = t.process_time()
    #     print('[{time}] constructing curvature class instance'.format(time=t.ctime()))
    #
    # # Initialize curvature class instance: This will handle all symbolic computations requried
    # """ TODO: add kwargs for simplifying curv. and pots. """
    # curv_instance = gravtools_pyt.curvatureObject(
    #     G, f, V, params=p, simpleGeometric=simpleGeometric, simplePotentials=simplePotentials)
    #
    # if not silent:
    #     print('[{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer))
    #
    # # Save record of curvature class (this can then be callable from PyTransScripts)
    # curv_file = open(curv_tmp, 'wb')
    # with curv_file as cf:
    #     pk.dump(curv_instance, cf)

    # if not silent:
    #     timer = t.process_time()
    #     print('[{time}] writing curvature expressions'.format(time=t.ctime()))
    #
    # # Pass curvature instance to fieldmetric function: Will write PyTrans files
    # fieldmetric(G, nF, nP, simple=simpleGeometric, silent=silent, curv_obj=curv_instance)
    #
    # # Construct (flattened) symbolic arrays to hold 2st, 2nd and 3rd order derivatives
    # vd = sym.symarray('vd', nF)
    # vdd = sym.symarray('vdd', nF * nF)
    # vddd = sym.symarray('vddd', nF * nF * nF)
    #
    # # Load precomputed arrays and remap to flattened configuration
    # dV_arr = curv_instance.dV
    # ddV_arr = curv_instance.ddV
    # dddV_arr = curv_instance.dddV

    # # We then flatten these arrays and (hopefully) match indices for the cpp writer
    # ddV_arr_flat = ddV_arr.flatten()
    # dddV_arr_flat = dddV_arr.T.flatten()  # Transpose of canonical matrix coords. + flatten seems to agree
    #
    # for i in range(nF):
    #     vd[i] = dV_arr[i]
    #     for j in range(nF):
    #         vdd[i + j * nF] = ddV_arr_flat[i + j * nF]
    #         for k in range(nF):
    #             vddd[i + j * nF + k * nF * nF] = dddV_arr_flat[i + j * nF + k * nF * nF]

    # if not silent:
    #     print('[{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer))

    v, vd, vdd, vddd = covd_sym.get_potential_sym_arrays()

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


def fieldmetric(nF, nP, G, V, recache=False, simple_fmet=False, simple_potential=False, simple_covd=False, silent=True):
    dir = os.path.dirname(__file__)
    filename1 = os.path.join(dir, 'CppTrans', 'fieldmetricProto.h')
    filename2 = os.path.join(dir, 'CppTrans', 'fieldmetric.h')
    e = open(filename1, 'r')
    h = open(filename2, 'w')

    covd_sym = sym_tools.CovDSym(nF, nP, G, V, recache, simple_fmet, simple_potential, simple_covd)

    G_array, Gamma_array, R_array, gradR_array = covd_sym.get_curvature_sym_arrays()

    for line in e:
        h.write(line)
        if line == "// #FP\n":
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
                print('[{time}] performing CSE for Christoffel symbols'.format(time=t.ctime()))
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
