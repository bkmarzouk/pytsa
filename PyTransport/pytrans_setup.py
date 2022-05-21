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


# This file contains python scripts used to setup the compiled pyt module
import sympy as sym
import subprocess
import os
import shutil
import time as t
from pytransport import sym_tools

_PyMethodDefs = ",".join([
    '{"H", (PyCFunction)MT_H,   METH_VARARGS, "Compute Hubble rate"}',
    '{"nF", (PyCFunction)MT_fieldNumber,   METH_NOARGS, "Get number of fields for model"}',
    '{"nP", (PyCFunction)MT_paramNumber,   METH_NOARGS, "Get number of params for model"}',
    '{"Epsilon", (PyCFunction)MT_Ep,   METH_VARARGS, "Compute Epsilon slow roll value"}',
    '{"Eta", (PyCFunction)MT_Eta,   METH_VARARGS, "Compute Eta slow roll value"}',
    '{"V", (PyCFunction)MT_V,   METH_VARARGS, "Compute potential"}',
    '{"dV", (PyCFunction)MT_dV,   METH_VARARGS, "Compute 1st derivative of potential"}',
    '{"ddV", (PyCFunction)MT_ddV,   METH_VARARGS, "Compute 2nd derivative of potential"}',
    '{"dotfieldsSR", (PyCFunction)MT_dotfieldsSR,   METH_VARARGS, "REMOVE"}',  # REMOVE DEP
    '{"massMatrix", (PyCFunction)MT_massMatrix,   METH_VARARGS, "Compute mass-matrix"}',
    '{"findEndOfInflation", (PyCFunction)MT_findEndOfInflation,   METH_VARARGS, "End of inflation (SR violation)"}',
    '{"backEvolve", (PyCFunction)MT_backEvolve,   METH_VARARGS, "Compute background evolution"}',
    '{"sigEvolve", (PyCFunction)MT_sigEvolve,   METH_VARARGS, "Evolve 2pt correlation functions"}',
    '{"alphaEvolve", (PyCFunction)MT_alphaEvolve,   METH_VARARGS, "Evolve 3pt correlation functions"}',
])


def _delta_ctime(a, b):
    """
    Translates time.ctime() outputs into time difference in seconds

    :param a: inital time
    :param b: final time
    :return: seconds difference
    """
    h = lambda t: int(t.split(":")[0][-2:])
    m = lambda t: int(t.split(":")[1])
    s = lambda t: int(t.split(":")[2][:2])

    i = h(a) * 60 * 60 + m(a) * 60 + s(a)
    j = h(b) * 60 * 60 + m(b) * 60 + s(b)

    return j - i


def set_template_headers(NC: bool):
    """
    Configures directory structure for pytransport

    :param NC: if True methods use non-canonical field metric
    :return: None, configure pyt.cpp file
    """

    # Root of all the cpp sources
    source_dir = os.path.dirname(__file__)

    # PyTransport compilation file
    pytrans_path = os.path.join(source_dir, 'pyt', 'PyTrans.cpp')

    # This is where the header files are kept (simpler versions avail if not NC)
    cppt_dir = os.path.join(source_dir, "cppt", "NC") if NC else os.path.join(source_dir, "cppt")

    # This is common between NC and not NC
    stepper_path = os.path.join(source_dir, 'cppt', 'stepper', 'rkf45.hpp')

    # Load lines from template
    with open(pytrans_path, "r") as f:
        lines = f.readlines()

    with open(pytrans_path, "w") as f:

        # Update header file names accordingly
        for line in lines:
            if line.endswith("//evolve\n"):
                header_fpath = os.path.join(cppt_dir, 'evolve.h')
                f.write('#include' + '"' + header_fpath + '"' + '//evolve' + '\n')
            elif line.endswith("//moments\n"):
                header_fpath = os.path.join(cppt_dir, 'moments.h')
                f.write('#include' + '"' + header_fpath + '"' + '//moments' + '\n')
            elif line.endswith("//model\n"):
                header_fpath = os.path.join(cppt_dir, 'model.h')
                f.write('#include' + '"' + header_fpath + '"' + '//model' + '\n')
            elif line.endswith("//stepper\n"):  # invariant to NC
                f.write('#include' + '"' + stepper_path + '"' + '//stepper' + '\n')
            else:
                f.write(line)


def set_paths():  # Remove dep
    pass


def compile_module(name, NC=False):
    """
    Compiles template module

    :param name: module name
    :param NC: Non-canonical field metric if True, else Euclidean
    """

    t_start = t.ctime()
    print('   [{time}] start'.format(time=t_start))

    set_template_headers(NC)
    cwd = os.path.dirname(__file__)
    location = os.path.join(cwd, 'pyt')
    setup_file_path = os.path.join(cwd, 'pyt', 'moduleSetup.py')

    with open(setup_file_path, "r") as f:
        module_setup_lines = f.readlines()

    # Rewrite module name in template file
    with open(setup_file_path, "w") as f:
        for line in module_setup_lines:
            if "PYT_MODNAME" in line:
                f.write("mod_name = '{}'  # PYT_MODNAME\n".format(name))
            else:
                f.write(line)

    # Load PyTransport C++ file
    pytrans_cpp_path = os.path.join(cwd, 'pyt', 'PyTrans.cpp')
    with open(pytrans_cpp_path, "r") as f:
        pytrans_cpp_lines = f.readlines()

    with open(pytrans_cpp_path, "w") as f:
        for line in pytrans_cpp_lines:
            if not line.endswith("//FuncDef\n") and not line.endswith("//initFunc\n") and not line.endswith(
                    "//modDef\n"):
                f.write(line)
            if line.endswith("//FuncDef\n"):
                f.write('static PyMethodDef ' + name + '_funcs[] = {' + _PyMethodDefs + ',   {NULL}};//FuncDef\n')

            if line.endswith("//modDef\n"):
                f.write('static struct PyModuleDef PyTransModule = {PyModuleDef_HEAD_INIT, "' +
                        name + '", PyTrans_docs, -1, ' + name + '_funcs}; //modDef\n')

            if line.endswith("//initFunc\n"):
                f.write('PyMODINIT_FUNC PyInit_' +
                        name + '(void)    {import_array();  return PyModule_Create(&PyTransModule);} //initFunc\n')

    t_start_compile = t.ctime()
    print('   [{time}] start compile phase'.format(time=t_start_compile))

    import numpy as np
    os.system("export CFLAGS='-I {}'".format(np.get_include()))
    # install_lib = os.path.join(os.path.dirname(__file__), "pyt")
    # subprocess.run(["python", setup_file_path, "install", "--prefix={}".format(install_lib)],
    #                cwd=location)
    subprocess.run(["python", setup_file_path, "install"],
                   cwd=location)

    t_end = t.ctime()
    print("\n-- Compiled source in {} seconds, total time {} seconds".format(_delta_ctime(t_start_compile, t_end),
                                                                             _delta_ctime(t_start, t_end)))

    try:
        shutil.rmtree(os.path.join(location, 'build'), ignore_errors=False)
    except:
        pass


def delete_module(name):
    """
    Deletes module based on compile name

    :param name: module name
    """
    location = os.path.join(os.path.dirname(__file__), 'pyt')
    [os.remove(os.path.join(location, f)) for f in os.listdir(location) if f.startswith("pyt" + name)]


# TODO: Tidy translator; wrap into single call?


def potential(V, nF, nP, simplify_fmet=False, simplify_pot=False, simplify_covd=False, G: sym.Matrix or None = None,
              silent=False, recache=False):
    """
    Wrapper for compile source for potential, field metric and derivatives / combinations thereof

    :param V: Symbolic potential
    :param nF: Number of fields
    :param nP: Number of parameters
    :param simplify_fmet: if True, simplifies symbolic expressions related to field space metric
    :param simplify_pot: if True, simplifies symbolic expressions related to the potential
    :param simplify_covd: if True, simplifies symbolic expressions for covariant derivatives
    :param G: Symbolic field space metric
    :param silent: if True, no output
    :param recache: uf True, recomputes results
    :return:
    """

    # Get symbolic expressions for covariant derivatives
    covd_sym = sym_tools.CovDSym(nF, nP, G, V, simplify_fmet=simplify_fmet, simplify_pot=simplify_pot,
                                 simplify=simplify_covd, recache=recache)

    fieldmetric(nF, nP, G, V, recache=recache, simple_fmet=simplify_fmet, simple_potential=simplify_pot,
                simple_covd=simplify_covd, silent=silent)

    v, vd, vdd, vddd = covd_sym.get_potential_sym_arrays()

    dir = os.path.dirname(__file__)
    filename1 = os.path.join(dir, 'cppt', 'potentialProto.h')
    filename2 = os.path.join(dir, 'cppt', 'potential.h')
    f = open(filename1, 'r')
    g = open(filename2, 'w')

    if not silent:
        timer = t.process_time()
        print('   [{time}] writing to potential.h'.format(time=t.ctime()))

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
                print('   [{time}] performing CSE for V'.format(time=t.ctime()))

            decls, new_expr = sym.cse(V, order='none')
            if not silent:
                print('   [{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer_cse))

            # emit declarations for CSE variables
            write_cse_decls(decls, g, nF, nP)

            # emit main expression
            emit_expr = sym.printing.cxxcode(new_expr[0])
            rw_expr = rewrite_indices(emit_expr, nF, nP)
            g.write('  sum=' + str(rw_expr) + ';\n')

        if line == "// dPot\n":
            if not silent:
                timer_cse = t.process_time()
                print('   [{time}] performing CSE for dV'.format(time=t.ctime()))

            decls, new_exprs = sym.cse(vd, order='none')
            if not silent:
                print('   [{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer_cse))

            # emit declarations for CSE variables
            write_cse_decls(decls, g, nF, nP)

            for i in range(nF):
                emit_expr = sym.printing.cxxcode(new_exprs[i])
                rw_expr = rewrite_indices(emit_expr, nF, nP)
                g.write('\n sum[' + str(i) + ']=' + str(rw_expr) + ';\n')

        if line == "// ddPot\n":

            if not silent:
                timer_cse = t.process_time()
                print('   [{time}] performing CSE for ddV'.format(time=t.ctime()))

            decls, new_exprs = sym.cse(vdd, order='none')
            if not silent:
                print('   [{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer_cse))

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
                print('   [{time}] performing CSE for dddV'.format(time=t.ctime()))

            decls, new_exprs = sym.cse(vddd, order='none')
            if not silent:
                print('   [{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer_cse))

            # emit declarations for CSE variables
            write_cse_decls(decls, g, nF, nP)

            for i in range(nF):
                for j in range(nF):
                    for k in range(nF):
                        emit_expr = sym.printing.cxxcode(new_exprs[i + nF * j + nF * nF * k])
                        rw_expr = rewrite_indices(emit_expr, nF, nP)
                        g.write('\n sum[' + str(i + nF * j + nF * nF * k) + ']=' + str(rw_expr) + ';\n')

    if not silent:
        print('   [{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer))

    g.close()
    f.close()


def write_cse_decls(decls, g, nF, nP):
    """
    Wrapper for performing CSE on a set of declarations

    :param decls: rules / declarations
    :param g: field space metric
    :param nF: number of fields
    :param nP: number of paramters
    """
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
    """
    Prepare and compile codes relevant to field space metric

    :param nF: Number of fields
    :param nP: Number of parameters
    :param G: Field space metric (symbolic)
    :param V: Potential (symbolic)
    :param recache: if True, recomputes / recaches symbolic expressions
    :param simple_fmet: if True, simplifies expressions related to field space
    :param simple_potential: if True, simplifies expressions related to potential
    :param simple_covd: if True, simplifies expressions of covariant derivatives
    :param silent: if True, not output is printed
    """
    dir = os.path.dirname(__file__)
    filename1 = os.path.join(dir, 'cppt', 'fieldmetricProto.h')
    filename2 = os.path.join(dir, 'cppt', 'fieldmetric.h')
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
                print('   [{time}] performing CSE for field metric'.format(time=t.ctime()))
            decls, new_expr = sym.cse(G_array, order='none')
            if not silent:
                print('   [{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer_cse))

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
                print('   [{time}] performing CSE for Christoffel symbols'.format(time=t.ctime()))
            decls, new_expr = sym.cse(Gamma_array, order='none')
            if not silent:
                print('   [{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer_cse))

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
                print('   [{time}] performing CSE for Riemann tensor'.format(time=t.ctime()))

            decls, new_expr = sym.cse(R_array, order='none')
            if not silent:
                print('   [{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer_cse))

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
                print('   [{time}] performing CSE for Riemann tensor'.format(time=t.ctime()))
            decls, new_expr = sym.cse(gradR_array, order='none')
            if not silent:
                print('   [{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer_cse))

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
