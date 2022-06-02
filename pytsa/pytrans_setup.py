# This file is part of PyTransport.

# PyTransport is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# PyTransport is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with PyTransport.  If not, see <http://www.gnu.org/licenses/>.

import sympy as sym
import subprocess
import os
import shutil
import time as t
import numpy as np
from . import sym_tools

# 1st in {} defines function name for python calls
# 2nd in {} defines function name in C++ extension
# 3rd in {} defines the input for the method; e.g. METH_NOARGS -> no arguments expected for function
# 4th in {} defines the doc string, user can then yield info using help(<module>.method). docs define in PyTrans.cpp
_PyMethodDefs = ",".join([
    '{"H", (PyCFunction)MT_H,   METH_VARARGS, doc_H}',
    '{"nF", (PyCFunction)MT_fieldNumber,   METH_NOARGS, doc_nF}',
    '{"nP", (PyCFunction)MT_paramNumber,   METH_NOARGS, doc_nP}',
    '{"Epsilon", (PyCFunction)MT_Ep,   METH_VARARGS, doc_Epsilon}',
    '{"Eta", (PyCFunction)MT_Eta,   METH_VARARGS, doc_Eta}',
    '{"V", (PyCFunction)MT_V,   METH_VARARGS, doc_V}',
    '{"dV", (PyCFunction)MT_dV,   METH_VARARGS, doc_dV}',
    '{"ddV", (PyCFunction)MT_ddV,   METH_VARARGS, doc_ddV}',
    '{"massMatrix", (PyCFunction)MT_massMatrix,   METH_VARARGS, doc_massMatrix}',
    '{"findEndOfInflation", (PyCFunction)MT_findEndOfInflation,   METH_VARARGS, doc_findEndOfInflation}',
    '{"backEvolve", (PyCFunction)MT_backEvolve,   METH_VARARGS, doc_backEvolve}',
    '{"sigEvolve", (PyCFunction)MT_sigEvolve,   METH_VARARGS, doc_sigEvolve}',
    '{"alphaEvolve", (PyCFunction)MT_alphaEvolve,   METH_VARARGS, doc_alphaEvolve}',
])

cwd = os.path.abspath(os.path.dirname(__file__))

# Installation location for models
MODELS_LOC = os.path.join(cwd, "models")

pytrans_path = os.path.join(cwd, 'pyt', 'PyTrans.cpp')
stepper_path = os.path.join(cwd, 'cppt', 'stepper', 'rkf45.hpp')


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


def _set_template_headers(non_canonical: bool):
    """
    Rewrites header files used in compile script to point to (non-)canonical file versions

    :param non_canonical: if True methods use non-canonical
    """

    # This is where the header files are kept (simpler versions avail if not NC)
    cppt_dir = os.path.join(cwd, "cppt", "NC") if non_canonical else os.path.join(cwd, "cppt")

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


def _rewrite_indices(expr, nF, nP):
    """
    Rewrite symbolic expressions containing subsripts as indexed vector components
    e.g. f_0 -> f[0]

    :param expr: symbolic expression in string format
    :param nF: number of fields
    :param nP: number of params
    :return: translated expression
    """
    new_expr = expr

    for l in range(max(nP, nF)):
        l = max(nP, nF) - 1 - l
        new_expr = new_expr.replace("_" + str(l), "[" + str(l) + "]")

    return new_expr


def _write_cse_decls(decls, g, nF, nP):
    """
    Wrapper for performing CSE on a set of declarations

    :param decls: rules / declarations
    :param g: file (open in w mode)
    :param nF: number of fields
    :param nP: number of paramters
    """
    # emit declarations for common subexpressions
    for rule in decls:
        symb = sym.printing.cxxcode(rule[0])
        expr = sym.printing.cxxcode(rule[1])
        new_expr = _rewrite_indices(expr, nF, nP)
        g.write('  auto ' + symb + ' = ' + new_expr + ';\n')


class Translator:

    def __init__(self, nF: int, nP: int, V: sym.Expr, G: sym.Matrix or None = None,
                 simplify_potential=True, simplify_metric=True, simplify_covd=True,
                 verbose=True, recache=False):
        """
        Translates SymPy expressions into c++ source ready for compilation

        :param nF: number of fields
        :param nP: number of params
        :param V: sym expression for potential
        :param G: sym expression for fieldspace metric, may be None then euclidean metric assumed
        :param simplify_potential: if True, simplifies potential and partial derivatives of
        :param simplify_metric: if True, simplifies metric and inv. metric
        :param simplify_covd: if True, simplifies covariant derivative expressions
        :param verbose: if True, prints status as output
        :param recache: if True, rebuilds symbolic expressions, otherwise uses cached results
        """
        self.nF = nF
        self.nP = nP
        self.V = V
        self.G = G
        self.canonical = G is None

        self._translate_fmet(simplify_potential, simplify_metric, simplify_covd, verbose, recache)
        self._translate_pot(simplify_potential, simplify_metric, simplify_covd, verbose, recache)
        _set_template_headers(G is not None)

    def _translate_fmet(self, simplify_potential: bool, simplify_metric: bool, simplify_covd: bool, verbose: bool,
                        recache: bool):
        """
        Translate field space components
        """

        nF, nP = self.nF, self.nP

        covd_sym = sym_tools.CovDSym(nF, nP, self.G, self.V, recache, simplify_fmet=simplify_metric,
                                     simplify_pot=simplify_potential,
                                     simplify=simplify_covd)

        if not self.canonical:

            fmet_template_path = os.path.join(cwd, 'cppt', 'fieldmetricProto.h')
            fmet_translated_path = os.path.join(cwd, 'cppt', 'fieldmetric.h')

            with open(fmet_template_path, "r") as f:
                template_lines = f.readlines()

            if verbose:
                print("-- Translating field space symbols")

            with open(fmet_translated_path, "w") as h:

                G_array, Gamma_array, R_array, gradR_array = covd_sym.get_curvature_sym_arrays()

                for line in template_lines:
                    h.write(line)
                    if line == "// #FP\n":
                        h.write('nF=' + str(nF) + ';\n' + 'nP=' + str(nP) + ';\n')

                    if line == "// metric\n":

                        if verbose:
                            timer_cse = t.process_time()
                            print('   [{time}] performing CSE for field metric'.format(time=t.ctime()))
                        decls, new_expr = sym.cse(G_array, order='none')
                        if verbose:
                            print('   [{time}] complete in {x} sec'.format(time=t.ctime(),
                                                                           x=t.process_time() - timer_cse))

                        # emit declarations for CSE variables
                        _write_cse_decls(decls, h, nF, nP)

                        for i in range(2 * nF):
                            for j in range(2 * nF):
                                # emit main expression
                                emit_expr = sym.printing.cxxcode(new_expr[(2 * nF) * i + j])
                                rw_expr = _rewrite_indices(emit_expr, nF, nP)
                                h.write('\n FM[' + str((2 * nF) * i + j) + ']=' + str(rw_expr) + ';\n')

                    if line == "// Christoffel\n":

                        if verbose:
                            timer_cse = t.process_time()
                            print('   [{time}] performing CSE for Christoffel symbols'.format(time=t.ctime()))

                        decls, new_expr = sym.cse(Gamma_array, order='none')
                        if verbose:
                            print('   [{time}] complete in {x} sec'.format(time=t.ctime(),
                                                                           x=t.process_time() - timer_cse))

                        # emit declarations for CSE variables
                        _write_cse_decls(decls, h, nF, nP)

                        for i in range(2 * nF):
                            for j in range(2 * nF):
                                for k in range(2 * nF):
                                    # emit main expression
                                    emit_expr = sym.printing.cxxcode(
                                        new_expr[(2 * nF) * (2 * nF) * i + (2 * nF) * j + k])
                                    rw_expr = _rewrite_indices(emit_expr, nF, nP)
                                    h.write(
                                        '\n CS[' + str((2 * nF) * (2 * nF) * i + (2 * nF) * j + k) + ']=' + str(
                                            rw_expr) + ';\n')

                    if line == "// Riemann\n":

                        if verbose:
                            timer_cse = t.process_time()
                            print('   [{time}] performing CSE for Riemann tensor'.format(time=t.ctime()))

                        decls, new_expr = sym.cse(R_array, order='none')
                        if verbose:
                            print('   [{time}] complete in {x} sec'.format(time=t.ctime(),
                                                                           x=t.process_time() - timer_cse))

                        # emit declarations for CSE variables
                        _write_cse_decls(decls, h, nF, nP)

                        for i in range(nF):
                            for j in range(nF):
                                for k in range(nF):
                                    for l in range(nF):
                                        # emit main expression
                                        emit_expr = sym.printing.cxxcode(
                                            new_expr[(nF) * (nF) * (nF) * i + (nF) * (nF) * j + (nF) * k + l])
                                        rw_expr = _rewrite_indices(emit_expr, nF, nP)
                                        h.write(
                                            '\n RM[' + str(
                                                (nF) * (nF) * (nF) * i + (nF) * (nF) * j + (nF) * k + l) + ']=' + str(
                                                rw_expr) + ';\n')

                    if line == "// Riemanncd\n":

                        if verbose:
                            timer_cse = t.process_time()
                            print('   [{time}] performing CSE for Riemann tensor'.format(time=t.ctime()))
                        decls, new_expr = sym.cse(gradR_array, order='none')
                        if verbose:
                            print('   [{time}] complete in {x} sec'.format(time=t.ctime(),
                                                                           x=t.process_time() - timer_cse))

                        # emit declarations for CSE variables
                        _write_cse_decls(decls, h, nF, nP)

                        for i in range(nF):
                            for j in range(nF):
                                for k in range(nF):
                                    for l in range(nF):
                                        for m in range(nF):
                                            # emit main expression
                                            emit_expr = sym.printing.cxxcode(new_expr[
                                                                                 (nF) * (nF) * (nF) * (nF) * i + (
                                                                                     nF) * (
                                                                                     nF) * (
                                                                                     nF) * j + (nF) * (nF) * k + (
                                                                                     nF) * l + m])
                                            rw_expr = _rewrite_indices(emit_expr, nF, nP)
                                            h.write('\n RMcd[' + str(
                                                (nF) * (nF) * (nF) * (nF) * i + (nF) * (nF) * (nF) * j + (nF) * (
                                                    nF) * k + (
                                                    nF) * l + m) + ']=' + str(rw_expr) + ';\n')

    def _translate_pot(self, simplify_potential: bool, simplify_metric: bool, simplify_covd: bool, verbose: bool,
                       recache: bool):
        """
        Translate potential
        """

        nF, nP = self.nF, self.nP

        covd_sym = sym_tools.CovDSym(nF, nP, self.G, self.V, recache, simplify_fmet=simplify_metric,
                                     simplify_pot=simplify_potential,
                                     simplify=simplify_covd)

        v, vd, vdd, vddd = covd_sym.get_potential_sym_arrays()

        pot_template_path = os.path.join(cwd, 'cppt', 'potentialProto.h')
        pot_translated_path = os.path.join(cwd, 'cppt', 'potential.h')

        with open(pot_template_path, "r") as _:
            f = _.readlines()

        if verbose:
            timer = t.process_time()
            print('   [{time}] writing to potential.h'.format(time=t.ctime()))

        with open(pot_translated_path, "w") as g:

            for line in f:

                g.write(line)

                if line == "// #Rewrite\n":
                    g.write('// Potential file rewriten at' + ' ' + t.strftime("%c") + '\n')

                if line == "// #FP\n":
                    g.write('nF=' + str(nF) + ';\n' + 'nP=' + str(nP) + ';\n')

                if line == "// Pot\n":

                    # extract common subexpressions from V
                    if verbose:
                        timer_cse = t.process_time()
                        print('   [{time}] performing CSE for V'.format(time=t.ctime()))

                    decls, new_expr = sym.cse(v, order='none')
                    if verbose:
                        print('   [{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer_cse))

                    # emit declarations for CSE variables
                    _write_cse_decls(decls, g, nF, nP)

                    # emit main expression
                    emit_expr = sym.printing.cxxcode(new_expr[0])
                    rw_expr = _rewrite_indices(emit_expr, nF, nP)
                    g.write('  sum=' + str(rw_expr) + ';\n')

                if line == "// dPot\n":
                    if verbose:
                        timer_cse = t.process_time()
                        print('   [{time}] performing CSE for dV'.format(time=t.ctime()))

                    decls, new_exprs = sym.cse(vd, order='none')
                    if verbose:
                        print('   [{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer_cse))

                    # emit declarations for CSE variables
                    _write_cse_decls(decls, g, nF, nP)

                    for i in range(nF):
                        emit_expr = sym.printing.cxxcode(new_exprs[i])
                        rw_expr = _rewrite_indices(emit_expr, nF, nP)
                        g.write('\n sum[' + str(i) + ']=' + str(rw_expr) + ';\n')

                if line == "// ddPot\n":

                    if verbose:
                        timer_cse = t.process_time()
                        print('   [{time}] performing CSE for ddV'.format(time=t.ctime()))

                    decls, new_exprs = sym.cse(vdd, order='none')
                    if verbose:
                        print('   [{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer_cse))

                    # emit declarations for CSE variables
                    _write_cse_decls(decls, g, nF, nP)

                    for i in range(nF):
                        for j in range(nF):
                            emit_expr = sym.printing.cxxcode(new_exprs[i + nF * j])
                            rw_expr = _rewrite_indices(emit_expr, nF, nP)
                            g.write('\n sum[' + str(i + nF * j) + ']=' + str(rw_expr) + ';\n')

                if line == "// dddPot\n":

                    if verbose:
                        timer_cse = t.process_time()
                        print('   [{time}] performing CSE for dddV'.format(time=t.ctime()))

                    decls, new_exprs = sym.cse(vddd, order='none')
                    if verbose:
                        print('   [{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer_cse))

                    # emit declarations for CSE variables
                    _write_cse_decls(decls, g, nF, nP)

                    for i in range(nF):
                        for j in range(nF):
                            for k in range(nF):
                                emit_expr = sym.printing.cxxcode(new_exprs[i + nF * j + nF * nF * k])
                                rw_expr = _rewrite_indices(emit_expr, nF, nP)
                                g.write('\n sum[' + str(i + nF * j + nF * nF * k) + ']=' + str(rw_expr) + ';\n')

            if verbose:
                print('   [{time}] complete in {x} sec'.format(time=t.ctime(), x=t.process_time() - timer))

    @staticmethod
    def install(name: str):
        """
        Compiles c++ source code into python module. Importable from pyts.models with module name 'name'

        :param name: Python module name
        """

        t_start = t.ctime()
        print('   [{time}] start'.format(time=t_start))
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

        os.system("export CFLAGS='-I {}'".format(np.get_include()))
        subprocess.run(["bash", "moduleSetup.sh"], cwd=location)

        if os.path.exists(os.path.join(location, "build")):
            subprocess.run(["rm", "-r", "build"], cwd=location)

        if os.path.exists(os.path.join(location, "dist")):
            subprocess.run(["rm", "-r", "dist"], cwd=location)

        if os.path.exists(os.path.join(location, f"{name}.egg-info")):
            subprocess.run(["rm", "-r", f"{name}.egg-info"], cwd=location)

        models_root = os.path.abspath(os.path.join(cwd, "models"))

        models_tree = [x[0] for x in os.walk(models_root)]

        models_eggs = [x for x in models_tree if x.endswith(".egg") and name in x]

        for egg_path in models_eggs:

            fnames = [f for f in os.listdir(egg_path) if f.endswith(".so") or f.endswith(".py")]

            for fn in fnames:

                src = os.path.join(egg_path, fn)
                dst = os.path.join(models_root, fn)

                if os.path.exists(dst):
                    os.remove(dst)

                shutil.copyfile(src, dst)

        t_end = t.ctime()
        print("\n-- Compiled source in {} seconds, total time {} seconds".format(_delta_ctime(t_start_compile, t_end),
                                                                                 _delta_ctime(t_start, t_end)))
