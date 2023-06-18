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
# 3rd in {} defines the input for the method;
# e.g. METH_NOARGS -> no arguments expected for function
# 4th in {} defines the doc string, user can then yield info using
# help(<module>.method). docs define in PyTrans.cpp
_PyMethodDefs = ",".join(
    [
        '{"H", (PyCFunction)MT_H,   METH_VARARGS, doc_H}',
        '{"nF", (PyCFunction)MT_fieldNumber,   METH_NOARGS, doc_nF}',
        '{"nP", (PyCFunction)MT_paramNumber,   METH_NOARGS, doc_nP}',
        '{"Epsilon", (PyCFunction)MT_Ep,   METH_VARARGS, doc_Epsilon}',
        '{"Eta", (PyCFunction)MT_Eta,   METH_VARARGS, doc_Eta}',
        '{"V", (PyCFunction)MT_V,   METH_VARARGS, doc_V}',
        '{"dV", (PyCFunction)MT_dV,   METH_VARARGS, doc_dV}',
        '{"ddV", (PyCFunction)MT_ddV,   METH_VARARGS, doc_ddV}',
        '{"massMatrix", (PyCFunction)MT_massMatrix,   '
        "METH_VARARGS, doc_massMatrix}",
        '{"findEndOfInflation", (PyCFunction)MT_findEndOfInflation,   '
        "METH_VARARGS, doc_findEndOfInflation}",
        '{"backEvolve", (PyCFunction)MT_backEvolve,   '
        "METH_VARARGS, doc_backEvolve}",
        '{"sigEvolve", (PyCFunction)MT_sigEvolve,   '
        "METH_VARARGS, doc_sigEvolve}",
        '{"alphaEvolve", (PyCFunction)MT_alphaEvolve,   '
        "METH_VARARGS, doc_alphaEvolve}",
    ]
)

cwd = os.path.abspath(os.path.dirname(__file__))

# Installation location for models
MODELS_LOC = os.path.join(cwd, "models")

_pytrans_path = os.path.join(cwd, "pyt", "_PyTrans.cpp")
pytrans_path = os.path.join(cwd, "pyt", "PyTrans.cpp")

cppt_dir = os.path.join(cwd, "cppt")
stepper_path = os.path.join(cppt_dir, "stepper", "rkf45.hpp")


def _delta_ctime(a, b):
    """
    Translates time.ctime() outputs into time difference in seconds

    :param a: inital time
    :param b: final time
    :return: seconds difference
    """

    def hours(t):
        return int(t.split(":")[0][-2:])

    def mins(t):
        return int(t.split(":")[1])

    def secs(t):
        return int(t.split(":")[2][:2])

    i = hours(a) * 60 * 60 + mins(a) * 60 + secs(a)
    j = hours(b) * 60 * 60 + mins(b) * 60 + secs(b)

    return j - i


def _set_template_headers(non_canonical: bool):
    """
    Rewrites header files used in compile script to point to
    (non-)canonical file versions

    :param non_canonical: if True methods use non-canonical
    """

    with open(_pytrans_path, "r") as f:
        lines = f.readlines()

    src_dir = os.path.join(cppt_dir, "NC") if non_canonical else cppt_dir
    step_dir = os.path.join(cppt_dir, "stepper")

    with open(pytrans_path, "w") as f:
        for cpp_line in lines:
            if "CPPT_DIR" in cpp_line:
                f.write(cpp_line.replace("CPPT_DIR", src_dir))
            elif "STEP_DIR" in cpp_line:
                f.write(cpp_line.replace("STEP_DIR", step_dir))
            else:
                f.write(cpp_line)


def _rewrite_indices(expr, nF, nP):
    """
    Rewrite symbolic expressions containing subsripts as
    indexed vector components
    e.g. f_0 -> f[0]

    :param expr: symbolic expression in string format
    :param nF: number of fields
    :param nP: number of params
    :return: translated expression
    """
    new_expr = expr

    for ll in range(max(nP, nF)):
        ll = max(nP, nF) - 1 - ll
        new_expr = new_expr.replace("_" + str(ll), "[" + str(ll) + "]")

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
        g.write("  auto " + symb + " = " + new_expr + ";\n")


def print_compute_time(msg, ref=None):
    time = t.ctime()
    if ref is None:
        print(f"[{time}] {msg}")
    else:
        secs = t.time() - ref
        print(f"[{time}] {msg} in {secs} seconds")


class Translator:
    def __init__(
        self,
        nF: int,
        nP: int,
        V: sym.Expr,
        G: sym.Matrix or None = None,
        simplify_potential=True,
        simplify_metric=True,
        simplify_covd=True,
        verbose=True,
        recache=False,
    ):
        """
        Translates SymPy expressions into c++ source ready for compilation

        :param nF: number of fields
        :param nP: number of params
        :param V: sym expression for potential
        :param G: sym expression for fieldspace metric,
                    may be None then euclidean metric assumed
        :param simplify_potential: if True,
                    simplifies potential and partial derivatives of
        :param simplify_metric: if True, simplifies metric and inv. metric
        :param simplify_covd: if True,
                                simplifies covariant derivative expressions
        :param verbose: if True, prints status as output
        :param recache: if True, rebuilds symbolic expressions,
                        otherwise uses cached results
        """
        self.nF = nF
        self.nP = nP
        self.V = V
        self.G = G
        self.canonical = G is None

        self._translate_fmet(
            simplify_potential,
            simplify_metric,
            simplify_covd,
            verbose,
            recache,
        )
        self._translate_pot(
            simplify_potential,
            simplify_metric,
            simplify_covd,
            verbose,
            recache,
        )
        _set_template_headers(G is not None)

    def _translate_fmet(
        self,
        simplify_potential: bool,
        simplify_metric: bool,
        simplify_covd: bool,
        verbose: bool,
        recache: bool,
    ):
        """
        Translate field space components
        """

        nF, nP = self.nF, self.nP

        covd_sym = sym_tools.CovDSym(
            nF,
            nP,
            self.G,
            self.V,
            recache,
            simplify_fmet=simplify_metric,
            simplify_pot=simplify_potential,
            simplify=simplify_covd,
        )

        if not self.canonical:
            fmet_template_path = os.path.join(
                cwd, "cppt", "fieldmetricProto.h"
            )
            fmet_translated_path = os.path.join(cwd, "cppt", "fieldmetric.h")

            with open(fmet_template_path, "r") as f:
                template_lines = f.readlines()

            if verbose:
                print("-- Translating field space symbols")

            with open(fmet_translated_path, "w") as h:
                (
                    G_array,
                    Gamma_array,
                    R_array,
                    gradR_array,
                ) = covd_sym.get_curvature_sym_arrays()

                for line in template_lines:
                    h.write(line)
                    if line == "// #FP\n":
                        h.write(
                            "nF=" + str(nF) + ";\n" + "nP=" + str(nP) + ";\n"
                        )

                    if line == "// metric\n":
                        timer_cse = t.time()
                        if verbose:
                            print_compute_time(
                                "Performing CSE for field metric"
                            )

                        decls, new_expr = sym.cse(G_array, order="none")

                        if verbose:
                            print_compute_time(
                                "CSE for field metric",
                                ref=timer_cse,
                            )

                        # emit declarations for CSE variables
                        _write_cse_decls(decls, h, nF, nP)

                        for ii in range(2 * nF):
                            for jj in range(2 * nF):
                                # emit main expression
                                emit_expr = sym.printing.cxxcode(
                                    new_expr[(2 * nF) * ii + jj]
                                )
                                rw_expr = _rewrite_indices(emit_expr, nF, nP)
                                h.write(
                                    "\n FM["
                                    + str((2 * nF) * ii + jj)
                                    + "]="
                                    + str(rw_expr)
                                    + ";\n"
                                )

                    if line == "// Christoffel\n":
                        timer_cse = t.process_time()

                        if verbose:
                            print_compute_time(
                                "Performing CSE for Christoffel symbols"
                            )

                        decls, new_expr = sym.cse(Gamma_array, order="none")

                        if verbose:
                            print_compute_time(
                                "CSE for Christoffel symbols", ref=timer_cse
                            )

                        # emit declarations for CSE variables
                        _write_cse_decls(decls, h, nF, nP)

                        for ii in range(2 * nF):
                            for jj in range(2 * nF):
                                for kk in range(2 * nF):
                                    # emit main expression
                                    emit_expr = sym.printing.cxxcode(
                                        new_expr[
                                            (2 * nF) * (2 * nF) * ii
                                            + (2 * nF) * jj
                                            + kk
                                        ]
                                    )
                                    rw_expr = _rewrite_indices(
                                        emit_expr, nF, nP
                                    )
                                    h.write(
                                        "\n CS["
                                        + str(
                                            (2 * nF) * (2 * nF) * ii
                                            + (2 * nF) * jj
                                            + kk
                                        )
                                        + "]="
                                        + str(rw_expr)
                                        + ";\n"
                                    )

                    if line == "// Riemann\n":
                        timer_cse = t.process_time()

                        if verbose:
                            print_compute_time(
                                "Performing CSE for Riemann tensor"
                            )

                        decls, new_expr = sym.cse(R_array, order="none")
                        if verbose:
                            print_compute_time(
                                " CSE for Riemann tensor",
                                timer_cse,
                            )

                        # emit declarations for CSE variables
                        _write_cse_decls(decls, h, nF, nP)

                        for ii in range(nF):
                            for jj in range(nF):
                                for kk in range(nF):
                                    for ll in range(nF):
                                        # emit main expression
                                        emit_expr = sym.printing.cxxcode(
                                            new_expr[
                                                (nF) * (nF) * (nF) * ii
                                                + (nF) * (nF) * jj
                                                + (nF) * kk
                                                + ll
                                            ]
                                        )
                                        rw_expr = _rewrite_indices(
                                            emit_expr, nF, nP
                                        )
                                        h.write(
                                            "\n RM["
                                            + str(
                                                (nF) * (nF) * (nF) * ii
                                                + (nF) * (nF) * jj
                                                + (nF) * kk
                                                + ll
                                            )
                                            + "]="
                                            + str(rw_expr)
                                            + ";\n"
                                        )

                    if line == "// Riemanncd\n":
                        timer_cse = t.process_time()

                        if verbose:
                            print_compute_time(
                                "Performing CSE for grad of Riemann tensor"
                            )

                        decls, new_expr = sym.cse(gradR_array, order="none")

                        if verbose:
                            print_compute_time(
                                "CSE for grad of Riemann tensor", ref=timer_cse
                            )

                        # emit declarations for CSE variables
                        _write_cse_decls(decls, h, nF, nP)

                        for ii in range(nF):
                            for jj in range(nF):
                                for kk in range(nF):
                                    for ll in range(nF):
                                        for mm in range(nF):
                                            # emit main expression
                                            emit_expr = sym.printing.cxxcode(
                                                new_expr[
                                                    (nF)
                                                    * (nF)
                                                    * (nF)
                                                    * (nF)
                                                    * ii
                                                    + (nF) * (nF) * (nF) * jj
                                                    + (nF) * (nF) * kk
                                                    + (nF) * ll
                                                    + mm
                                                ]
                                            )
                                            rw_expr = _rewrite_indices(
                                                emit_expr, nF, nP
                                            )
                                            h.write(
                                                "\n RMcd["
                                                + str(
                                                    (nF)
                                                    * (nF)
                                                    * (nF)
                                                    * (nF)
                                                    * ii
                                                    + (nF) * (nF) * (nF) * jj
                                                    + (nF) * (nF) * kk
                                                    + (nF) * ll
                                                    + mm
                                                )
                                                + "]="
                                                + str(rw_expr)
                                                + ";\n"
                                            )

    def _translate_pot(
        self,
        simplify_potential: bool,
        simplify_metric: bool,
        simplify_covd: bool,
        verbose: bool,
        recache: bool,
    ):
        """
        Translate potential
        """

        nF, nP = self.nF, self.nP

        covd_sym = sym_tools.CovDSym(
            nF,
            nP,
            self.G,
            self.V,
            recache,
            simplify_fmet=simplify_metric,
            simplify_pot=simplify_potential,
            simplify=simplify_covd,
        )

        v, vd, vdd, vddd = covd_sym.get_potential_sym_arrays()

        pot_template_path = os.path.join(cwd, "cppt", "potentialProto.h")
        pot_translated_path = os.path.join(cwd, "cppt", "potential.h")

        with open(pot_template_path, "r") as _:
            f = _.readlines()

        timer_translate = t.process_time()

        if verbose:
            print_compute_time("Writing potential.h")

        with open(pot_translated_path, "w") as g:
            for line in f:
                g.write(line)

                if line == "// #Rewrite\n":
                    g.write(
                        "// Potential file rewriten at"
                        + " "
                        + t.strftime("%c")
                        + "\n"
                    )

                if line == "// #FP\n":
                    g.write("nF=" + str(nF) + ";\n" + "nP=" + str(nP) + ";\n")

                if line == "// Pot\n":
                    timer_cse = t.process_time()

                    # extract common subexpressions from V
                    if verbose:
                        print_compute_time("Performing CSE for V")

                    decls, new_expr = sym.cse(v, order="none")

                    if verbose:
                        print_compute_time("CSE for V", timer_cse)

                    # emit declarations for CSE variables
                    _write_cse_decls(decls, g, nF, nP)

                    # emit main expression
                    emit_expr = sym.printing.cxxcode(new_expr[0])
                    rw_expr = _rewrite_indices(emit_expr, nF, nP)
                    g.write("  sum=" + str(rw_expr) + ";\n")

                if line == "// dPot\n":
                    timer_cse = t.process_time()

                    if verbose:
                        print_compute_time("Performing CSE for dV")

                    decls, new_exprs = sym.cse(vd, order="none")

                    if verbose:
                        print_compute_time("CSE for dV", ref=timer_cse)

                    # emit declarations for CSE variables
                    _write_cse_decls(decls, g, nF, nP)

                    for i in range(nF):
                        emit_expr = sym.printing.cxxcode(new_exprs[i])
                        rw_expr = _rewrite_indices(emit_expr, nF, nP)
                        g.write(
                            "\n sum[" + str(i) + "]=" + str(rw_expr) + ";\n"
                        )

                if line == "// ddPot\n":
                    timer_cse = t.process_time()

                    if verbose:
                        print_compute_time("Performing CSE for ddV")

                    decls, new_exprs = sym.cse(vdd, order="none")

                    if verbose:
                        print_compute_time("CSE for ddV", ref=timer_cse)

                    # emit declarations for CSE variables
                    _write_cse_decls(decls, g, nF, nP)

                    for i in range(nF):
                        for j in range(nF):
                            emit_expr = sym.printing.cxxcode(
                                new_exprs[i + nF * j]
                            )
                            rw_expr = _rewrite_indices(emit_expr, nF, nP)
                            g.write(
                                "\n sum["
                                + str(i + nF * j)
                                + "]="
                                + str(rw_expr)
                                + ";\n"
                            )

                if line == "// dddPot\n":
                    timer_cse = t.process_time()

                    if verbose:
                        print_compute_time("Performing CSE for dddV")

                    decls, new_exprs = sym.cse(vddd, order="none")

                    if verbose:
                        print_compute_time("CSE for dddV", ref=timer_cse)

                    # emit declarations for CSE variables
                    _write_cse_decls(decls, g, nF, nP)

                    for i in range(nF):
                        for j in range(nF):
                            for k in range(nF):
                                emit_expr = sym.printing.cxxcode(
                                    new_exprs[i + nF * j + nF * nF * k]
                                )
                                rw_expr = _rewrite_indices(emit_expr, nF, nP)
                                g.write(
                                    "\n sum["
                                    + str(i + nF * j + nF * nF * k)
                                    + "]="
                                    + str(rw_expr)
                                    + ";\n"
                                )

            if verbose:
                print(
                    "   [{time}] complete in {x} sec".format(
                        time=t.ctime(), x=t.process_time() - timer_translate
                    )
                )

    @staticmethod
    def install(name: str):
        """
        Compiles c++ source code into python module. Importable
        from pyts.models with module name 'name'

        :param name: Python module name
        """

        t_start = t.time()
        print_compute_time("Start model installation")
        cwd = os.path.dirname(__file__)
        location = os.path.join(cwd, "pyt")

        # Make copy of setup file that will install module
        _setup_file_path = os.path.join(cwd, "pyt", "_moduleSetup.py")
        setup_file_path = os.path.join(cwd, "pyt", "moduleSetup.py")
        shutil.copyfile(_setup_file_path, setup_file_path)

        with open(setup_file_path, "r") as f:
            module_setup_lines = f.readlines()

        # Rewrite module name in template file
        with open(setup_file_path, "w") as f:
            for line in module_setup_lines:
                if "PYT_MODNAME" in line:
                    f.write("mod_name = '{}'  # PYT_MODNAME\n".format(name))
                else:
                    f.write(line)

        # Make copy of cpp file that will be compiled
        assert isinstance(_pytrans_path, str) and os.path.exists(
            _pytrans_path
        ), _pytrans_path

        # Load PyTransport C++ file
        with open(pytrans_path, "r") as f:
            pytrans_cpp_lines = f.readlines()

        with open(pytrans_path, "w") as f:
            for line in pytrans_cpp_lines:
                if (
                    not line.endswith("//FuncDef\n")
                    and not line.endswith("//initFunc\n")
                    and not line.endswith("//modDef\n")
                ):
                    f.write(line)
                if line.endswith("//FuncDef\n"):
                    f.write(
                        "static PyMethodDef "
                        + name
                        + "_funcs[] = {"
                        + _PyMethodDefs
                        + ",   {NULL}};//FuncDef\n"
                    )

                if line.endswith("//modDef\n"):
                    f.write(
                        "static struct PyModuleDef PyTransModule = "
                        '{PyModuleDef_HEAD_INIT, "'
                        + name
                        + '", PyTrans_docs, -1, '
                        + name
                        + "_funcs}; //modDef\n"
                    )

                if line.endswith("//initFunc\n"):
                    f.write(
                        "PyMODINIT_FUNC PyInit_"
                        + name
                        + "(void)    {import_array();  "
                        "return PyModule_Create(&PyTransModule);} "
                        "//initFunc\n"
                    )

        t_start_compile = t.time()
        print_compute_time("Start compile phase")

        os.system("export CFLAGS='-I {}'".format(np.get_include()))
        subprocess.run(["bash", "moduleSetup.sh"], cwd=location)

        if os.path.exists(os.path.join(location, "build")):
            subprocess.run(["rm", "-r", "build"], cwd=location)

        if os.path.exists(os.path.join(location, "dist")):
            subprocess.run(["rm", "-r", "dist"], cwd=location)

        if os.path.exists(os.path.join(location, f"{name}.egg-info")):
            subprocess.run(["rm", "-r", f"{name}.egg-info"], cwd=location)

        # remove temporary file versions
        os.remove(setup_file_path)
        os.remove(pytrans_path)

        models_root = os.path.abspath(os.path.join(cwd, "models"))

        models_tree = [x[0] for x in os.walk(models_root)]

        models_eggs = [
            x for x in models_tree if x.endswith(".egg") and name in x
        ]

        for egg_path in models_eggs:
            fnames = [
                f
                for f in os.listdir(egg_path)
                if f.endswith(".so") or f.endswith(".py")
            ]

            for fn in fnames:
                src = os.path.join(egg_path, fn)
                dst = os.path.join(models_root, fn)

                if os.path.exists(dst):
                    os.remove(dst)

                shutil.copyfile(src, dst)

        print_compute_time("Compiled sources", ref=t_start_compile)

        print_compute_time("Model built in", ref=t_start)
