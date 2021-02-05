from hashlib import sha224
import numpy as np
import sympy as sym
import os
import pickle as pk

model_cache = os.environ["model_cache"]
fmet_cache = os.environ["fmet_cache"]
pot_cache = os.environ["pot_cache"]


def fmet_path(metric: sym.Matrix):
    """
    Creates hash string from symbolic metric and builds path

    :param metric: sympy representation of metric
    :return: path to metric calculations
    """
    assert isinstance(metric, sym.Matrix), "Must pass sympy matrix representation of metric: {}".format(metric)
    assert metric.is_symmetric(), "Metric must be symmetric: {}".format(metric)

    expr_str = str(metric)

    hash_str = sha224(bytes(expr_str, encoding='utf-8')).hexdigest()

    fname = os.path.join(fmet_cache, hash_str + ".fmet")

    return fname


def pot_path(potential: sym.Expr):
    """
    Creates hash string from symbolic potential and builds path

    :param potential: sympy expression for potential
    :return: path to potential calculations
    """
    assert isinstance(potential, sym.Expr), "Must pass sympy expression for the potential: {}".format(potential)

    expr_str = str(potential)

    hash_str = sha224(bytes(expr_str, encoding='utf-8')).hexdigest()

    fname = os.path.join(pot_cache, hash_str + ".pot")

    return fname


def model_path(metric: sym.Matrix, potential: sym.Expr):
    """
    Creates hash string from sympy metric and potential and builds path

    :param metric: sympy representation of metric
    :param potential: sympy expression for potential
    :return: path to joint metric + potential calculations
    """

    fpath = fmet_path(metric)
    vpath = pot_path(potential)

    hash_str_fmet = os.path.split(fpath)[-1].split(".")[0]
    hash_str_pot = os.path.split(vpath)[-1].split(".")[0]

    fname = os.path.join(model_cache, hash_str_fmet + "_" + hash_str_pot + ".model")

    return fname


def cache_fmet(metric: sym.Matrix, fmet_derived, recache=False):
    path = fmet_path(metric)

    if os.path.exists(path):
        if recache:
            os.remove(path)
        else:
            raise OSError("Metric caclulations already cached: {}".format(path))

    with open(path, "wb") as f:
        pk.dump(fmet_derived, f)


def load_fmet(metric: sym.Matrix, delete=False):
    path = fmet_path(metric)

    if delete and os.path.exists(path):
        os.remove(path)

    if os.path.exists(path):
        with open(path, "rb") as f:
            fmet_derived = pk.load(f)
    else:
        raise OSError("Field metric derived object does not exist: {}".format(path))

    return fmet_derived


def cache_pot(potential: sym.Expr, pot_derived, recache=False):
    path = pot_path(potential)

    if os.path.exists(path):
        if recache:
            os.remove(path)
        else:
            raise OSError("Potential calculation already cached: {}".format(path))

    with open(path, "wb") as f:
        pk.dump(pot_derived, f)


def load_pot(potential: sym.Expr, delete=False):
    path = pot_path(potential)

    if delete and os.path.exists(path):
        os.remove(path)

    if os.path.exists(path):
        with open(path, "rb") as f:
            pot_derived = pk.load(f)
    else:
        raise OSError("Potential derived object does not exist: {}".format(path))

    return pot_derived


def cache_model(metric: sym.Matrix, potential: sym.Expr, delete=False):
    path = model_path(metric, potential)

    if delete and os.path.exists(path):
        os.remove(path)

    if os.path.exists(path):
        with open(path, "rb") as f:
            model_derived = pk.load(f)
    else:
        raise OSError("Model derived from metric and potential does not exists: {}".format(path))

    return model_derived


nF = 2  # number of fields needed to define the double quadratic potential
nP = 3  # number of parameters needed to define the double quadtartic potential
f = sym.symarray('f', nF)  # an array representing the nF fields present for this model
p = sym.symarray('p',
                 nP)  # an array representing the nP parameters needed to define this model (that we might wish to change) if we don't
# wish to change them they could be typed explicitly in the potential below

V = 1. / 2. * p[0] ** 2.0 * f[0] ** 2.0 + 1. / 2. * p[1] ** 2.0 * f[
    1] ** 2.0  # this is the potential written in sympy notation
G = sym.Matrix([[p[2] ** 2.0, 0], [0, p[2] ** 2.0 * sym.sin(f[0]) ** 2.0]])

# fpath = fmet_path(G)
# vpath = pot_path(V)
# mpath = model_path(G, V)
#
# for p in [fpath, vpath, mpath]:
#     print(p)
