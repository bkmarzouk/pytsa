from hashlib import sha224
import sympy as sym
import os
import pickle as pk

model_cache = os.environ["model_cache"]
fmet_cache = os.environ["fmet_cache"]
pot_cache = os.environ["pot_cache"]
covd_cache = os.environ["covd_cache"]


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


def covd_path(metric: sym.Matrix, potential: sym.Expr):
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

    fname = os.path.join(covd_cache, hash_str_fmet + "_" + hash_str_pot + ".covd")

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


def cache_covd(metric: sym.Matrix, potential: sym.Expr, covd_derived, recache=False):
    path = covd_path(metric, potential)

    if os.path.exists(path):
        if recache:
            os.remove(path)
        else:
            raise OSError("Covariant derivatives already cached: {}".format(path))

    with open(path, "wb") as f:
        pk.dump(covd_derived, f)


def load_covd(metric: sym.Matrix, potential: sym.Expr, delete=False):
    path = covd_path(metric, potential)

    if delete and os.path.exists(path):
        os.remove(path)

    if os.path.exists(path):
        with open(path, "rb") as f:
            covd_derived = pk.load(f)
    else:
        raise OSError("Covariant derivatives do not exists: {}".format(path))

    return covd_derived
