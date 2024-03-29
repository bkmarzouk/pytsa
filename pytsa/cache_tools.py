from hashlib import sha224, md5
import sympy as sym
import os
import pickle as pk
import numpy as np

sym_cache = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "sym_cache")
)
model_cache = os.path.join(sym_cache, "model")
fmet_cache = os.path.join(sym_cache, "fmet")
pot_cache = os.path.join(sym_cache, "pot")
covd_cache = os.path.join(sym_cache, "covd")


def hash_alpha_beta(alpha: float or int, beta: float or int):
    arr = np.array([alpha, beta], dtype=np.float64)
    return md5(np.copy(arr, order="C")).hexdigest()


def hash_pars(*pars: object):
    s = ""

    for par in pars:
        s += str(par)

    # Cheap and nasty way of converting arbitrary parameters into an
    # identifying string

    return sha224(bytes(s, encoding="utf-8")).hexdigest()


def fmet_path(metric: sym.Matrix, simplify=True):
    """
    Creates hash string from symbolic metric and builds path

    :param metric: sympy representation of metric
    :param simplify: if True, simplified representation
    :return: path to metric calculations
    """

    assert isinstance(
        metric, sym.Matrix
    ), "Must pass sympy matrix representation of metric: {}".format(metric)
    assert metric.is_symmetric(), "Metric must be symmetric: {}".format(metric)

    expr_str = str(metric)
    if simplify:
        expr_str += "_simple"

    hash_str = sha224(bytes(expr_str, encoding="utf-8")).hexdigest()

    fname = os.path.join(fmet_cache, hash_str + ".fmet")

    return fname


def pot_path(potential: sym.Expr, simplify=True):
    """
    Creates hash string from symbolic potential and builds path

    :param potential: sympy expression for potential
    :param simplify: if True, simplified potential
    :return: path to potential calculations
    """
    assert isinstance(
        potential, sym.Expr
    ), "Must pass sympy expression for the potential: {}".format(potential)

    expr_str = str(potential)

    if simplify:
        expr_str += "_simple"

    hash_str = sha224(bytes(expr_str, encoding="utf-8")).hexdigest()

    fname = os.path.join(pot_cache, hash_str + ".pot")

    return fname


def covd_path(
    metric: sym.Matrix,
    potential: sym.Expr,
    simplify=True,
    simplify_fmet=True,
    simplify_pot=True,
):
    """
    Creates hash string from sympy metric and potential and builds path

    :param metric: sympy representation of metric
    :param potential: sympy expression for potential
    :param simplify: if True, simplified covd results
    :param simplify_fmet: if True, simplified metric
    :param simplify_pot: if True, simplified potential
    :return: path to joint metric + potential calculations
    """

    fpath = fmet_path(metric, simplify=simplify_fmet)
    vpath = pot_path(potential, simplify=simplify_pot)

    hash_str_fmet = os.path.split(fpath)[-1].split(".")[0]
    hash_str_pot = os.path.split(vpath)[-1].split(".")[0]

    if simplify:
        hash_str_pot += "_simple"

    fname = os.path.join(
        covd_cache, hash_str_fmet + "_" + hash_str_pot + ".covd"
    )

    return fname


def cache_fmet(metric: sym.Matrix, fmet_derived, recache=False, simplify=True):
    path = fmet_path(metric, simplify=simplify)

    if os.path.exists(path):
        if recache:
            os.remove(path)
        else:
            raise OSError(
                "Metric calculations already cached: {}".format(path)
            )

    with open(path, "wb") as f:
        pk.dump(fmet_derived, f)


def load_fmet(metric: sym.Matrix, delete=False, simplify=True):
    path = fmet_path(metric, simplify=simplify)

    if delete and os.path.exists(path):
        os.remove(path)

    if os.path.exists(path):
        with open(path, "rb") as f:
            fmet_derived = pk.load(f)
    else:
        raise OSError(
            "Field metric derived object does not exist: {}".format(path)
        )

    return fmet_derived


def cache_pot(potential: sym.Expr, pot_derived, recache=False, simplify=True):
    path = pot_path(potential, simplify=simplify)

    if os.path.exists(path):
        if recache:
            os.remove(path)
        else:
            raise OSError(
                "Potential calculation already cached: {}".format(path)
            )

    with open(path, "wb") as f:
        pk.dump(pot_derived, f)


def load_pot(potential: sym.Expr, delete=False, simplify=True):
    path = pot_path(potential, simplify=simplify)

    if delete and os.path.exists(path):
        os.remove(path)

    if os.path.exists(path):
        with open(path, "rb") as f:
            pot_derived = pk.load(f)
    else:
        raise OSError(
            "Potential derived object does not exist: {}".format(path)
        )

    return pot_derived


def cache_covd(
    metric: sym.Matrix,
    potential: sym.Expr,
    covd_derived,
    recache=False,
    simplify=True,
    simplify_fmet=True,
    simplify_pot=True,
):
    path = covd_path(
        metric,
        potential,
        simplify=simplify,
        simplify_pot=simplify_pot,
        simplify_fmet=simplify_fmet,
    )

    if os.path.exists(path):
        if recache:
            os.remove(path)
        else:
            raise OSError(
                "Covariant derivatives already cached: {}".format(path)
            )

    with open(path, "wb") as f:
        pk.dump(covd_derived, f)


def load_covd(
    metric: sym.Matrix,
    potential: sym.Expr,
    delete=False,
    simplify_fmet=True,
    simplify_pot=True,
    simplify=True,
):
    path = covd_path(
        metric,
        potential,
        simplify=simplify,
        simplify_fmet=simplify_fmet,
        simplify_pot=simplify_pot,
    )

    if delete and os.path.exists(path):
        os.remove(path)

    if os.path.exists(path):
        with open(path, "rb") as f:
            covd_derived = pk.load(f)
    else:
        raise OSError("Covariant derivatives do not exists: {}".format(path))

    return covd_derived
