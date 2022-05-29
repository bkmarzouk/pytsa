import numpy as np
import scipy.interpolate.interpolate
from scipy.interpolate import UnivariateSpline, interp1d

k_order = 1
use_interp1d = True
s_method = interp1d if use_interp1d else UnivariateSpline

if use_interp1d:
    k_order = {
        1: 'linear',
        2: 'quadratic',
        3: 'cubic'
    }[k_order]


def update_default_k(k):
    """
    Updates global spline orders for fitting function (default values)

    :param k: spline order (see e.g. UnivariateSpline)
    """
    global k_order
    k_order = k


def monotonic(array: np.ndarray, col_idx=0, strict=True):
    """
    Checks that column (col_idx) of a 2D array is strictly monotonically increasing

    :param array: 2d numpy array
    :param strict: if True, strictly monotonic
    :return: True if monotonic, else false
    """

    col = array.T[col_idx]

    if strict:
        return np.all(np.concatenate([col[:-1] < np.roll(col, -1)[:-1], np.array([col[-2] < col[-1]])]))

    return np.all(np.concatenate([col[:-1] <= np.roll(col, -1)[:-1], np.array([col[-2] <= col[-1]])]))


def approx_row_closest(pivot_value, array: np.ndarray, col_idx):
    """
    Finds the closest data point in the array to a target pivot value, and returns the containing row data

    :param pivot_value: value to approximate array at
    :param array: input array
    :param col_idx: column containing pivot value
    :return: row from array approximately matching the pivot value
    """

    lower = None
    upper = None
    ii = None

    compare_to = array.T[col_idx]

    for ii in range(1, len(array)):

        if compare_to[ii] >= pivot_value > compare_to[ii - 1]:
            lower = compare_to[ii - 1]
            upper = compare_to[ii]
            break

    if lower is None:
        print("COULD NOT FIND SUFFICIENT VALUE")

    if abs(lower - pivot_value) < abs(upper - pivot_value):
        return array[ii - 1]
    else:
        return array[ii]

    idx_h = None
    val_h = None

    idx_l = None
    val_l = None

    for idx, v in enumerate(array.T[col_idx]):
        if v > pivot_value:
            idx_h = idx
            val_h = v
            break

    backwards = array.T[col_idx][::-1]

    for idx, v in enumerate(backwards):
        if v < pivot_value:
            idx_l = -(idx + 1)
            val_l = v
            break

    if abs(val_l - pivot_value) < abs(val_h - pivot_value):
        approx_idx = idx_l
    else:
        approx_idx = idx_h

    return array[approx_idx]


def approx_row_spline(pivot_value, pivot_width, array: np.ndarray, col_idx):
    """
    Finds the closest data point in the array to a target pivot value, and returns the containing row data

    :param pivot_value: value to approximate array at
    :param array: input array
    :param col_idx: column containing pivot value
    :return: row from array approximately matching the pivot value
    """

    assert monotonic(array, col_idx), "Can only formulate approximate row data if independent col is strictly " \
                                      "monotonically increasing!"

    pivot_width /= 2

    l_idx = np.where(array.T[col_idx] > pivot_value - pivot_width)
    h_idx = np.where(array.T[col_idx] < pivot_value + pivot_width)
    indices = np.intersect1d(l_idx, h_idx)

    window_array = array.copy()[indices]

    ind_data = window_array.T[col_idx]
    dep_data = [window_array.T[ii] for ii in range(window_array.shape[1]) if ii != col_idx]

    if use_interp1d:
        spl_data = [s_method(ind_data, _dep, kind=k_order) for _dep in dep_data]
    else:
        spl_data = [s_method(ind_data, _dep, k=k_order) for _dep in dep_data]

    out = np.array([s(pivot_value) for s in spl_data])
    out = np.insert(out, col_idx, pivot_value)

    return out
