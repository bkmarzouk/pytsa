import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import kde


def plot_kde(
    df: pd.DataFrame,
    key_x: str,
    key_y: str,
    pct_x=None,
    pct_y=None,
    x_cuts=None,
    y_cuts=None,
    nbins=500,
    ax=None,
    cmap="Spectral",
):
    """
    Plot 2d Kernel density estimates from columned pandas data

    :param df: Pandas data frame
    :param key_x: Column key for x plotting
    :param key_y: Column key for y plotting
    :param pct_x: Percentile cut on x data
    :param pct_y: Percentile cut on y data
    :param x_cuts: Domain cut on x data
    :param y_cuts: Domain cut on y data
    :param nbins: Number of bins in kernel estimate
    :param ax: Matplotlib axis object: If None,
               automatically renders figure in current axis
    """
    # Extract x and y
    x = df.dropna()[key_x]
    y = df.dropna()[key_y]

    locs = None

    if x_cuts is not None and y_cuts is not None:
        locs = (y > y_cuts[0]) * (y < y_cuts[1])
        locs *= (x > x_cuts[0]) * (x < x_cuts[1])
    elif x_cuts is not None:
        locs = (x > x_cuts[0]) * (x < x_cuts[1])
    elif y_cuts is not None:
        locs = (y > y_cuts[0]) * (y < y_cuts[1])
    else:
        pass

    if pct_x is not None and pct_y is not None:
        locs = (y > np.percentile(y, 50 - pct_y / 2)) * (
            y < np.percentile(y, 50 + pct_y / 2)
        )
        locs *= (x > np.percentile(x, 50 - pct_x / 2)) * (
            x < np.percentile(x, 50 + pct_x / 2)
        )
    elif pct_x is not None:
        locs = (x > np.percentile(x, 50 - pct_x / 2)) * (
            x < np.percentile(x, 50 + pct_x / 2)
        )
    elif pct_y is not None:
        locs = (y > np.percentile(y, 50 - pct_y / 2)) * (
            y < np.percentile(y, 50 + pct_y / 2)
        )
    else:
        pass

    x = x[locs] if locs is not None else x
    y = y[locs] if locs is not None else y

    # Define the borders
    deltaX = (max(x) - min(x)) / 10
    deltaY = (max(y) - min(y)) / 10
    xmin = min(x) - deltaX
    xmax = max(x) + deltaX
    ymin = min(y) - deltaY
    ymax = max(y) + deltaY

    # Create meshgrid
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]

    # Evaluate a gaussian kde on a regular grid of nbins x
    # nbins over data extents
    k = kde.gaussian_kde([x, y])
    xi, yi = np.mgrid[
        x.min() : x.max() : nbins * 1j, y.min() : y.max() : nbins * 1j
    ]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    # Make the plot
    method = plt if ax is None else ax
    method.pcolormesh(xi, yi, zi.reshape(xi.shape), shading="auto", cmap=cmap)
