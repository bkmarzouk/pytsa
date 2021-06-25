import argparse
import os
import dill as pk
import numpy as np

from PyTransport.Sampler.configs.setup_sampler import build_catalogue, Setup  # Import Setup for dill handling
from PyTransport.Sampler.methods import pyt_methods


class _Data(object):

    def __init__(self):
        self.background = None
        self.Nend = None

    def _add_value(self, attr_name, value):
        assert not hasattr(self, attr_name), attr_name
        setattr(self, attr_name, value)

    def add_fields_dotfields(self, fields_dotfields):
        self._add_value("fields_dotfields", fields_dotfields)

    def add_params(self, params):
        self._add_value("params", params)

    def add_background(self, background):
        self._add_value("background", background)
        self._add_value("Nend", None if background is None else background[-1][0])


def main(pool, setup, args):
    n_samples = args.n_samples

    build_catalogue(setup, args)

    # Print configuration to command line
    print("-- Starting sampling routine: {}".format(setup.PyT))

    pool_data = []

    # Construct initial numpy array of samples (required to support latin sampling)
    catalogue_path = build_catalogue(setup, args, path_only=True)

    for idx in range(n_samples):
        pd = _Data()
        pd.PyT = setup.PyT.__name__
        pd.model_number = idx
        pd.path = catalogue_path
        pd.tols = setup.tols
        pd.N_min = setup.N_min
        pd.N_adi = setup.N_adiabitc
        pd.N_sub = setup.N_sub_evo
        pool_data.append(pd)

    # get results from background (returns data objs. with updated ics, params and back vals
    r = pool.map(pyt_methods.compute_background, pool_data)

    # TODO: cache and reload pool data to following list

    # redefine pool data with successful trajectories
    pool_data = [dat for dat in r if isinstance(dat.background, np.ndarray)]

