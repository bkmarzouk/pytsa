import numpy as np
from pytransport.sampler.configs.setup_sampler import build_catalogue, Setup  # Import Setup for dill handling
from pytransport.sampler.methods import pyt_methods


class _2pt:

    def __int__(self):
        self.ns = None
        self.running = None
        self.error_code = None

    def add_result(self, ns, running, error_code):
        assert self.ns is None, self.ns
        self.ns = ns
        assert self.running is None, self.running
        self.running = running
        self.error_code = error_code


class _3pt:

    def __init__(self, alpha: int or float, beta: int or float):
        self.alpha = alpha
        self.beta = beta
        self.result = None
        self.error_code = None

    def add_result(self, result, error_code):
        assert self.result is None, self.result
        self.result = result
        self.error_code = error_code


class _Sample:
    def __int__(self, index, cache_loc):

        self.index = index
        self.cache_loc = cache_loc


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


def main(pool, setup: Setup, args):
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
