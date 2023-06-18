import os.path
import numpy as np
import pickle as pk
from numpy.random import bit_generator, Generator


class RandomStates(object):
    def __init__(self, n_states, entropy=None):
        """
        Spawns sequence of random states from an initial seeding value
        :param entropy: initial seed value, if None assign from system entropy
        :param n_states: number of child seeds to spawn
        """
        seq = bit_generator.SeedSequence(entropy=entropy)

        self.entropy = (
            seq.entropy
        )  # We log the entropy after seq init in case input entropy is None
        self.n_states = n_states  # Number of states to spawn
        self.states = seq.spawn(n_states)  # builds new uncorrelated streams

    def get_state(self, idx: int) -> Generator:
        """
        Gets random number generator state
        :param idx: index of state
        :return: Generator instance
        """
        assert idx < self.n_states, [idx, self.n_states]
        return np.random.default_rng(self.states[idx])

    @classmethod
    def from_cache(cls, cache_loc: str, n_states: int, entropy=None):
        """
        Attempts loading states from cache. If not found, builds cache file

        :param path: Path to random states object
        :param n_states: number of states to spawn
        :param entropy: system entropy to initialize sequence with;
                should be None unless for test cases
        :return: RandomStates object
        """

        path = os.path.join(cache_loc, "rng_states.pk")

        if os.path.exists(path):
            with open(path, "rb") as f:
                s: RandomStates = pk.load(f)

                # check states match
                assert s.n_states == n_states, [s.n_states, n_states]

                if entropy is not None:
                    assert s.entropy == entropy, [s.entropy, entropy]

        new_state = cls(n_states, entropy=entropy)

        with open(path, "wb") as f:
            pk.dump(new_state, f)

        return new_state
