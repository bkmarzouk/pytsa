from numpy.random import PCG64, SeedSequence, default_rng
import pickle as pk
from dataclasses import dataclass


class _Generators:

    def __init__(self, n_generators: int, base_seed: int or None, gen_dir: str or None):

        try:

            self._generators = self._load(gen_dir)
            self.n_generators = len(self._generators)

            assert n_generators == len(self._generators), [n_generators, len(self._generators)]

        except (FileNotFoundError, TypeError) as e:

            seed_sequence = SeedSequence(base_seed)

            if base_seed is None:
                print("-- Initializing sampler rng(s) using system entropy...")
            else:
                print("-- Initializing sampler rng(s) with sequence seed: {}".format(base_seed))
            print("   -> Entropy: {}".format(seed_sequence.entropy))

            seed_sequence = seed_sequence

            child_seqs = seed_sequence.spawn(n_generators)

            self._generators = [PCG64(seq) for seq in child_seqs]

            self.n_generators = n_generators

            if not isinstance(e, TypeError):
                self._cache(gen_dir)

    def _load(self, gen_dir):

        gen_path = os.path.join(gen_dir, "gen.states")

        with open(gen_path, "rb") as f:
            return pk.load(f)

    def _cache(self, gen_dir):
        gen_path = os.path.join(gen_dir, "gen.states")

        with open(gen_path, "wb") as f:
            pk.dump(self._generators, f)

    def rng(self, *args, **kwargs):
        assert 0, "override method"


class RNGAPriori(_Generators):

    def __init__(self, n_sample_params: int, base_seed: int or None = None, gen_dir: str or None = None):
        super(RNGAPriori, self).__init__(n_sample_params, base_seed, gen_dir)

    def rng(self, param_index, n):
        assert param_index < self.n_generators, [param_index, self.n_generators]

        return default_rng(self._generators[param_index].jumped(n))


class RNGLatin(_Generators):

    def __init__(self, n_sample_params: int, base_seed: int or None = None, gen_dir: str or None = None):
        super(RNGLatin, self).__init__(n_sample_params + 1, base_seed, gen_dir)

    def rng(self, param_index, n):
        assert param_index < self.n_generators - 1, [param_index, self.n_generators - 1]

        return default_rng(self._generators[param_index].jumped(n))

    def grid_rng(self):
        return default_rng(self._generators[-1])


if __name__ == "__main__":
    import os
    import scipy.stats as stats

    seed = 121212
    test = RNGAPriori(5, seed)

    g1 = stats.uniform(0, 1)
    g2 = stats.uniform(0, 1)

    g1.random_state = test.rng(0, 0)
    g2.random_state = test.rng(1, 0)

    print(g1.rvs(1))
    print(g2.rvs(1))
    #
    # del test
    #
    # test = load_states("test.pk")
    #
    # g1.random_state = test.rng(0)
    # g2.random_state = test.rng(0)
    # g3.random_state = test.rng(0)
    #
    # print(g1.rvs(1))
    # print(g2.rvs(1))
    # print(g3.rvs(1))
    #
    # os.remove("test.pk")
    #
    # p = _ParamStates(5, 10)
    #
    # g1.random_state = p.rng(0, 11)
    # g2.random_state = p.rng(0, 11)
    #
    # print(g1.rvs(1))
    # print(g2.rvs(1))
    #
    # print(p.rng(0, 100).uniform(0, 1))
    # print(p.rng(0, 100).uniform(0, 1))

# class _RandomStates(object):
#
#     def __init__(self, n_states, seed=None):
#         """
#         Spawns sequence of random states from an initial seeding value
#         :param seed: initial seed value, if None assigned from system entropy
#         :param n_states: number of child seeds to spawn
#         """
#
#         self.n_states = n_states
#         rng = np.random.default_rng(seed)
#         seq = rng.bit_generator._seed_seq
#         self.states = seq.spawn(self.n_states)  # builds new uncorrelated states
#
#     def get_rng(self, idx: int):
#         """
#         Gets random number generator with an index i, belonging to n states
#         :param idx: index of random state
#         :return: Generator instance
#         """
#         assert idx < self.n_states, [idx, self.n_states]
#         return np.random.default_rng(self.states[idx])
#
# class APrioriStates(_RandomStates):
#
#     def __init__(self, n_samples, n_params, n_states):
#
#
#
#         super().__init__(n_states)
