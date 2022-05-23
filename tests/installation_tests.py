import unittest
import numpy as np
import pyt_quad
import pyt_dquad_euclidean
import pyt_dquad_2sphere

_GLOBAL_TOLS = np.array([1e-4, 1e-4])  # tols for testing; relatively weak for speed. Tests may fail if edited


class ModelDefs:

    # Typehints for (expected) compiled module behaviour

    @staticmethod
    def nF() -> int:
        pass

    @staticmethod
    def nP() -> int:
        pass

    @staticmethod
    def H(*args, **kwargs) -> float:
        pass

    @staticmethod
    def Epsilon(*args, **kwargs) -> float:
        pass

    @staticmethod
    def Eta(*args, **kwargs) -> float:
        pass

    @staticmethod
    def findEndOfInflation(*args, **kwargs) -> float:
        pass

    @staticmethod
    def V(*args, **kwargs) -> float:
        pass

    @staticmethod
    def dV(*args, **kwargs) -> np.ndarray:
        pass

    @staticmethod
    def ddV(*args, **kwargs) -> np.ndarray:
        pass

    @staticmethod
    def backEvolve(*args, **kwargs) -> np.ndarray:
        pass

    @staticmethod
    def sigEvolve(*args, **kwargs) -> np.ndarray:
        pass

    @staticmethod
    def alphaEvolve(*args, **kwargs) -> np.ndarray:
        pass

    @staticmethod
    def massMatrix(*args, **kwargs) -> np.ndarray:
        pass


class _Expectations:
    nF: int
    nP: int
    Hi: float
    Nend: float
    Epsi: float
    Etai: float
    Vi: float
    dVi: np.ndarray
    ddVi: np.ndarray
    back: np.ndarray
    sig: np.ndarray
    alpha: np.ndarray
    masses: np.ndarray


class TestQuadratic(unittest.TestCase):
    # Test quadratic model

    pyt_model = pyt_quad
    ics = np.array([20., 1e-4], dtype=float)
    pars = np.array([1e-3], dtype=float)
    n_end = 101.17024620322667
    n_arr = np.linspace(0, n_end, 1000)
    eps_arr = np.array([7.499812504629944e-05, 1.1057319168321953], dtype=float)

    def get_back(self):
        return self.pyt_model.backEvolve(self.n_arr, self.ics, self.pars, _GLOBAL_TOLS, True)

    def test_nf(self):
        assert self.pyt_model.nF() == 1

    def test_np(self):
        assert self.pyt_model.nP() == 1

    def test_n_end(self):
        assert np.isclose(self.pyt_model.findEndOfInflation(self.ics, self.pars, _GLOBAL_TOLS), self.n_end)

    def test_eps(self):
        back = self.get_back()

        eps_start = self.pyt_model.Epsilon(back[0][1:1 + 2 * self.pyt_model.nF()], self.pars)
        eps_end = self.pyt_model.Epsilon(back[-1][1:1 + 2 * self.pyt_model.nF()], self.pars)

        assert np.allclose(self.eps_arr, np.array([eps_start, eps_end]))


class TestDoubleQuadratic(unittest.TestCase):
    # Test double quadratic model

    pyt_model = pyt_dquad_euclidean

    def test_nf(self):
        assert self.pyt_model.nF() == 2

    def test_np(self):
        assert self.pyt_model.nP() == 2


class TestDoubleQuadratic2Sphere(unittest.TestCase):
    # Test double quadratic model

    pyt_model = pyt_dquad_2sphere

    def test_nf(self):
        assert self.pyt_model.nF() == 2

    def test_np(self):
        assert self.pyt_model.nP() == 3


if __name__ == "__main__":
    unittest.main()
