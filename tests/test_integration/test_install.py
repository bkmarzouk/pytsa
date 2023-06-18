import unittest
import numpy as np

try:
    from pytsa.models import dquad_2sphere as test_model
except ImportError:
    assert 0, "Unable to do tests, no module"


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


class TestDoubleQuadratic2Sphere(unittest.TestCase):
    # Test compiled model file

    tols = np.array([1e-5, 1e-5])
    ics = np.array(
        [
            -1.5890304023505983,
            -18.389416025688845,
            0.00232248417286116,
            0.00418422599903375,
        ]
    )
    pars = np.array(
        [0.01329597888476721, 0.005246064625353842, 5.195022895180804]
    )

    def _get_back(self):
        N_evo = np.linspace(0, 3000, 20 * 3000)

        back = test_model.backEvolve(
            N_evo, self.ics, self.pars, self.tols, True, -1
        )

        return back

    def test_import(self):
        assert test_model is not None, test_model

    def test_nf(self):
        assert test_model.nF() == 2

    def test_np(self):
        assert test_model.nP() == 3

    def test_back(self):
        back = self._get_back()

        expt = np.array(
            [
                8.05363423e02,
                -1.28122309e-01,
                1.34744008e-01,
                2.48930507e-04,
                -1.14051532e-03,
            ]
        )

        assert np.allclose(back[-1], expt)

    def test_eps(self):
        back = self._get_back()

        eps = test_model.Epsilon(back[-1][1:], self.pars)

        expt = 1.1929047258540073

        assert np.isclose(eps, expt), [eps, expt]


if __name__ == "__main__":
    unittest.main()
