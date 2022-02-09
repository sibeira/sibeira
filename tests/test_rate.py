import unittest
from bebim.rate import *
from bebim.integrator import Rate1DIntegralElectron, Rate2DIntegralIon


class TestRate(unittest.TestCase):

    @staticmethod
    def get_impact_energy(mass, velocity):
        return 0.5 * mass * velocity ** 2 / scipy.constants.elementary_charge

    @staticmethod
    def get_velocity(mass, impact_energy):
        return numpy.sqrt(2.0 * impact_energy * scipy.constants.elementary_charge / mass)

    @staticmethod
    def get_deuterium_mass():
        return Beam('D', 0).get_mass()

    def test_normalisation_1d(self):
        rate = Rate1DIntegralElectron(10, 10, lambda x: 1.0/self.get_velocity(scipy.constants.electron_mass, x))
        numpy.testing.assert_array_almost_equal(rate.get_coefficient(), 1, decimal=4,
                                                err_msg='1D normalisation factor')

    def test_normalisation_2d(self):
        rate = Rate2DIntegralIon(10, 10, lambda x: 1.0/self.get_velocity(self.get_deuterium_mass(), x))
        numpy.testing.assert_array_almost_equal(rate.get_coefficient(), 1, decimal=4,
                                                err_msg='2D normalisation factor')
