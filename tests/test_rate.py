import scipy.constants
import unittest
from bebim.rate import *
from bebim.integrator import RateIntegrator


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
        rate = RateIntegrator('electron impact ionisation', 'Li', 40, 100, 1)
        normalisation_factor = rate.integrate(rate.integrand_normalisation)
        numpy.testing.assert_array_almost_equal(normalisation_factor, 1, decimal=4,
                                                err_msg='1D normalisation factor')

    def test_normalisation_2d(self):
        rate = RateIntegrator('charge exchange', 'Li', 40, 100, 2)
        normalisation_factor = rate.integrate(rate.integrand_normalisation)
        numpy.testing.assert_array_almost_equal(normalisation_factor, 2.0*scipy.constants.pi, decimal=4,
                                                err_msg='2D normalisation factor')
