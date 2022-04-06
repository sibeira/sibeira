import numpy
import scipy.constants
import unittest

from sibeira.rate_profile import RateProfile
from sibeira.integrator import RateIntegrator
from sibeira.beam import Beam


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


class TestRateProfileSpline(unittest.TestCase):

    def test_spline_for_0(self):
        r = RateProfile('Li', 40)
        s = r.get_spline([10, 20, 50, 100], [1, 3, 5, 10])
        numpy.testing.assert_array_almost_equal([0.0], s([0.]), decimal=15, err_msg='spline test for 0')

    def test_spline_for_1p1(self):
        r = RateProfile('Li', 40)
        s = r.get_spline([10, 20, 50, 100], [1, 3, 5, 10])
        numpy.testing.assert_array_almost_equal([3.135669805050988e-08], s([1.1]), decimal=15, err_msg='spline test for 1.1')

    def test_spline_for_11(self):
        r = RateProfile('Li', 40)
        s = r.get_spline([10, 20, 50, 100], [1, 3, 5, 10])
        numpy.testing.assert_array_almost_equal([1.2481484271391177], s([11.]), decimal=4, err_msg='spline test for 11')
