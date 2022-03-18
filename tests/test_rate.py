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

    def test_normalisation_3d(self):
        rate = RateIntegrator('charge exchange', 'Li', 40, 100, 3)
        normalisation_factor = rate.integrate(rate.integrand_normalisation)
        numpy.testing.assert_array_almost_equal(normalisation_factor, 2.0*scipy.constants.pi**2, decimal=4,
                                                err_msg='3D normalisation factor')


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

