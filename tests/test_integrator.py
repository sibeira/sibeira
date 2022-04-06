import unittest

import numpy
import scipy.constants

from sibeira.integrator import RateIntegrator
from sibeira.beam import Beam


class TestIntegrator(unittest.TestCase):
    def test_target_mass_electron_impact_ionisation(self):
        r = RateIntegrator('electron impact ionisation', 'Li', 0, 0)
        numpy.testing.assert_array_almost_equal(r.get_target_mass(), scipy.constants.electron_mass,
                                                err_msg='Wrong target mass for electron impact ionisation')

    def test_target_mass_charge_exchange(self):
        plasma_ion_mass = Beam('D', 0).get_mass()
        r = RateIntegrator('charge exchange', 'Li', 0, 0)
        numpy.testing.assert_array_almost_equal(r.get_target_mass(), plasma_ion_mass,
                                                err_msg='Wrong target mass for charge exchange reaction')

    def test_target_reaction_invalid(self):
        self.assertRaises(ValueError, RateIntegrator, 'unknown reaction', 'Li', 0, 0)

    def test_target_mass_invalid(self):
        self.assertRaises(ValueError, RateIntegrator, 'charge exchange', 'unknown species', 0, 0)

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