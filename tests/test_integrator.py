import unittest

import numpy
import scipy.constants

from integrator import RateIntegrator
from sibeira.beam import Beam


class TestIntegrator(unittest.TestCase):
    def test_target_mass_electron_impact_ionisation(self):
        r = RateIntegrator('electron impact ionisation', 0, 0)
        numpy.testing.assert_array_almost_equal(r.get_target_mass(), scipy.constants.electron_mass,
                                                err_msg='Wrong target mass for electron impact ionisation')

    def test_target_mass_charge_exchange(self):
        plasma_ion_mass = Beam('D', 0).get_mass()
        r = RateIntegrator('charge exchange', 0, 0)
        numpy.testing.assert_array_almost_equal(r.get_target_mass(), plasma_ion_mass,
                                                err_msg='Wrong target mass for charge exchange reaction')

    def test_target_mass_invalid(self):
        self.assertRaises(ValueError, RateIntegrator, 'unknown species', 0, 0)

