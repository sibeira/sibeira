from tests.test_integrator import *


class TestIntegratorNormalisation3D(TestIntegratorNormalisation):
    def test_normalisation_3d(self):
        rate = RateIntegrator('charge exchange', 'Li', 40, 100, 3)
        normalisation_factor = rate.integrate(rate.integrand_normalisation)
        numpy.testing.assert_array_almost_equal(normalisation_factor, 2.0*scipy.constants.pi**2, decimal=4,
                                                err_msg='3D normalisation factor')


class TestIntegratorSpeed3D(TestIntegratorSpeed):
    def test_mean_speed_3d(self):
        temperature = 100
        rate = RateIntegrator('electron impact ionisation', 'Li', 40, temperature, 3)
        rate.cross_section = self.one
        mean_speed = rate.get_coefficient()
        expected = numpy.sqrt(8.0 * scipy.constants.e * temperature / scipy.constants.pi / scipy.constants.electron_mass)
        numpy.testing.assert_approx_equal(mean_speed, expected, significant=5,
                                          err_msg='1D normalisation factor')


if __name__ == '__main__':
    unittest.main()
