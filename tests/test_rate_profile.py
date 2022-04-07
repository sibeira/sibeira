import unittest

import numpy
import tempfile

from sibeira.rate_profile import RateProfile


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


class TestRateProfileIO(unittest.TestCase):
    def test_import_profile_not_exist(self):
        r = RateProfile('Li', 40)
        with self.assertRaises(ValueError) as i:
            r.import_profile('invalid profile', 0)
        self.assertEqual('The profile is not found: invalid profile (Tabata 0D)', str(i.exception))

    def test_import_profile_not_exist_tabata_off(self):
        r = RateProfile('Li', 40)
        with self.assertRaises(ValueError) as i:
            r.import_profile('invalid profile', -1)
        self.assertEqual('The profile is not found: invalid profile (Tabata OFF)', str(i.exception))


if __name__ == '__main__':
    unittest.main()
