import unittest

import numpy

from sibeira.rate_profile import RateProfile


class TestRateProfileSpline(unittest.TestCase):

    def test_spline_for_0(self):
        r = RateProfile('Li', 40)
        s = r.get_spline([10, 20, 50, 100], [1, 3, 5, 10])
        result = r.resolve_log_spline(s, 0.)
        numpy.testing.assert_array_almost_equal([0.0], result, decimal=15, err_msg='spline test for 0')

    def test_spline_for_1p1(self):
        r = RateProfile('Li', 40)
        s = r.get_spline([10, 20, 50, 100], [1, 3, 5, 10])
        result = r.resolve_log_spline(s, 1.1)
        numpy.testing.assert_array_almost_equal([3.135669805050988e-08], result, decimal=15, err_msg='spline test for 1.1')

    def test_spline_for_11(self):
        r = RateProfile('Li', 40)
        s = r.get_spline([10, 20, 50, 100], [1, 3, 5, 10])
        result = r.resolve_log_spline(s, 11.)
        numpy.testing.assert_array_almost_equal([1.2481484271391177], result, decimal=4, err_msg='spline test for 11')


if __name__ == '__main__':
    unittest.main()
