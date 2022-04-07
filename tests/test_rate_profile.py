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
        with self.assertRaises(FileNotFoundError) as i:
            r.import_profile('invalid profile', 0)
        self.assertEqual('There is no profile for Li', str(i.exception))

    def test_export_and_import(self):
        reference_profile = 'test_data'
        profile_name = 'test'
        tabata_integration_dimension = -1

        directory = tempfile.mkdtemp()
        r = RateProfile('Li', 40)
        r.export_profile(profile_name, tabata_integration_dimension, reference_profile, directory)
        result_profile = r.import_profile(profile_name, tabata_integration_dimension, directory)

        self.assertEqual(reference_profile, result_profile, msg='Profile I/O does not work.')

    def test_export_and_import_fail_tabata(self):
        reference_profile = 'test_data'
        profile_name = 'test'

        directory = tempfile.mkdtemp()
        r = RateProfile('Li', 40)
        r.export_profile(profile_name, 1, reference_profile, directory)

        with self.assertRaises(KeyError) as i:
            r.import_profile(profile_name, -1, directory)

        self.assertEqual('The profile is not found: test (Tabata OFF)', i.exception.args[0],
                         msg='Profile I/O database error for invalid Tabata')

    def test_export_and_import_fail_profile(self):
        reference_profile = 'test_data'

        directory = tempfile.mkdtemp()
        r = RateProfile('Li', 40)
        r.export_profile('another_test', 1, reference_profile, directory)

        with self.assertRaises(KeyError) as i:
            r.import_profile('test', 1, directory)

        self.assertEqual('The profile is not found: test (Tabata 1D)', i.exception.args[0],
                         msg='Profile I/O database error for invalid Tabata')

    def test_export_and_import_database(self):
        reference_profile1 = 'test_data1'
        reference_profile2 = 'test_data2'
        reference_profile3 = 'test_data3'

        directory = tempfile.mkdtemp()
        r = RateProfile('Li', 40)
        r.export_profile('test', 1, reference_profile1, directory)
        r.export_profile('test', -1, reference_profile2, directory)
        r.export_profile('another', 0, reference_profile3, directory)
        result_profile1 = r.import_profile('test', 1, directory)
        result_profile2 = r.import_profile('test', -1, directory)
        result_profile3 = r.import_profile('another', 0, directory)

        self.assertEqual(reference_profile1, result_profile1, msg='Profile I/O database does not work. Case1')
        self.assertEqual(reference_profile2, result_profile2, msg='Profile I/O database does not work. Case2')
        self.assertEqual(reference_profile3, result_profile3, msg='Profile I/O database does not work. Case3')

    def test_get_beam_energy_as_string_int(self):
        reference = '40.0'
        r = RateProfile('Li', 40000)
        result = r.get_beam_energy_as_string()
        self.assertEqual(reference, result)

    def test_get_beam_energy_as_string_double(self):
        reference = '66.666'
        r = RateProfile('Li', 66666)
        result = r.get_beam_energy_as_string()
        self.assertEqual(reference, result)

    def test_get_file_name(self):
        reference = 'data/Li.npy'
        r = RateProfile('Li', 40)
        result = r.get_file_name('data')
        self.assertEqual(reference, result, msg='Test file name')


if __name__ == '__main__':
    unittest.main()
