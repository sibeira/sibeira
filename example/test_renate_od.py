import unittest
import os.path

from example.renate_od import get_export_name, get_title_name, get_scenario_latex, get_scenario_path


class TestRenateODExportName(unittest.TestCase):
    directory_path = os.path.dirname(os.path.abspath(__file__)) + '/figures/'

    def test_linear_plot(self):
        reference = self.directory_path + 'test_123keV_12345_678'
        self.assertEqual(reference, get_export_name('plot', 12345, 678, 'test', 123), 'Test linear plot path')

    def test_logarithmic_plot(self):
        reference = self.directory_path + 'test_123keV_12345_678_log'
        self.assertEqual(reference, get_export_name('log', 12345, 678, 'test', 123), 'Test logarithmic plot path')

    def test_other_plot(self):
        reference = self.directory_path + 'test_123keV_12345_678_lorem_ipsum'
        self.assertEqual(reference, get_export_name('lorem_ipsum', 12345, 678, 'test', 123), 'Test other plot path')

    def test_2d_integration_plot(self):
        reference = self.directory_path + 'test_123keV_12345_678_dolor_sit_amet_2d'
        self.assertEqual(reference, get_export_name('dolor_sit_amet', 12345, 678, 'test', 123, 2),
                         'Test 2D integration plot path')

    def test_3d_integration_plot(self):
        reference = self.directory_path + 'test_123keV_12345_678_dolor_sit_amet'
        self.assertEqual(reference, get_export_name('dolor_sit_amet', 12345, 678, 'test', 123, 3),
                         'Test 3D integration plot path')

    def test_3d_integration_plot_with_pure_electron_case(self):
        reference = self.directory_path + 'test_123keV_12345_678_electron_log'
        self.assertEqual(reference, get_export_name('log', 12345, 678, 'test', 123, 3, 'just electron'),
                         'Test pure electron case path with 3D integration with logarithmic plot')


class TestRenateODTitleName(unittest.TestCase):
    def test_default(self):
        reference = 'COMPASS #12345 (678 ms, test, 123 keV)'
        self.assertEqual(reference, get_title_name(12345, 678, 'test', 123), 'Test linear plot path')

    def test_just_ion(self):
        reference = 'COMPASS #12345 (678 ms, test, 123 keV, $n_e = 0$)'
        self.assertEqual(reference, get_title_name(12345, 678, 'test', 123, 'just ion'), 'Test logarithmic plot path')

    def test_just_electron(self):
        reference = 'COMPASS #12345 (678 ms, test, 123 keV, $n_i = 0$)'
        self.assertEqual(reference, get_title_name(12345, 678, 'test', 123, 'just electron'), 'Test other plot path')

    def test_invalid(self):
        self.assertRaises(ValueError, get_title_name, 12345, 678, 'test', 123, 'invalid')


class TestRenateODScenarioLatex(unittest.TestCase):
    def test_total(self):
        reference = ''
        self.assertEqual(reference, get_scenario_latex('total'), 'Scenario LaTeX test for combnined models')

    def test_just_electron(self):
        reference = ', $n_i = 0$'
        self.assertEqual(reference, get_scenario_latex('just electron'), 'Scenario LaTeX test for pure electron models')

    def test_just_ion(self):
        reference = ', $n_e = 0$'
        self.assertEqual(reference, get_scenario_latex('just ion'), 'Scenario LaTeX test for pure ion models')

    def test_invalid(self):
        self.assertRaises(ValueError, get_scenario_latex, 'Scenario LaTeX test invalid case')


class TestRenateODScenarioPath(unittest.TestCase):
    def test_total(self):
        reference = ''
        self.assertEqual(reference, get_scenario_path('total'), 'Scenario path test for combnined models')

    def test_just_electron(self):
        reference = '_electron'
        self.assertEqual(reference, get_scenario_path('just electron'), 'Scenario path test for pure electron models')

    def test_just_ion(self):
        reference = '_ion'
        self.assertEqual(reference, get_scenario_path('just ion'), 'Scenario path test for pure ion models')

    def test_invalid(self):
        self.assertRaises(ValueError, get_scenario_path, 'Scenario path test invalid case')


if __name__ == '__main__':
    unittest.main()
