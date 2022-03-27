import unittest

from renate_od import get_export_name, get_scenario_latex, get_scenario_path


class TestRenateODExportName(unittest.TestCase):
    def test_export_name(self):
        print(get_export_name('plot', 12345, 678, 'test', 123))
        print(get_export_name('log', 12345, 678, 'test', 123))
        print(get_export_name('lorem_ipsum', 12345, 678, 'test', 123))
        print(get_export_name('dolor_sit_amet', 12345, 678, 'test', 123, 3))
        print(get_export_name('dolor_sit_amet', 12345, 678, 'test', 123, 3, 'just electron'))


class TestRenateODScenarioLatex(unittest.TestCase):
    def test_total(self):
        reference = ''
        self.assertEqual(get_scenario_path('total'), reference, 'Scenario LaTeX test for combnined models')

    def test_just_electron(self):
        reference = ''
        self.assertEqual(get_scenario_path('just electron'), reference, 'Scenario LaTeX test for pure electron models')

    def test_just_ion(self):
        reference = ''
        self.assertEqual(get_scenario_path('just electron'), reference, 'Scenario LaTeX test for pure ion models')

    def test_invalid(self):
        self.assertRaises(ValueError, get_scenario_path, 'Scenario LaTeX test invalid case')


class TestRenateODScenarioPath(unittest.TestCase):
    def test_total(self):
        reference = ''
        self.assertEqual(get_scenario_path('total'), reference, 'Scenario path test for combnined models')

    def test_just_electron(self):
        reference = ''
        self.assertEqual(get_scenario_path('just electron'), reference, 'Scenario path test for pure electron models')

    def test_just_ion(self):
        reference = ''
        self.assertEqual(get_scenario_path('just electron'), reference, 'Scenario path test for pure ion models')

    def test_invalid(self):
        self.assertRaises(ValueError, get_scenario_path, 'Scenario path test invalid case')


if __name__ == '__main__':
    unittest.main()
