import unittest

import numpy

from renate_od import get_export_name, get_scenario_latex, plot_attenuation_profile


class TestRenateOD(unittest.TestCase):
    def test_export_name(self):
        print(get_export_name('plot', 12345, 678, 'test', 123))
        print(get_export_name('log', 12345, 678, 'test', 123))
        print(get_export_name('lorem_ipsum', 12345, 678, 'test', 123))
        print(get_export_name('dolor_sit_amet', 12345, 678, 'test', 123, 3))
        print(get_export_name('dolor_sit_amet', 12345, 678, 'test', 123, 3, 'electron'))

    def test_scenario_latex(self):
        print(get_scenario_latex('total'))
        print(get_scenario_latex('just electron'))
        print(get_scenario_latex('just ion'))
        print(get_scenario_latex('lorem ipsum'))

    def test_scenario_path(self):
        print(get_scenario_latex('total'))
        print(get_scenario_latex('just electron'))
        print(get_scenario_latex('just ion'))
        print(get_scenario_latex('lorem ipsum'))


class TestRenateODPlot(unittest.TestCase):
    def test_plot_attenuation_profile(self):
        r = numpy.linspace(0.6, 0.74, 11)
        p1 = numpy.random.random_sample(r.shape)
        p2 = numpy.random.random_sample(r.shape)
        p3 = numpy.random.random_sample(r.shape)
        plot_attenuation_profile(12345, 678, 'test', 123, 3, r, [p1, p2, p3], ['p1', 'p2', 'p3'])


if __name__ == '__main__':
    unittest.main()
