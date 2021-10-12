import numpy
import unittest
import matplotlib.pyplot
from cross_section import CrossSection


def plot_profile(an_energy, a_cross_section, reference_energy, reference_cross_section):
    matplotlib.pyplot.semilogx(an_energy, a_cross_section, '.')
    matplotlib.pyplot.semilogx(reference_energy, reference_cross_section)
    matplotlib.pyplot.show()


def read_reference(species):
    reference_data = numpy.loadtxt('reference_cross_section_' + species + '.dat')
    energy = reference_data[:, 0]
    a_cross_section = reference_data[:, 1]
    return energy, a_cross_section


class TestCrossSection(unittest.TestCase):

    def test_hydrogen_with_Q(self):
        reference_energy, reference_cross_section = read_reference('H_with_Q')

        an_energy = reference_energy
        c = CrossSection(an_energy, 'H', 0)
        a_cross_section = c.calculate() / 1e-20

        numpy.testing.assert_array_almost_equal(a_cross_section, reference_cross_section, decimal=4,
                                                err_msg='Cross section test for hydrogen')

    def test_hydrogen(self):
        reference_energy, reference_cross_section = read_reference('H')

        an_energy = reference_energy
        c = CrossSection(an_energy, 'H', 0, with_Q=False)
        a_cross_section = c.calculate() / 1e-20

        numpy.testing.assert_array_almost_equal(a_cross_section, reference_cross_section, decimal=4,
                                                err_msg='Cross section test for hydrogen (Q=1)')


if __name__ == '__main__':
    unittest.main()
