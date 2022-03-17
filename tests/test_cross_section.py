import numpy
import unittest
from bebim.cross_section import CrossSection
from tests.read_reference import read_reference


class TestCrossSection(unittest.TestCase):
    def test_hydrogen_with_Q(self):
        reference_energy, reference_cross_section = read_reference('H_with_Q')

        energy = reference_energy
        c = CrossSection('H', 0)
        cross_section = c.calculate(energy) / 1e-20

        numpy.testing.assert_array_almost_equal(cross_section, reference_cross_section, decimal=4,
                                                err_msg='Cross section test for hydrogen')

    def test_hydrogen(self):
        reference_energy, reference_cross_section = read_reference('H')

        energy = reference_energy
        c = CrossSection('H', 0, with_Q=False)
        cross_section = c.calculate(energy) / 1e-20

        numpy.testing.assert_array_almost_equal(cross_section, reference_cross_section, decimal=4,
                                                err_msg='Cross section test for hydrogen (Q=1)')

    def test_helium_with_Q(self):
        reference_energy, reference_cross_section = read_reference('He_with_Q')

        energy = reference_energy
        c = CrossSection('He', 0)
        cross_section = c.calculate(energy) / 1e-20

        numpy.testing.assert_array_almost_equal(cross_section, reference_cross_section, decimal=3,
                                                err_msg='Cross section test for helium')

    def test_helium(self):
        reference_energy, reference_cross_section = read_reference('He')

        energy = reference_energy
        c = CrossSection('He', 0, with_Q=False)
        cross_section = c.calculate(energy) / 1e-20

        numpy.testing.assert_array_almost_equal(cross_section, reference_cross_section, decimal=3,
                                                err_msg='Cross section test for helium (Q=1)')


if __name__ == '__main__':
    unittest.main()
