import numpy
import unittest
import matplotlib.pyplot


class TestCrossSection(unittest.TestCase):

    def test_hydrogen(self):
        reference_energy, reference_cross_section = self.read_reference('H')

        cross_section = reference_cross_section
        energy = reference_energy

        numpy.testing.assert_equal(cross_section, reference_cross_section, 'Cross section test for hydrogen')
        self.plot(energy, cross_section, reference_energy, reference_cross_section)

    @staticmethod
    def read_reference(species):
        reference_data = numpy.loadtxt('reference_cross_section_' + species + '.dat')
        energy = reference_data[:, 0]
        cross_section = reference_data[:, 1]
        return energy, cross_section

    @staticmethod
    def plot(energy, cross_section, reference_energy, reference_cross_section):
        matplotlib.pyplot.semilogx(energy, cross_section, '.')
        matplotlib.pyplot.semilogx(reference_energy, reference_cross_section)
        matplotlib.pyplot.show()


if __name__ == '__main__':
    unittest.main()
