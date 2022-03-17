import unittest
import numpy

from tabata_ctf.cross_section import CrossSection


class TestTabataCrossSection(unittest.TestCase):
    def test_fake_species(self):
        self.assertRaises(ValueError, CrossSection, 'fake_species')

    def test_broken_database(self):
        c = CrossSection('Li')
        c.tabata_data.__delitem__('Et')
        self.assertRaises(KeyError, c.calculate, 1000)

    def test_low_energy(self):
        energy = numpy.logspace(0, 6, 50)
        c = CrossSection('Na', 'double')
        value = c.calculate(energy)[0]
        self.assertFalse(numpy.isnan(value))


if __name__ == '__main__':
    unittest.main()
