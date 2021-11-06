import numpy
import unittest
from beam import Beam


class TestBeam(unittest.TestCase):
    def test_string_energy(self):
        self.assertRaises(TypeError, Beam, 'Li', 'test', 0, 0)

    def test_negative_energy(self):
        self.assertRaises(ValueError, Beam, 'Li', -1.3, 0, 0)

    def test_fake_species(self):
        self.assertRaises(ValueError, Beam, 'fake_species', 10, 0, 0)

    def test_mass_lithium(self):
        b = Beam('Li', 10, 0, 0)
        numpy.testing.assert_array_almost_equal(b.get_mass(), 1.165e-26, decimal=30)

    def test_speed_lithium(self):
        b = Beam('Li', 10000, 0, 0)
        numpy.testing.assert_array_almost_equal(b.get_speed(), 524446, decimal=0)


if __name__ == '__main__':
    unittest.main()
