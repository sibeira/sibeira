from species import *
import unittest


class TestSpecies(unittest.TestCase):
    def test_hydrogen(self):
        o = OrbitalConstants('H')
        self.assertEqual(o.get('B'), RYDBERG, 'U for hydrogen')
        self.assertEqual(o.get('U'), RYDBERG, 'U for hydrogen')
        self.assertEqual(o.get('N'), 1, 'N for hydrogen')
        self.assertEqual(o.get('Q'), 0.5668, 'Q for hydrogen')
        self.assertEqual(o.get('n'), 1, 'n for hydrogen')

    def test_helium(self):
        o = OrbitalConstants('He')
        self.assertEqual(o.get('B'), 24.587, 'U for hydrogen')
        self.assertEqual(o.get('U'), 3*RYDBERG, 'U for hydrogen')
        self.assertEqual(o.get('N'), 2, 'N for hydrogen')
        self.assertEqual(o.get('Q'), 0.8841, 'Q for hydrogen')
        self.assertEqual(o.get('n'), 1, 'n for hydrogen')


if __name__ == '__main__':
    unittest.main()
