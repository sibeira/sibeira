from species import *
import unittest


class TestSpecies(unittest.TestCase):
    def test_hydrogen(self):
        o = OrbitalConstants('H')
        self.assertEqual(RYDBERG, o.get('B'), 'B for hydrogen')
        self.assertEqual(RYDBERG, o.get('U'), 'U for hydrogen')
        self.assertEqual(1, o.get('N'),  'N for hydrogen')
        self.assertEqual(0.5668, o.get('Q'),  'Q for hydrogen')
        self.assertEqual(1, o.get('n'), 'n for hydrogen')

    def test_helium(self):
        o = OrbitalConstants('He')
        self.assertEqual(24.587389011, o.get('B'), 'B for helium')
        self.assertEqual(3*RYDBERG, o.get('U'), 'U for helium')
        self.assertEqual(2, o.get('N'), 'N for helium')
        self.assertEqual(0.8841, o.get('Q'), 'Q for helium')
        self.assertEqual(1, o.get('n'), 'n for helium')

    def test_caesium(self):
        o = OrbitalConstants('Cs')
        self.assertEqual(3.89390572743, o.get('B'), 'B for caesium')
        self.assertEqual(3.89390572743, o.get('U'), 'U for caesium')
        self.assertEqual(1, o.get('N'), 'N for caesium')
        self.assertEqual(1, o.get('Q'), 'Q for caesium')
        self.assertEqual(1, o.get('n'), 'n for caesium')

    def test_lithium_ion(self):
        o = IonOrbitalConstants('Li')
        self.assertEqual(75.6400970, o.get('B'), 'B for Li+')
        self.assertEqual(14.56*RYDBERG, o.get('U'), 'U for Li+')
        self.assertEqual(2, o.get('N'), 'N for Li+')
        self.assertEqual(1, o.get('Q'), 'Q for Li+')
        self.assertEqual(1, o.get('n'), 'n for Li+')
        self.assertEqual(3, o.get('Z'), 'Z for Li+')
        o.fill_U()
        self.assertAlmostEqual(2.0*7.223*RYDBERG, o.get('U'), places=1, msg='U from variation theory for Li+')


if __name__ == '__main__':
    unittest.main()
