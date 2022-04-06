import numpy
import scipy.constants
import unittest

from sibeira.beam import Beam


class TestRate(unittest.TestCase):

    @staticmethod
    def get_impact_energy(mass, velocity):
        return 0.5 * mass * velocity ** 2 / scipy.constants.elementary_charge

    @staticmethod
    def get_velocity(mass, impact_energy):
        return numpy.sqrt(2.0 * impact_energy * scipy.constants.elementary_charge / mass)

    @staticmethod
    def get_deuterium_mass():
        return Beam('D', 0).get_mass()


if __name__ == '__main__':
    unittest.main()

