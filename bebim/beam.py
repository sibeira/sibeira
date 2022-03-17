import numpy
import scipy.constants
import scipy.integrate
import scipy.stats


class Beam:
    def __init__(self, species, beam_energy, ionisation_level=0):
        self.beam_energy = beam_energy
        self.species = species
        self.mass = self.get_mass()
        self.speed = self.get_speed()
        self.ionisation_level = ionisation_level

    def get_mass(self):
        if self.species == 'D':
            return 2.0141017781212 * scipy.constants.physical_constants['atomic mass constant'][0]
        elif self.species == 'Li':
            return 7.016003436645 * scipy.constants.physical_constants['atomic mass constant'][0]
        elif self.species == 'Na':
            return 22.989769282019 * scipy.constants.physical_constants['atomic mass constant'][0]
        elif self.species == 'K':
            return 38.963706486449 * scipy.constants.physical_constants['atomic mass constant'][0]
        elif self.species == 'Rb':
            return 84.911789737954 * scipy.constants.physical_constants['atomic mass constant'][0]
        elif self.species == 'Cs':
            return 132.905451961080 * scipy.constants.physical_constants['atomic mass constant'][0]
        else:
            raise ValueError('Illegal input species: ' + self.species)

    def get_speed(self):
        if not (type(self.beam_energy) == int or float):
            raise TypeError('Energy is not numeric: ' + str(self.beam_energy))
        if self.beam_energy < 0.0:
            raise ValueError('Energy cannot be a negative value! (' + str(self.beam_energy) + ')')
        return numpy.sqrt(2.0 * self.beam_energy * scipy.constants.elementary_charge / self.mass)
