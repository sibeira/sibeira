import numpy
import scipy.constants
import scipy.integrate
import scipy.stats
from cross_section import CrossSection
from tabata_ctf.cross_section import CrossSection as CXCrossSection


class Beam:
    def __init__(self, species, beam_energy,
                 electron_temperature=numpy.nan, electron_density=numpy.nan, ionisation_level=0):
        self.beam_energy = beam_energy
        self.species = species
        self.mass = self.get_mass()
        self.speed = self.get_speed()
        self.electron_temperature = electron_temperature
        self.electron_density = electron_density
        self.ionisation_level = ionisation_level

        self.beb = CrossSection(1000, self.species, self.ionisation_level)
        self.beb.set_polynomial()
        self.tabata = CXCrossSection(1000, self.species)
        self.tabata.set_polynomial()

    def set_profiles(self, electron_temperature=numpy.nan, electron_density=numpy.nan):
        if ~numpy.isnan(electron_temperature):
            self.electron_temperature = electron_temperature
        if ~numpy.isnan(electron_density):
            self.electron_density = electron_density

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

    ######################

    def integrand_1d_coefficient(self, v):
        velocity = v * self.velocity_normalisation_factor
        kinetic_energy = 0.5 * scipy.constants.electron_mass * velocity ** 2 / scipy.constants.elementary_charge
        return scipy.stats.maxwell.pdf(v) * velocity * (self.beb.f(kinetic_energy) + self.tabata.f(kinetic_energy))

    def integrand_1d_coefficient_tabata(self, v):
        velocity = v * self.velocity_normalisation_factor
        kinetic_energy = 0.5 * scipy.constants.electron_mass * velocity ** 2 / scipy.constants.elementary_charge
        return scipy.stats.maxwell.pdf(v) * velocity * self.tabata.f(kinetic_energy)

    def get_1d_coefficient(self):
        self.velocity_normalisation_factor = numpy.sqrt(
            self.electron_temperature * scipy.constants.elementary_charge / scipy.constants.electron_mass)

        c = self.get_1d_normalisation()
        return c * scipy.integrate.quad(self.integrand_1d_coefficient, 0, numpy.inf)[0]

    def get_1d_coefficient_tabata(self):
        self.velocity_normalisation_factor = numpy.sqrt(
            self.electron_temperature * scipy.constants.elementary_charge / scipy.constants.electron_mass)

        c = self.get_1d_normalisation()
        return c * scipy.integrate.quad(self.integrand_1d_coefficient_tabata, 0, numpy.inf)[0]

    @staticmethod
    def integrand_1d_normalisation(v):
        return scipy.stats.maxwell.pdf(v)

    def get_1d_normalisation(self):
        try:
            return self.value_1d_normalisation
        except AttributeError:
            self.value_1d_normalisation = \
                scipy.integrate.quad(self.integrand_1d_normalisation, 0, numpy.inf)[0]
            return self.value_1d_normalisation

    def get_attenuation(self):
        if self.electron_density == 0:
            return 0.0
        rate_coefficient = self.get_1d_coefficient()
        return rate_coefficient * self.electron_density

    ######################

    @staticmethod
    def get_third_side_length(a, b, alpha, beta):
        return numpy.sqrt(a ** 2 + b ** 2 + 2 * a * b * numpy.cos(alpha) * numpy.cos(beta))

    def integrand_3d_coefficient(self, beta, alpha, v):
        velocity = self.get_third_side_length(v * self.velocity_normalisation_factor, self.speed, alpha, beta)
        kinetic_energy = 0.5 * scipy.constants.electron_mass * velocity ** 2 / scipy.constants.elementary_charge
        return scipy.stats.maxwell.pdf(v) * velocity * (self.beb.f(kinetic_energy) + self.tabata.f(kinetic_energy))

    @staticmethod
    def integrand_3d_normalisation(beta, alpha, v):
        return scipy.stats.maxwell.pdf(v)

    def get_3d_normalisation(self):
        try:
            return self.value_3d_normalisation
        except AttributeError:
            self.value_3d_normalisation = \
                scipy.integrate.tplquad(self.integrand_3d_normalisation, 0, numpy.inf, -numpy.pi, numpy.pi,
                                        -numpy.pi / 2, numpy.pi / 2)[0]
            return self.value_3d_normalisation

    def get_3d_coefficient(self):
        self.velocity_normalisation_factor = numpy.sqrt(
            self.electron_temperature * scipy.constants.elementary_charge / scipy.constants.electron_mass)
        return scipy.integrate.tplquad \
               (self.integrand_3d_coefficient, 0, numpy.inf, -numpy.pi, numpy.pi, -numpy.pi / 2, numpy.pi / 2)[0] / \
               self.get_3d_normalisation()

    def get_attenuation_3d(self):
        if self.electron_density == 0:
            return 0.0
        rate_coefficient = self.get_3d_coefficient()
        return rate_coefficient * self.electron_density

    def get_attenuation_nrl(self, is_with_tabata=False):
        c = CrossSection(self.electron_temperature, self.species, self.ionisation_level)
        t = c.get_t()
        r = 1e-11 * numpy.sqrt(t) / c.B ** 1.5 / (6.0 + t) * numpy.exp(-1.0 / t)
        if is_with_tabata:
            r += self.get_1d_coefficient_tabata()
        return r * self.electron_density
