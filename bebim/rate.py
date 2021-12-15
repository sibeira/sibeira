import numpy
import scipy.constants
import scipy.integrate
import scipy.stats
from bebim.beam import Beam
from bebim.cross_section import CrossSection
from tabata_ctf.cross_section import CrossSection as CXCrossSection


class Rate(Beam):
    def __init__(self, species, beam_energy, ionisation_level=0):
        super().__init__(species, beam_energy, ionisation_level)
        self.beb = CrossSection(1000, self.species, self.ionisation_level)
        self.beb.set_polynomial()
        self.tabata = CXCrossSection(1000, self.species)
        self.tabata.set_polynomial()

    def set_profiles(self, electron_temperature=numpy.nan, electron_density=numpy.nan):
        if ~numpy.isnan(electron_temperature):
            self.electron_temperature = electron_temperature
        if ~numpy.isnan(electron_density):
            self.electron_density = electron_density

    def integrand_1d_coefficient(self, v):
        velocity = v * self.velocity_normalisation_factor
        kinetic_energy = 0.5 * scipy.constants.electron_mass * velocity ** 2 / scipy.constants.elementary_charge
        return scipy.stats.maxwell.pdf(v) * velocity * (self.beb.f(kinetic_energy) + self.tabata.f(kinetic_energy))

    def get_1d_coefficient(self):
        self.velocity_normalisation_factor = numpy.sqrt(
            self.electron_temperature * scipy.constants.elementary_charge / scipy.constants.electron_mass)

        c = self.get_1d_normalisation()
        return c * scipy.integrate.quad(self.integrand_1d_coefficient, 0, numpy.inf)[0]

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
        return r * self.electron_density
