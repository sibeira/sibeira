import numpy
import scipy.constants
import scipy.integrate
import scipy.stats

from beam import Beam


class RateNDIntegral:
    def __init__(self, electron_temperature, ion_velocity):
        self.normalisation_factor = 1

    @staticmethod
    def get_third_side_length(a, b, alpha):
        return numpy.sqrt(a ** 2 + b ** 2 - 2 * a * b * numpy.cos(alpha))

    @staticmethod
    def maxwell(v):
        return scipy.stats.maxwell.pdf(v)

    @staticmethod
    def get_impact_energy(mass, velocity):
        return 0.5 * mass * velocity ** 2 / scipy.constants.elementary_charge

    def get_coefficient(self):
        return self.integrate(self.integrand) / self.integrate(self.integrand_normalisation)


class RateNDIntegralElectron(RateNDIntegral):
    def __init__(self, electron_temperature, ion_velocity, cross_section):
        super().__init__(electron_temperature, ion_velocity)
        self.electron_velocity_normalisation_factor = self.get_electron_velocity_normalisation_factor(electron_temperature)
        self.cross_section = cross_section

    @staticmethod
    def get_electron_velocity_normalisation_factor(electron_temperature):
        return numpy.sqrt(electron_temperature * scipy.constants.elementary_charge / scipy.constants.electron_mass)

    def integrand(self, v):
        velocity = v * self.electron_velocity_normalisation_factor
        impact_energy = self.get_impact_energy(scipy.constants.electron_mass, velocity)
        return self.maxwell(v) * velocity * self.cross_section(impact_energy)


class RateNDIntegralIon(RateNDIntegral):
    def __init__(self, electron_temperature, ion_velocity, cross_section, cross_section_of_double_ionisation=lambda x: 0):
        self.deuterium_mass = self.get_deuterium_mass()
        self.ion_velocity = ion_velocity
        self.ion_velocity_normalisation_factor = self.get_ion_velocity_normalisation_factor(electron_temperature)
        self.cross_section = cross_section
        self.cross_section_of_double_ionisation = cross_section_of_double_ionisation
        super().__init__(electron_temperature, ion_velocity)

    @staticmethod
    def get_deuterium_mass():
        return Beam('D', 0).get_mass()

    def get_ion_velocity_normalisation_factor(self, electron_temperature):
        return numpy.sqrt(electron_temperature * scipy.constants.elementary_charge / self.deuterium_mass)

    def integrand(self, alpha, v):
        velocity = self.get_third_side_length(v * self.ion_velocity_normalisation_factor, self.ion_velocity, alpha)
        impact_energy = self.get_impact_energy(self.deuterium_mass, velocity)
        return self.maxwell(v) * velocity * \
            (self.cross_section(impact_energy) +
             2.0 * self.cross_section_of_double_ionisation(impact_energy))


class Rate1DIntegralElectron(RateNDIntegralElectron):
    @staticmethod
    def integrate(function):
        return scipy.integrate.quad(function, 0, numpy.inf)[0]

    def integrand_normalisation(self, v):
        return self.maxwell(v)


class Rate2DIntegralIon(RateNDIntegralIon):
    @staticmethod
    def integrate(function):
        return scipy.integrate.dblquad(function, 0, numpy.inf, -numpy.pi, numpy.pi)[0]

    def integrand_normalisation(self, alpha, v):
        return self.maxwell(v)


class Rate3DIntegralIon(RateNDIntegralIon):
    @staticmethod
    def integrate(function):
        return scipy.integrate.tplquad(function, 0, numpy.inf, -numpy.pi, numpy.pi, -numpy.pi/2, numpy.pi/2)[0]

    def integrand_normalisation(self, beta, alpha, v):
        return self.maxwell(v)
