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
        self.tabata_double = CXCrossSection(1000, self.species, degree='double')
        self.tabata_double.set_polynomial()

    def set_profiles(self, electron_temperature=numpy.nan, electron_density=numpy.nan):
        if ~numpy.isnan(electron_temperature):
            self.electron_temperature = electron_temperature
        if ~numpy.isnan(electron_density):
            self.electron_density = electron_density

    def get_attenuation_nrl(self, is_with_tabata=False, tabata_integration_dimension=2):
        if self.electron_density == 0:
            return 0.0
        c = CrossSection(self.electron_temperature, self.species, self.ionisation_level)
        t = c.get_t()
        r = 1e-11 * numpy.sqrt(t) / c.B ** 1.5 / (6.0 + t) * numpy.exp(-1.0 / t)
        if is_with_tabata:
            if tabata_integration_dimension == 2:
                r += Rate2DIntegralIon(self.electron_temperature, self.speed, self.tabata, self.tabata_double).get_coefficient()
            elif tabata_integration_dimension == 3:
                r += Rate3DIntegralIon(self.electron_temperature, self.speed, self.tabata, self.tabata_double).get_coefficient()
            else:
                raise ValueError
        return r * self.electron_density

    def get_attenuation_beb(self, is_with_tabata=False, tabata_integration_dimension=2):
        if self.electron_density == 0:
            return 0.0
        r = Rate1DIntegralElectron(self.electron_temperature, self.speed, self.beb).get_coefficient()
        if is_with_tabata:
            if tabata_integration_dimension == 2:
                r += Rate2DIntegralIon(self.electron_temperature, self.speed, self.tabata, self.tabata_double).get_coefficient()
            elif tabata_integration_dimension == 3:
                r += Rate3DIntegralIon(self.electron_temperature, self.speed, self.tabata, self.tabata_double).get_coefficient()
            else:
                raise ValueError
        return r * self.electron_density


class RateNDIntegral:
    def __init__(self, electron_temperature, ion_velocity):
        self.normalisation_factor = 1
        self.deuterium_mass = self.get_deuterium_mass()
        self.ion_velocity = ion_velocity
        self.electron_velocity_normalisation = self.get_electron_velocity_normalisation_factor(electron_temperature)
        self.ion_velocity_normalisation = self.get_ion_velocity_normalisation_factor(electron_temperature)

    @staticmethod
    def get_third_side_length(a, b, alpha):
        return numpy.sqrt(a ** 2 + b ** 2 - 2 * a * b * numpy.cos(alpha))

    @staticmethod
    def maxwell(v):
        return scipy.stats.maxwell.pdf(v)

    @staticmethod
    def get_deuterium_mass():
        return Beam('D', 0).get_mass()

    @staticmethod
    def get_electron_velocity_normalisation_factor(electron_temperature):
        return numpy.sqrt(electron_temperature * scipy.constants.elementary_charge / scipy.constants.electron_mass)

    def get_ion_velocity_normalisation_factor(self, electron_temperature):
        return numpy.sqrt(electron_temperature * scipy.constants.elementary_charge / self.deuterium_mass)

    def get_coefficient(self):
        return self.integrate(self.integrand) / self.integrate(self.integrand_normalisation)


class Rate1DIntegralElectron(RateNDIntegral):
    def __init__(self, electron_temperature, ion_velocity, beb):
        super().__init__(electron_temperature, ion_velocity)
        self.beb = beb

    @staticmethod
    def integrate(function):
        return scipy.integrate.quad(function, 0, numpy.inf)[0]

    def integrand_normalisation(self, v):
        return self.maxwell(v)

    def integrand(self, v):
        velocity = v * self.electron_velocity_normalisation
        kinetic_energy = 0.5 * scipy.constants.electron_mass * velocity ** 2 / scipy.constants.elementary_charge
        return self.maxwell(v) * velocity * self.beb.f(kinetic_energy)


class Rate2DIntegralIon(RateNDIntegral):
    def __init__(self, electron_temperature, ion_velocity, tabata, tabata_double):
        self.tabata = tabata
        self.tabata_double = tabata_double
        super().__init__(electron_temperature, ion_velocity)

    @staticmethod
    def integrate(function):
        return scipy.integrate.dblquad(function, 0, numpy.inf, -numpy.pi, numpy.pi)[0]

    def integrand_normalisation(self, alpha, v):
        return self.maxwell(v)

    def integrand(self, alpha, v):
        velocity = self.get_third_side_length(v * self.ion_velocity_normalisation, self.ion_velocity, alpha)
        kinetic_energy = 0.5 * self.deuterium_mass * velocity ** 2.0 / scipy.constants.elementary_charge
        return self.maxwell(v) * velocity * \
               (self.tabata.f(kinetic_energy) + 2.0 * self.tabata_double.f(kinetic_energy))


class Rate3DIntegralIon(RateNDIntegral):
    def __init__(self, electron_temperature, ion_velocity, tabata, tabata_double):
        self.tabata = tabata
        self.tabata_double = tabata_double
        super().__init__(electron_temperature, ion_velocity)

    @staticmethod
    def integrate(function):
        return scipy.integrate.tplquad(function, 0, numpy.inf, -numpy.pi, numpy.pi, -numpy.pi/2, numpy.pi/2)[0]

    def integrand_normalisation(self, beta, alpha, v):
        return self.maxwell(v)

    def integrand(self, beta, alpha, v):
        velocity = self.get_third_side_length(v * self.ion_velocity_normalisation, self.ion_velocity, alpha)
        kinetic_energy = 0.5 * self.deuterium_mass * velocity ** 2.0 / scipy.constants.elementary_charge
        return self.maxwell(v) * velocity * \
               (self.tabata.f(kinetic_energy) + 2.0 * self.tabata_double.f(kinetic_energy))
