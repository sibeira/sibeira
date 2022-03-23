import numpy
import scipy.constants
import scipy.integrate
import scipy.stats

from bebim.beam import Beam
from bebim.cross_section import CrossSection
from tabata_ctf.cross_section import CrossSection as CXCrossSection


class RateIntegrator:
    def __init__(self, reaction_name, beam_species, beam_energy, temperature, dimension=2):
        self.reaction_name = reaction_name
        self.temperature = temperature
        self.beam_species = beam_species
        self.beam_energy = beam_energy
        self.beam_speed = self.get_projectile_velocity()
        self.target_mass = self.get_target_mass()
        self.cross_section = self.get_cross_section()
        self.normalisation_factor = 1
        self.maxwell_normalisation_factor = self.get_maxwell_normalisation_factor(temperature)
        self.set_integrator(dimension)

    def get_coefficient(self):
        return self.integrate(self.integrand) / self.integrate(self.integrand_normalisation)

    @staticmethod
    def get_third_side_length(a, b, alpha, beta=0):
        return numpy.sqrt(a ** 2 + b ** 2 - 2 * a * b * numpy.cos(alpha) * numpy.cos(beta))

    @staticmethod
    def maxwell(v):
        return scipy.stats.maxwell.pdf(v)

    def get_maxwell_normalisation_factor(self, temperature):
        return numpy.sqrt(temperature * scipy.constants.elementary_charge / self.target_mass)

    def get_impact_energy(self, velocity):
        return 0.5 * self.target_mass * velocity ** 2 / scipy.constants.elementary_charge

    def get_projectile_velocity(self):
        return Beam(self.beam_species, self.beam_energy).get_speed()

    def get_target_mass(self):
        if self.reaction_name == 'electron impact ionisation':
            return scipy.constants.electron_mass
        elif self.reaction_name == 'charge exchange':
            return Beam('D', 0).get_mass()
        else:
            raise ValueError('The ionisation reaction is unknown: ' + self.reaction_name)

    def get_cross_section(self):
        if self.reaction_name == 'electron impact ionisation':
            beb = CrossSection(self.beam_species)
            return beb.get_polynomial()
        elif self.reaction_name == 'charge exchange':
            tabata = CXCrossSection(self.beam_species)
            return tabata.get_polynomial()
            #tabata_double = CXCrossSection(self.beam_species, degree='double')
            #return lambda x: tabata.get_polynomial()(x) + 2.0 * tabata_double.get_polynomial()(x)
        else:
            raise ValueError('The ionisation reaction is unknown: ' + self.reaction_name)

    def set_integrator(self, dimension):
        if dimension == 1:
            self.integrand = self.integrand_1d
            self.integrand_normalisation = self.integrand_normalisation_1d
            self.integrate = self.integrate_1d
        elif dimension == 2:
            self.integrand = self.integrand_2d
            self.integrand_normalisation = self.integrand_normalisation_2d
            self.integrate = self.integrate_2d
        elif dimension == 3:
            self.integrand = self.integrand_3d
            self.integrand_normalisation = self.integrand_normalisation_3d
            self.integrate = self.integrate_3d
        else:
            raise ValueError('Invalid integration dimension: ' + str(dimension))

    def integrand_all(self, v, alpha=0, beta=0):
        velocity = self.get_third_side_length(v * self.maxwell_normalisation_factor, self.beam_speed, alpha)
        impact_energy = self.get_impact_energy(velocity)
        return self.maxwell(v) * velocity * self.cross_section(impact_energy)

    def integrand_normalisation_all(self, v, alpha=0, beta=0):
        return self.maxwell(v)

    def integrand_1d(self, v):
        velocity = v * self.maxwell_normalisation_factor
        impact_energy = self.get_impact_energy(velocity)
        return self.maxwell(v) * velocity * self.cross_section(impact_energy)

    @staticmethod
    def integrate_1d(function):
        return scipy.integrate.quad(function, 0, numpy.inf)[0]

    def integrand_normalisation_1d(self, v):
        return self.maxwell(v)

    def integrand_2d(self, alpha, v):
        velocity = self.get_third_side_length(v * self.maxwell_normalisation_factor, self.beam_speed, alpha)
        impact_energy = self.get_impact_energy(velocity)
        return self.maxwell(v) * velocity * self.cross_section(impact_energy)

    @staticmethod
    def integrate_2d(function):
        return scipy.integrate.dblquad(function, 0, numpy.inf, -numpy.pi, numpy.pi)[0]

    def integrand_normalisation_2d(self, alpha, v):
        return self.maxwell(v)

    @staticmethod
    def integrate_3d(function):
        return scipy.integrate.tplquad(function, 0, numpy.inf, -numpy.pi, numpy.pi, -numpy.pi/2, numpy.pi/2)[0]

    def integrand_normalisation_3d(self, beta, alpha, v):
        return self.maxwell(v)

    def integrand_3d(self, beta, alpha, v):
        velocity = self.get_third_side_length(v * self.maxwell_normalisation_factor, self.beam_speed, alpha, beta)
        impact_energy = self.get_impact_energy(velocity)
        return self.maxwell(v) * velocity * self.cross_section(impact_energy)
