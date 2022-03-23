import numpy
import scipy.constants
import scipy.integrate
import scipy.stats

from sibeira.integrator import RateIntegrator


def maxwell(m, T, v):
    m__2kT = m / (2.0 * T * scipy.constants.elementary_charge)
    return numpy.sqrt((m__2kT / numpy.pi) ** 3) * numpy.exp(- m__2kT * v ** 2)


class RateIntegratorExtra(RateIntegrator):
    def maxwell(self, velocity):
        return maxwell(self.target_mass, self.temperature, velocity)

    def integrand_1d(self, velocity):
        impact_energy = self.get_impact_energy(velocity)
        return self.maxwell(velocity) * velocity * self.cross_section(impact_energy)

    def integrand_2d(self, alpha, v):
        velocity = self.get_third_side_length(v, self.beam_speed, alpha)
        impact_energy = self.get_impact_energy(velocity)
        return self.maxwell(v) * velocity * self.cross_section(impact_energy)

    def integrand_3d(self, beta, alpha, v):
        velocity = self.get_third_side_length(v, self.beam_speed, alpha, beta)
        impact_energy = self.get_impact_energy(velocity)
        return self.maxwell(v) * velocity * self.cross_section(impact_energy)
