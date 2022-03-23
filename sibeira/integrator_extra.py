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

    def integrand_normalisation_1d(self, velocity):
        return self.maxwell(velocity)
