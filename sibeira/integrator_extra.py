import numpy
import scipy.constants
import scipy.integrate
import scipy.stats

from sibeira.integrator import RateIntegrator


class RateIntegratorExtra(RateIntegrator):
    @staticmethod
    def maxwell_v2(m, T, v):
        m__2kT = m / (2.0 * T * scipy.constants.elementary_charge)
        return numpy.sqrt((m__2kT / numpy.pi) ** 3) * numpy.exp(- m__2kT * v ** 2)

    def integrand_1d(self, velocity):
        impact_energy = self.get_impact_energy(velocity)
        return self.maxwell_v2(velocity) * velocity * self.cross_section(impact_energy)

    def integrand_normalisation_1d(self, v):
        return self.maxwell_v2(v)
