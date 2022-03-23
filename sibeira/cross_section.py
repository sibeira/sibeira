import numpy
import scipy.interpolate
from sibeira.species import *


class CrossSection:
    def __init__(self, species, ionisation_level=0, with_Q=True):
        if ionisation_level == 0:
            self.o = OrbitalConstants(species)
        elif ionisation_level == 1:
            self.o = IonOrbitalConstants(species)
        else:
            raise 'Invalid ionisation level'

        if with_Q:
            self.Q = self.o.get('Q')
        else:
            self.Q = 1

        self.B = self.o.get('B')
        self.U = self.o.get('U')
        self.N = self.o.get('N')
        self.n = self.o.get('n')
        self.u = self.get_u()
        self.S = self.get_S()

    def get_t(self, energy):
        return energy/self.B

    def get_u(self):
        return self.U/self.B

    def get_S(self):
        a0 = 0.52918e-10
        return 4.0 * numpy.pi * a0**2 * self.N * (RYDBERG / self.B)**2

    def calculate(self, energy):
        t = self.get_t(energy)
        ln_t = numpy.log(t)
        return self.S / (t + (self.u + 1.0) / self.n) * (
            self.Q * ln_t / 2.0 * (1.0 - t**(-2)) +
            (2.0 - self.Q) * (1.0 - 1.0/t - (ln_t / (t + 1.0)))
        )

    def get_polynomial(self):
        energy = numpy.logspace(0.75, 5, 50)
        cross_section = self.calculate(energy)
        return scipy.interpolate.interp1d(energy, cross_section, fill_value="extrapolate")
