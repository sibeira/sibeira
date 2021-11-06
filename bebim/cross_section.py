import numpy
import scipy.interpolate
from species import *


class CrossSection:
    def __init__(self, energy, species, ionisation_level=0, with_Q=True):
        self.energy = energy
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
        self.t = self.get_t()
        self.u = self.get_u()
        self.S = self.get_S()
        self.f = lambda: None

    def get_t(self):
        return self.energy/self.B

    def get_u(self):
        return self.U/self.B

    def get_S(self):
        a0 = 0.52918e-10
        return 4.0 * numpy.pi * a0**2 * self.N * (RYDBERG / self.B)**2

    def calculate(self):
        ln_t = numpy.log(self.t)
        return self.S / (self.t + (self.u + 1.0) / self.n) * (
            self.Q * ln_t / 2.0 * (1.0 - self.t**(-2)) +
            (2.0 - self.Q) * (1.0 - 1.0/self.t - (ln_t / (self.t + 1.0)))
        )

    def set_polynomial(self):
        energy = numpy.logspace(0.75, 5, 50)
        self.energy = energy
        self.t = self.get_t()
        cross_section = self.calculate()
        self.f = scipy.interpolate.interp1d(energy, cross_section, fill_value="extrapolate")

