import numpy
from bebim.beam import Beam
from bebim.cross_section import CrossSection
from integrator import RateIntegrator


class Rate(Beam):
    def __init__(self, species, beam_energy, ionisation_level=0):
        super().__init__(species, beam_energy, ionisation_level)

    def set_profiles(self, electron_temperature=numpy.nan, electron_density=numpy.nan):
        if ~numpy.isnan(electron_temperature):
            self.electron_temperature = electron_temperature
        if ~numpy.isnan(electron_density):
            self.electron_density = electron_density

    def get_attenuation_nrl(self, is_with_tabata=False, tabata_integration_dimension=2):
        if self.electron_density == 0:
            return 0.0
        c = CrossSection(self.species, self.ionisation_level)
        t = c.get_t(self.electron_temperature)
        r = 1e-11 * numpy.sqrt(t) / c.B ** 1.5 / (6.0 + t) * numpy.exp(-1.0 / t)
        if is_with_tabata:
            r += RateIntegrator('charge exchange',
                                self.species, self.speed, self.electron_temperature, tabata_integration_dimension)\
                .get_coefficient()
        return r * self.electron_density

    def get_attenuation_beb(self, is_with_tabata=False, tabata_integration_dimension=2):
        if self.electron_density == 0:
            return 0.0
        r = RateIntegrator('electron impact ionisation', self.species, self.speed, self.electron_temperature, 1)\
            .get_coefficient()
        if is_with_tabata:
            r += RateIntegrator('charge exchange',
                                self.species, self.speed, self.electron_temperature, tabata_integration_dimension)\
                .get_coefficient()
        return r * self.electron_density
