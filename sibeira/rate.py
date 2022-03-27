import numpy

from sibeira.beam import Beam
from nrl.rate import get_nrl_rate
from sibeira.integrator import RateIntegrator


class Rate(Beam):
    def __init__(self, species, beam_energy, ionisation_level=0):
        super().__init__(species, beam_energy, ionisation_level)

    def set_profiles(self, electron_temperature=numpy.nan):
        if ~numpy.isnan(electron_temperature):
            self.electron_temperature = electron_temperature

    def get_full_rate_with_nrl(self, is_with_tabata=False, tabata_integration_dimension=2):
        r = get_nrl_rate(self.species, self.ionisation_level, self.electron_temperature)
        if is_with_tabata:
            r += RateIntegrator('charge exchange',
                                self.species, self.beam_energy, self.electron_temperature, tabata_integration_dimension)\
                .get_coefficient()
        return r

    def get_full_rate_with_beb(self, is_with_tabata=False, tabata_integration_dimension=2):
        r = RateIntegrator('electron impact ionisation', self.species, self.beam_energy, self.electron_temperature, 1)\
            .get_coefficient()
        if is_with_tabata:
            r += RateIntegrator('charge exchange',
                                self.species, self.beam_energy, self.electron_temperature, tabata_integration_dimension)\
                .get_coefficient()
        return r

    def get_full_rate_with_tabata(self, tabata_integration_dimension=2):
        return RateIntegrator('charge exchange',
                              self.species, self.beam_energy, self.electron_temperature, tabata_integration_dimension)\
            .get_coefficient()
