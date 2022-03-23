import numpy
import scipy.interpolate

from bebim.beam import Beam
from bebim.nrl import get_nrl_rate
from bebim.integrator import RateIntegrator


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
        r = get_nrl_rate(self.species, self.ionisation_level, self.electron_temperature)
        if is_with_tabata:
            r += RateIntegrator('charge exchange',
                                self.species, self.beam_energy, self.electron_temperature, tabata_integration_dimension)\
                .get_coefficient()
        return r * self.electron_density

    def get_attenuation_beb(self, is_with_tabata=False, tabata_integration_dimension=2):
        if self.electron_density == 0:
            return 0.0
        r = RateIntegrator('electron impact ionisation', self.species, self.beam_energy, self.electron_temperature, 1)\
            .get_coefficient()
        if is_with_tabata:
            r += RateIntegrator('charge exchange',
                                self.species, self.beam_energy, self.electron_temperature, tabata_integration_dimension)\
                .get_coefficient()
        return r * self.electron_density


class RateProfile(Rate):
    def __init__(self, species, beam_energy, ionisation_level=0):
        super().__init__(species, beam_energy, ionisation_level)
        self.reference_energies = [10., 20., 50., 100., 200., 500., 1000.]
        self.nrl_spline = None
        self.beb_spline = None

    @staticmethod
    def get_spline(energy, cross_section):
        f = scipy.interpolate.interp1d(numpy.log(energy), numpy.log(cross_section), kind='cubic', fill_value='extrapolate')
        return lambda x: [0. if i == 0 else numpy.exp(f(numpy.log(i))) for i in x]

    def set_reference_energies(self, reference_energies):
        self.reference_energies = reference_energies

    def set_nrl_profile(self, is_with_tabata=False, tabata_integration_dimension=2):
        reference_rates = numpy.zeros_like(self.reference_energies, dtype=float)
        for i in range(len(reference_rates)):
            print('NRL  ' + str(int(i/len(reference_rates) * 100)) + '%', end='\r')
            self.set_profiles(self.reference_energies[i], 1)
            reference_rates[i] = self.get_attenuation_nrl(is_with_tabata, tabata_integration_dimension)
        print('NRL 100%')
        self.nrl_spline = self.get_spline(self.reference_energies, reference_rates)

    def get_nrl_profile(self, is_with_tabata=False, tabata_integration_dimension=2):
        self.set_nrl_profile(is_with_tabata, tabata_integration_dimension)
        return self.nrl_spline

    def set_beb_profile(self, is_with_tabata=False, tabata_integration_dimension=2):
        reference_rates = numpy.zeros_like(self.reference_energies, dtype=float)
        for i in range(len(reference_rates)):
            print('BEB  ' + str(int(i/len(reference_rates) * 100)) + '%', end='\r')
            self.set_profiles(self.reference_energies[i], 1)
            reference_rates[i] = self.get_attenuation_beb(is_with_tabata, tabata_integration_dimension)
        print('BEB 100%')
        self.beb_spline = self.get_spline(self.reference_energies, reference_rates)

    def get_beb_profile(self, is_with_tabata=False, tabata_integration_dimension=2):
        self.set_beb_profile(is_with_tabata, tabata_integration_dimension)
        return self.beb_spline
