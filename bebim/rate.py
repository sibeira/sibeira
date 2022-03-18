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
    reference_energies = [10, 20, 50, 100, 200, 500, 1000]

    def __init__(self, species, beam_energy, ionisation_level=0):
        super().__init__(species, beam_energy, ionisation_level)
        self.nrl_spline = None
        self.beb_spline = None

    @staticmethod
    def get_spline(energy, cross_section):
        f = scipy.interpolate.Akima1DInterpolator(numpy.log(energy), numpy.log(cross_section))
        return lambda x: numpy.exp(f(numpy.log(x)))

    def set_nrl_profile(self, is_with_tabata=False, tabata_integration_dimension=2):
        reference_rates = numpy.empty_like(self.reference_energies)
        for i in range(len(reference_rates)):
            self.set_profiles(self.reference_energies[i], 1)
            reference_rates[i] = self.get_attenuation_nrl(is_with_tabata, tabata_integration_dimension)
        self.nrl_spline = self.get_spline(self.reference_energies, reference_rates)

    def get_nrl_profile(self, is_with_tabata=False, tabata_integration_dimension=2):
        self.set_nrl_profile(is_with_tabata, tabata_integration_dimension)
        return self.nrl_spline

    def set_beb_profile(self, is_with_tabata=False, tabata_integration_dimension=2):
        reference_rates = numpy.empty_like(self.reference_energies)
        for i in range(len(reference_rates)):
            self.set_profiles(self.reference_energies[i], 1)
            reference_rates[i] = self.get_attenuation_beb(is_with_tabata, tabata_integration_dimension)
        self.beb_spline = self.get_spline(self.reference_energies, reference_rates)

    def get_beb_profile(self, is_with_tabata=False, tabata_integration_dimension=2):
        self.set_beb_profile(is_with_tabata, tabata_integration_dimension)
        return self.beb_spline
