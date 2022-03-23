import numpy
import scipy.interpolate
import scipy.integrate

from sibeira.rate import Rate


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

    def get_attenuation(self, major_radius, temperatures, densities, profile_name,
                        is_with_tabata=False, tabata_integration_dimension=2):
        if profile_name == 'beb':
            profile = self.get_beb_profile(is_with_tabata, tabata_integration_dimension)
        elif profile_name == 'nrl':
            profile = self.get_nrl_profile(is_with_tabata, tabata_integration_dimension)
        else:
            raise(ValueError('Invalid profile: ' + profile_name))
        rate = profile(temperatures) * densities / self.get_speed()
        return numpy.exp(scipy.integrate.cumulative_trapezoid(rate, major_radius, initial=0))
