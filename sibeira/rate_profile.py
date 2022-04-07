import os
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
        self.tabata_spline = None

    @staticmethod
    def resolve_log_spline(f, x):
        if numpy.isscalar(x):
            return 0 if x == 0 else numpy.exp(f(numpy.log(x)))
        return [0. if i == 0 else numpy.exp(f(numpy.log(i))) for i in x]

    @staticmethod
    def get_spline(energy, cross_section):
        return scipy.interpolate.interp1d(numpy.log(energy), numpy.log(cross_section),
                                          kind='cubic', fill_value='extrapolate')

    def set_reference_energies(self, reference_energies):
        self.reference_energies = reference_energies

    def set_nrl_profile(self, tabata_integration_dimension=-1):
        reference_rates = numpy.zeros_like(self.reference_energies, dtype=float)
        for i in range(len(reference_rates)):
            print('NRL  ' + str(int(i / len(reference_rates) * 100)) + '%', end='\r')
            self.set_profiles(self.reference_energies[i])
            reference_rates[i] = self.get_full_rate_with_nrl(tabata_integration_dimension)
        print('NRL 100%')
        self.nrl_spline = self.get_spline(self.reference_energies, reference_rates)

    def get_nrl_profile(self, tabata_integration_dimension=-1):
        self.set_nrl_profile(tabata_integration_dimension)
        return self.nrl_spline

    def set_beb_profile(self, tabata_integration_dimension=-1):
        reference_rates = numpy.zeros_like(self.reference_energies, dtype=float)
        for i in range(len(reference_rates)):
            print('BEB  ' + str(int(i / len(reference_rates) * 100)) + '%', end='\r')
            self.set_profiles(self.reference_energies[i])
            reference_rates[i] = self.get_full_rate_with_beb(tabata_integration_dimension)
        print('BEB 100%')
        self.beb_spline = self.get_spline(self.reference_energies, reference_rates)

    def get_beb_profile(self, tabata_integration_dimension=-1):
        self.set_beb_profile(tabata_integration_dimension)
        return self.beb_spline

    def set_tabata_profile(self, tabata_integration_dimension=2):
        reference_rates = numpy.zeros_like(self.reference_energies, dtype=float)
        for i in range(len(reference_rates)):
            print('Tabata  ' + str(int(i / len(reference_rates) * 100)) + '%', end='\r')
            self.set_profiles(self.reference_energies[i])
            reference_rates[i] = self.get_full_rate_with_tabata(tabata_integration_dimension)
        print('Tabata 100%')
        self.tabata_spline = self.get_spline(self.reference_energies, reference_rates)

    def get_tabata_profile(self, tabata_integration_dimension=2):
        self.set_tabata_profile(tabata_integration_dimension)
        return self.tabata_spline

    def get_attenuation(self, radial_coordinates, temperatures, densities, profile_name,
                        tabata_integration_dimension=-1):
        try:
            profile = self.import_profile(profile_name, tabata_integration_dimension)
        except (FileNotFoundError, EOFError, KeyError):
            if profile_name == 'beb':
                profile = self.get_beb_profile(tabata_integration_dimension)
            elif profile_name == 'nrl':
                profile = self.get_nrl_profile(tabata_integration_dimension)
            elif profile_name == 'tabata':
                profile = self.get_tabata_profile(tabata_integration_dimension)
            else:
                raise (ValueError('Invalid profile: ' + profile_name))
            self.export_profile(profile_name, tabata_integration_dimension, profile)
        rate = self.resolve_log_spline(profile, temperatures) * densities / self.speed
        return numpy.exp(scipy.integrate.cumulative_trapezoid(rate, radial_coordinates, initial=0))

    def import_profile(self, profile_name, tabata_integration_dimension, destination_directory='data'):
        try:
            path = self.get_file_name(destination_directory)
            beam_energy_as_string = self.get_beam_energy_as_string()
            dimension_as_string = str(tabata_integration_dimension)
            profile_database = numpy.load(path, allow_pickle=True).item()
            return profile_database[beam_energy_as_string][profile_name][dimension_as_string]
        except FileNotFoundError:
            raise (FileNotFoundError('There is no profile for ' + self.species))
        except KeyError:
            raise (KeyError('The profile is not found: ' + profile_name +
                            ' (Tabata ' + (str(tabata_integration_dimension) + 'D)'
                                           if tabata_integration_dimension >= 0 else 'OFF)')))

    def export_profile(self, profile_name, tabata_integration_dimension, profile, destination_directory='data'):
        path = self.get_file_name(destination_directory)
        try:
            profile_database = numpy.load(path, allow_pickle=True).item()
        except (FileNotFoundError, EOFError):
            profile_database = {}
        beam_energy_as_string = self.get_beam_energy_as_string()
        dimension_as_string = str(tabata_integration_dimension)
        self.add_to_database(profile_database, profile, beam_energy_as_string, dimension_as_string, profile_name)
        if not (os.path.exists(destination_directory)):
            os.mkdir(destination_directory)
        numpy.save(path, profile_database)

    @staticmethod
    def add_to_database(profile_database, profile, beam_energy_as_string, dimension_as_string, profile_name):
        try:
            profile_database[beam_energy_as_string]
        except KeyError:
            profile_database[beam_energy_as_string] = {}
        try:
            profile_database[beam_energy_as_string][profile_name]
        except KeyError:
            profile_database[beam_energy_as_string][profile_name] = {}
        profile_database[beam_energy_as_string][profile_name][dimension_as_string] = profile

    def get_beam_energy_as_string(self):
        return str(self.beam_energy / 1000.)

    def get_file_name(self, destination):
        return destination + '/' + self.species + '.npy'
