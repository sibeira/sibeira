import os
import numpy

from rate import Rate


class RateProfileIO(Rate):
    default_destination_directory = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')

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

    def export_profile(self, profile_name, tabata_integration_dimension, profile,
                       destination_directory=default_destination_directory):
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

    def get_beam_energy_as_string(self):
        return str(self.beam_energy / 1000.)

    def get_file_name(self, destination):
        return destination + '/' + self.species + '.npy'
