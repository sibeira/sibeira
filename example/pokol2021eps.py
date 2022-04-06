import numpy
import pandas
import matplotlib.pyplot

import sys
sys.path.append('/media/matyi/home2/Documents/git/taiga-project/taiga/preproc/renate_od/')

from beam import Beam
from beamlet import set_beamlet, BeamletGeometry
from manager import RenateODManager


class RenateODManagerMock(RenateODManager):
    def __init__(self, beamlet_geometry, temperature, density, species, energy):
        self.beamlet_geometry = beamlet_geometry
        self.species = species
        self.energy = str(energy)
        self.temperature = temperature
        self.density = density
        self.beamlet = self.get_beamlet()
        self.relative_attenuation = self.get_relative_attenuation()

    @staticmethod
    def get_components():
        return pandas.DataFrame(
            {'q': [-1, 1], 'Z': [0, 1], 'A': [0, 2]},
            index=['electron', 'ion'])

    def get_profiles(self):
        tuples = [('beamlet grid', 'distance', 'm'),
                  ('electron', 'density', 'm-3'),
                  ('electron', 'temperature', 'eV'),
                  ('ion1', 'density', 'm-3'),
                  ('ion1', 'temperature', 'eV')]

        header = pandas.MultiIndex.from_tuples(tuples, names=['type', 'property', 'unit'])
        distance = self.beamlet_geometry.rad
        density = numpy.full_like(distance, self.density)
        temperature = numpy.full_like(distance, self.temperature)

        profiles_data = numpy.transpose(numpy.array([distance, density, temperature, density, temperature]))
        return pandas.DataFrame(profiles_data, columns=header)


def plot_attenuation(ax, radial_coordinate, relative_attenuation, label):

    log_n = numpy.log10(relative_attenuation/relative_attenuation[0])
    ax.plot(radial_coordinate, log_n, '-', linewidth=2, label=label)


def run_test_eps2021(temperature, species='H', energy=30,  density=1e19):
    beamlet_geometry = BeamletGeometry()
    beamlet_geometry.rad = numpy.linspace(0.0, 2.0, 30)
    beamlet_geometry.set_with_value(0, 'z', 'rad')
    beamlet_geometry.set_with_value(0, 'tor', 'rad')

    b = Beam(species=species, beam_energy=float(energy) * 1000.0)


    b.set_profiles(electron_temperature=temperature, electron_density=density)
    v_beam = b.get_speed()
    timestep = (beamlet_geometry.rad[1] - beamlet_geometry.rad[0]) / v_beam

    relative_attenuation_1d = numpy.zeros_like(beamlet_geometry.rad)
    relative_attenuation_nrl = numpy.zeros_like(beamlet_geometry.rad)
    relative_attenuation_nrl0 = numpy.zeros_like(beamlet_geometry.rad)

    a0 = b.get_attenuation_nrl(is_with_tabata=False)
    #a1 = b.get_attenuation_beb()
    a1 = b.get_attenuation_nrl(is_with_tabata=True)
    a2 = b.get_attenuation_beb_tabata_3d()

    for i in range(beamlet_geometry.rad.size):
        if i == 0:
            relative_attenuation_nrl0[i] = (1.0 - a0 * timestep)
            relative_attenuation_nrl[i] = (1.0 - a1 * timestep)
            relative_attenuation_1d[i] = (1.0 - a2 * timestep)
        else:
            relative_attenuation_nrl0[i] = (1.0 - a0 * timestep) * relative_attenuation_nrl0[i - 1]
            relative_attenuation_nrl[i] = (1.0 - a1 * timestep) * relative_attenuation_nrl[i - 1]
            relative_attenuation_1d[i] = (1.0 - a2 * timestep) * relative_attenuation_1d[i - 1]


    r = RenateODManagerMock(beamlet_geometry, temperature, density, species, energy)
    radial_coordinates, relative_attenuation_rod = r.get_attenuation_profile()
    #relative_attenuation_rod = relative_attenuation_rod.fillna(0)

    fig, ax = matplotlib.pyplot.subplots()
    fig.set_size_inches(5, 2)

    plot_attenuation(ax, radial_coordinates, relative_attenuation_rod, 'ROD')
    plot_attenuation(ax, radial_coordinates, relative_attenuation_1d, 'BEBIM')
    plot_attenuation(ax, radial_coordinates, relative_attenuation_nrl, 'NRL+T')
    plot_attenuation(ax, radial_coordinates, relative_attenuation_nrl0, 'NRL')
    matplotlib.pyplot.grid('minor')
    ax.legend()
    matplotlib.pyplot.title('T=' + str(temperature/1000.0) + ' keV, ' + species + ' @ ' + str(energy) + ' keV')
    matplotlib.pyplot.show()


if __name__ == "__main__":
    run_test_eps2021(100, 'Na')
    run_test_eps2021(1000, 'Na')
    run_test_eps2021(20000, 'Na')
    run_test_eps2021(100, 'Li')
    run_test_eps2021(1000, 'Li')
    run_test_eps2021(20000, 'Li')
