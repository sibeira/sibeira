import sys
import os
import matplotlib.pyplot
import numpy
import scipy.integrate

from rate import Rate

sys.path.append(os.environ['RENATE_OD'])

from manager import RenateODManager
from beamlet import set_beamlet, BeamletGeometry
from profiles import Profiles
from utils import *


def get_export_name(mode):
    return 'atomic' + ('' if mode == 'plot' else '_' + mode)


def test_export_name():
    print(get_export_name('plot'))
    print(get_export_name('log'))
    print(get_export_name('lorem_ipsum'))


def plot_attenuation_profile(shot_number, time, species, energy, radial_coordinate,
                             relative_attenuation_from_renate_od,
                             relative_attenuation_from_renate_od_electron,
                             relative_attenuation_from_beb,
                             relative_attenuation_from_nrl,
                             relative_attenuation_from_beb_tabata,
                             relative_attenuation_from_nrl_tabata,
                             mode='plot'):
    fig, ax = matplotlib.pyplot.subplots()
    fig.set_size_inches(5, 2)
    if mode == 'log':
        plot = getattr(ax, 'semilogy')
    else:
        plot = getattr(ax, mode)
    export_name = get_export_name(mode)

    plot(radial_coordinate, relative_attenuation_from_renate_od,
         '-', linewidth=2, facecolor='tab:blue', label='RENATE-OD')
    plot(radial_coordinate, relative_attenuation_from_beb_tabata,
         '-', linewidth=1.5, facecolor='tab:red', label='BEB+Tabata 3D')
    plot(radial_coordinate, relative_attenuation_from_nrl_tabata,
         '-', linewidth=1.5, facecolor='tab:orange', label='NRL+Tabata 3D')
    plot(radial_coordinate, relative_attenuation_from_renate_od_electron,
         '--', linewidth=1, facecolor='tab:blue', label='RENATE-OD (no electron)')
    plot(radial_coordinate, relative_attenuation_from_beb,
         '--', linewidth=1, facecolor='tab:red', label='BEB')
    plot(radial_coordinate, relative_attenuation_from_nrl,
         '--', linewidth=1.5, facecolor='tab:orange', label='NRL')
    ax.legend()

    matplotlib.pyplot.xlim(0.6, 0.7399)
    matplotlib.pyplot.minorticks_on()
    matplotlib.pyplot.grid(which='major')
    matplotlib.pyplot.xlabel('$R$ [m]', labelpad=-10.5, loc='right')
    matplotlib.pyplot.ylabel('neutral beam attenuation')
    matplotlib.pyplot.title('COMPASS #' + shot_number + ' (' + time + ' ms, '+species + ', ' + energy + ' keV)')

    matplotlib.pyplot.savefig(export_name + '.png')
    matplotlib.pyplot.savefig(export_name + '.pdf')
    matplotlib.pyplot.show()


def plot_attenuation_comparison(shot_number='17178', time='1097', species='Li', energy='80'):
    z = 0
    tor = 0
    beamlet_geometry = set_beamlet(z, tor)

    r = RenateODManager(beamlet_geometry, shot_number, time, species, energy, 'just electron')
    radial_coordinates, relative_attenuation_rod = r.get_attenuation_profile()
    relative_attenuation_rod_no_ion = relative_attenuation_rod.fillna(0)

    r = RenateODManager(beamlet_geometry, shot_number, time, species, energy)
    radial_coordinates, relative_attenuation_rod = r.get_attenuation_profile()
    relative_attenuation_rod = relative_attenuation_rod.fillna(0)

    p = Profiles(is_export=False)
    if not numpy.array_equal(radial_coordinates.to_numpy, beamlet_geometry.rad):
        assert ValueError
    temperatures = p.get_temperature()
    densities = p.get_density()

    rate_beb = numpy.zeros_like(beamlet_geometry.rad)
    rate_nrl = numpy.zeros_like(beamlet_geometry.rad)
    rate_beb_tabata = numpy.zeros_like(beamlet_geometry.rad)
    rate_nrl_tabata = numpy.zeros_like(beamlet_geometry.rad)

    rate = Rate(species=species, beam_energy=float(energy) * 1000.0)
    v_beam = rate.get_speed()

    for i in range(beamlet_geometry.rad.size):
        print(i)
        rate.set_profiles(electron_temperature=temperatures[i], electron_density=densities[i])
        rate_beb[i] = rate.get_attenuation_beb() / v_beam
        rate_nrl[i] = rate.get_attenuation_nrl() / v_beam
        rate_beb_tabata[i] = rate.get_attenuation_beb(is_with_tabata=True, tabata_integration_dimension=2) / v_beam
        rate_nrl_tabata[i] = rate.get_attenuation_nrl(is_with_tabata=True, tabata_integration_dimension=2) / v_beam

    relative_attenuation_from_beb = numpy.exp(scipy.integrate.cumulative_trapezoid(rate_beb, beamlet_geometry.rad, initial=0))
    relative_attenuation_from_nrl = numpy.exp(scipy.integrate.cumulative_trapezoid(rate_nrl, beamlet_geometry.rad, initial=0))
    relative_attenuation_from_beb_tabata =\
        numpy.exp(scipy.integrate.cumulative_trapezoid(rate_beb_tabata, beamlet_geometry.rad, initial=0))
    relative_attenuation_from_nrl_tabata = \
        numpy.exp(scipy.integrate.cumulative_trapezoid(rate_nrl_tabata, beamlet_geometry.rad, initial=0))

    plot_attenuation_profile(shot_number, time, species, energy, radial_coordinates,
                             relative_attenuation_rod, relative_attenuation_rod_no_ion,
                             relative_attenuation_from_beb, relative_attenuation_from_nrl,
                             relative_attenuation_from_beb_tabata, relative_attenuation_from_nrl_tabata)
    plot_attenuation_profile(shot_number, time, species, energy, radial_coordinates,
                             relative_attenuation_rod, relative_attenuation_rod_no_ion,
                             relative_attenuation_from_beb, relative_attenuation_from_nrl,
                             relative_attenuation_from_beb_tabata, relative_attenuation_from_nrl_tabata,
                             mode='log')


if __name__ == "__main__":
    a_shot_number = '17178'
    a_time = '1097'
    plot_attenuation_comparison(shot_number=a_shot_number, time=a_time, species='Na')
