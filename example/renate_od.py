import sys
import os
import matplotlib.pyplot
import numpy
import scipy.integrate

from rate import RateProfile

sys.path.append(os.environ['RENATE_OD'])

from manager import RenateODManager
from beamlet import set_beamlet, BeamletGeometry
from profiles import Profiles
from utils import *


def get_export_name(mode, shot_number, time, species, energy, dimension=2):
    directory_path = os.path.dirname(os.path.abspath(__file__))
    return directory_path + '/figs/' + species + '_' + str(energy) + 'keV_' + str(shot_number) + '_' + str(time) + \
           ('' if mode == 'plot' else '_' + mode) + ('_3d' if dimension == 3 else '')


def test_export_name():
    print(get_export_name('plot', 12345, 678, 'test', 123))
    print(get_export_name('log', 12345, 678, 'test', 123))
    print(get_export_name('lorem_ipsum', 12345, 678, 'test', 123))
    print(get_export_name('dolor_sit_amet', 12345, 678, 'test', 123, 3))


def plot_attenuation_profile(shot_number, time, species, energy, dimension, radial_coordinate,
                             relative_attenuation_from_renate_od,
                             relative_attenuation_from_renate_od_electron,
                             relative_attenuation_from_beb,
                             relative_attenuation_from_nrl,
                             relative_attenuation_from_beb_tabata,
                             relative_attenuation_from_nrl_tabata,
                             mode='plot'):
    fig, ax = matplotlib.pyplot.subplots()
    fig.set_size_inches(6, 2)
    if mode == 'log':
        plot = getattr(ax, 'semilogy')
    else:
        plot = getattr(ax, mode)
    export_name = get_export_name(mode, shot_number, time, species, energy, dimension)

    plot(radial_coordinate, relative_attenuation_from_renate_od,
         '-', linewidth=2, color='tab:blue', label='RENATE-OD')
    plot(radial_coordinate, relative_attenuation_from_beb_tabata,
         '-', linewidth=1.5, color='tab:red', label='BEB+Tabata 3D')
    plot(radial_coordinate, relative_attenuation_from_nrl_tabata,
         '-', linewidth=1.5, color='tab:orange', label='NRL+Tabata 3D')
    plot(radial_coordinate, relative_attenuation_from_renate_od_electron,
         '--', linewidth=1, color='tab:blue', label='RENATE-OD ($n_i=0$)')
    plot(radial_coordinate, relative_attenuation_from_beb,
         '--', linewidth=1, color='tab:red', label='BEB')
    plot(radial_coordinate, relative_attenuation_from_nrl,
         '--', linewidth=1.5, color='tab:orange', label='NRL')
    ax.legend(bbox_to_anchor=(1.0, 0.5), loc="center left", borderaxespad=0, frameon=False)

    matplotlib.pyplot.xlim(0.6, 0.7399)
    matplotlib.pyplot.minorticks_on()
    matplotlib.pyplot.grid(which='major')
    matplotlib.pyplot.xlabel('$R$ [m]')
    ax.xaxis.set_label_coords(1.0, -0.055)
    matplotlib.pyplot.ylabel('neutral beam attenuation')
    matplotlib.pyplot.title('COMPASS #' + str(shot_number) + ' (' + str(time) + ' ms, '
                            + species + ', ' + str(energy) + ' keV)')

    matplotlib.pyplot.subplots_adjust(right=0.7)
    matplotlib.pyplot.savefig(export_name + '.png')
    matplotlib.pyplot.savefig(export_name + '.pdf')
    matplotlib.pyplot.show()


def test_plot_attenuation_profile():
    r = numpy.linspace(0.6, 0.74, 11)
    p1 = numpy.random.random_sample(r.shape)
    p2 = numpy.random.random_sample(r.shape)
    p3 = numpy.random.random_sample(r.shape)
    p4 = numpy.random.random_sample(r.shape)
    p5 = numpy.random.random_sample(r.shape)
    p6 = numpy.random.random_sample(r.shape)
    plot_attenuation_profile(12345, 678, 'test', 123, 3, r, p1, p2, p3, p4, p5, p6)


def run_attenuation_comparison(shot_number, time, species, energy, dimension=2):
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

    rate = RateProfile(species=species, beam_energy=float(energy) * 1000.0)
    v_beam = rate.get_speed()

    rate_beb = rate.get_beb_profile()(temperatures) * densities / v_beam
    rate_nrl = rate.get_nrl_profile()(temperatures) * densities / v_beam
    rate_beb_tabata = rate.get_beb_profile(True, 2)(temperatures) * densities / v_beam
    rate_nrl_tabata = rate.get_nrl_profile(True, 2)(temperatures) * densities / v_beam

    relative_attenuation_from_beb = numpy.exp(
        scipy.integrate.cumulative_trapezoid(rate_beb, beamlet_geometry.rad, initial=0))
    relative_attenuation_from_nrl = numpy.exp(
        scipy.integrate.cumulative_trapezoid(rate_nrl, beamlet_geometry.rad, initial=0))
    relative_attenuation_from_beb_tabata = \
        numpy.exp(scipy.integrate.cumulative_trapezoid(rate_beb_tabata, beamlet_geometry.rad, initial=0))
    relative_attenuation_from_nrl_tabata = \
        numpy.exp(scipy.integrate.cumulative_trapezoid(rate_nrl_tabata, beamlet_geometry.rad, initial=0))

    plot_attenuation_profile(shot_number, time, species, energy, dimension, radial_coordinates,
                             relative_attenuation_rod, relative_attenuation_rod_no_ion,
                             relative_attenuation_from_beb, relative_attenuation_from_nrl,
                             relative_attenuation_from_beb_tabata, relative_attenuation_from_nrl_tabata)
    plot_attenuation_profile(shot_number, time, species, energy, dimension, radial_coordinates,
                             relative_attenuation_rod, relative_attenuation_rod_no_ion,
                             relative_attenuation_from_beb, relative_attenuation_from_nrl,
                             relative_attenuation_from_beb_tabata, relative_attenuation_from_nrl_tabata,
                             mode='log')


if __name__ == "__main__":
#    run_attenuation_comparison(shot_number='17178', time='1097', species='Li', energy='40', dimension=3)
#    run_attenuation_comparison(shot_number='17178', time='1097', species='Na', energy='80', dimension=3)
    run_attenuation_comparison(shot_number='17178', time='1097', species='Li', energy='40')
    run_attenuation_comparison(shot_number='17178', time='1097', species='Li', energy='60')
    run_attenuation_comparison(shot_number='17178', time='1097', species='Li', energy='80')
    run_attenuation_comparison(shot_number='17178', time='1097', species='Na', energy='40')
    run_attenuation_comparison(shot_number='17178', time='1097', species='Na', energy='60')
    run_attenuation_comparison(shot_number='17178', time='1097', species='Na', energy='80')
