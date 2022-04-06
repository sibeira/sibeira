import sys
import os
import matplotlib.pyplot
import numpy


def get_export_name(mode, shot_number, time, species, energy, dimension=2, scenario='total'):
    directory_path = os.path.dirname(os.path.abspath(__file__))
    return directory_path + '/figures/' + species + '_' + str(energy) + 'keV_' + str(shot_number) + '_' + str(time) + \
           get_scenario_path(scenario) + \
           ('' if mode == 'plot' else '_' + mode) + \
           ('_3d' if dimension == 3 else '')


def get_title_name(shot_number, time, species, energy, scenario='total'):
    return 'COMPASS #' + str(shot_number) + ' (' + str(time) + ' ms, ' + species + ', ' + str(energy) + ' keV' + \
        get_scenario_latex(scenario) + ')'


def get_scenario_latex(scenario):
    if scenario == 'total':
        return ''
    elif scenario == 'just electron':
        return ', $n_i = 0$'
    elif scenario == 'just ion':
        return ', $n_e = 0$'
    else:
        raise(ValueError('Invalid scenario: ' + scenario))


def get_scenario_path(scenario):
    if scenario == 'total':
        return ''
    elif scenario == 'just electron':
        return '_electron'
    elif scenario == 'just ion':
        return '_ion'
    else:
        raise (ValueError('Invalid scenario: ' + scenario))


def plot_attenuation_profile(shot_number, time, species, energy, dimension, radial_coordinate,
                             relative_attenuation_profiles, profile_names,
                             mode='plot', scenario='total'):
    fig, ax = matplotlib.pyplot.subplots()
    fig.set_size_inches(6, 2)
    if mode == 'log':
        plot = getattr(ax, 'semilogy')
    else:
        plot = getattr(ax, mode)
    export_name = get_export_name(mode, shot_number, time, species, energy, dimension, scenario)
    title_name = get_title_name(shot_number, time, species, energy, scenario)
    plot(radial_coordinate, relative_attenuation_profiles[0],
         linewidth=2, label=profile_names[0])

    for i in range(1, len(relative_attenuation_profiles)):
        plot(radial_coordinate, relative_attenuation_profiles[i],
             linewidth=1.5, label=profile_names[i])
    ax.legend(bbox_to_anchor=(1.0, 0.5), loc="center left", borderaxespad=0, frameon=False)

    matplotlib.pyplot.xlim(0.6, 0.7399)
    matplotlib.pyplot.minorticks_on()
    matplotlib.pyplot.grid(which='major')
    matplotlib.pyplot.xlabel('$R$ [m]')
    ax.xaxis.set_label_coords(1.0, -0.055)
    matplotlib.pyplot.ylabel('neutral beam attenuation')
    matplotlib.pyplot.title(title_name)

    matplotlib.pyplot.subplots_adjust(right=0.7)
    matplotlib.pyplot.savefig(export_name + '.png')
    matplotlib.pyplot.savefig(export_name + '.pdf')
    matplotlib.pyplot.show()


def get_renate_od_attenuation_profile(beamlet_geometry, shot_number, time, species, energy, scenario='default'):
    sys.path.append(os.environ['RENATE_OD'])
    from manager import RenateODManager

    r = RenateODManager(beamlet_geometry, shot_number, time, species, energy, scenario)
    radial_coordinates, relative_attenuation_rod = r.get_attenuation_profile()
    return radial_coordinates, relative_attenuation_rod.fillna(0)


def run_attenuation_comparison(shot_number, time, species, energy, dimension=2):
    from sibeira.rate_profile import RateProfile

    sys.path.append(os.environ['RENATE_OD'])
    from beamlet import set_beamlet
    from profiles import Profiles

    z = 0
    tor = 0
    beamlet_geometry = set_beamlet(z, tor)

    radial_coordinates, relative_attenuation_rod_just_electron = \
        get_renate_od_attenuation_profile(beamlet_geometry, shot_number, time, species, energy, 'just electron')
    radial_coordinates, relative_attenuation_rod_just_ion = \
        get_renate_od_attenuation_profile(beamlet_geometry, shot_number, time, species, energy, 'just ion')
    radial_coordinates, relative_attenuation_rod =\
        get_renate_od_attenuation_profile(beamlet_geometry, shot_number, time, species, energy)

    p = Profiles(is_export=False)
    if not numpy.array_equal(radial_coordinates.to_numpy, beamlet_geometry.rad):
        assert ValueError
    temperatures = p.get_temperature()
    densities = p.get_density()

    rate = RateProfile(species=species, beam_energy=float(energy) * 1000.0)

    relative_attenuation_from_beb = rate.get_attenuation(beamlet_geometry.rad, temperatures, densities, 'beb')
    relative_attenuation_from_nrl = rate.get_attenuation(beamlet_geometry.rad, temperatures, densities, 'nrl')
    relative_attenuation_from_beb_tabata = rate.get_attenuation(beamlet_geometry.rad, temperatures, densities, 'beb', True, 2)
    relative_attenuation_from_nrl_tabata = rate.get_attenuation(beamlet_geometry.rad, temperatures, densities, 'nrl', True, 2)
    relative_attenuation_from_tabata = rate.get_attenuation(beamlet_geometry.rad, temperatures, densities, 'tabata', True, 2)

    plot_attenuation_profile(shot_number, time, species, energy, dimension, radial_coordinates,
                             [relative_attenuation_rod,
                              relative_attenuation_from_beb_tabata, relative_attenuation_from_nrl_tabata],
                             ['RENATE-OD', 'BEB + Tabata', 'NRL + Tabata'])
    plot_attenuation_profile(shot_number, time, species, energy, dimension, radial_coordinates,
                             [relative_attenuation_rod,
                              relative_attenuation_from_beb_tabata, relative_attenuation_from_nrl_tabata],
                             ['RENATE-OD', 'BEB + Tabata', 'NRL + Tabata'], mode='log')

    plot_attenuation_profile(shot_number, time, species, energy, dimension, radial_coordinates,
                             [relative_attenuation_rod_just_electron,
                              relative_attenuation_from_beb, relative_attenuation_from_nrl],
                             ['RENATE-OD', 'BEB', 'NRL'], scenario='just electron')
    plot_attenuation_profile(shot_number, time, species, energy, dimension, radial_coordinates,
                             [relative_attenuation_rod_just_electron,
                              relative_attenuation_from_beb, relative_attenuation_from_nrl],
                             ['RENATE-OD', 'BEB', 'NRL'], mode='log', scenario='just electron')

    plot_attenuation_profile(shot_number, time, species, energy, dimension, radial_coordinates,
                             [relative_attenuation_rod_just_ion,
                              relative_attenuation_from_tabata],
                             ['RENATE-OD', 'Tabata'], scenario='just ion')
    plot_attenuation_profile(shot_number, time, species, energy, dimension, radial_coordinates,
                             [relative_attenuation_rod_just_ion,
                              relative_attenuation_from_tabata],
                             ['RENATE-OD', 'Tabata'], mode='log', scenario='just ion')


if __name__ == "__main__":
    run_attenuation_comparison(shot_number='17178', time='1097', species='Li', energy='40')
    run_attenuation_comparison(shot_number='17178', time='1097', species='Li', energy='60')
    run_attenuation_comparison(shot_number='17178', time='1097', species='Li', energy='80')
    run_attenuation_comparison(shot_number='17178', time='1097', species='Na', energy='40')
    run_attenuation_comparison(shot_number='17178', time='1097', species='Na', energy='60')
    run_attenuation_comparison(shot_number='17178', time='1097', species='Na', energy='80')
