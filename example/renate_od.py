import sys
import os
import matplotlib
import numpy

from rate import Rate

sys.path.append(os.environ['RENATE_OD'])

from manager import RenateODManager
from beamlet import set_beamlet, BeamletGeometry
#from efit import EFITManager
from profiles import Profiles
from utils import *


def plot_attenuation_profile(shot_number, time, species, energy, radial_coordinate,
                             relative_attenuation, relative_attenuation2, relative_attenuation3, relative_attenuation4):
    fig, ax = matplotlib.pyplot.subplots()
    fig.set_size_inches(5, 2)
    ax.plot(radial_coordinate, relative_attenuation, '-', linewidth=2, label='RENATE-OD')
    #ax.plot(radial_coordinate, relative_attenuation2, '-', linewidth=1, label='BEB')
    #ax.plot(radial_coordinate, relative_attenuation3, '--', linewidth=1.5, label='NRL')
    ax.plot(radial_coordinate, relative_attenuation4, 'r--', linewidth=1, label='BEB+Tabata 3D')
    ax.legend()

    matplotlib.pyplot.xlim(0.6, 0.7399)
    matplotlib.pyplot.minorticks_on()
    matplotlib.pyplot.grid(which='major')
    matplotlib.pyplot.xlabel('$R$ [m]', labelpad=-10.5, loc='right')
    matplotlib.pyplot.ylabel('neutral beam attenuation')
    matplotlib.pyplot.title('COMPASS #' + shot_number + ' (' + time + ' ms, '+species + ', ' + energy + ' keV)')

    matplotlib.pyplot.savefig('atomic.png')
    matplotlib.pyplot.savefig('atomic.pdf')
    matplotlib.pyplot.show()

    fig, ax = matplotlib.pyplot.subplots()
    fig.set_size_inches(5, 2)
    ax.semilogy(radial_coordinate, relative_attenuation, '-', linewidth=2, label='RENATE-OD')
    ax.semilogy(radial_coordinate, relative_attenuation2, '-', linewidth=1, label='BEB')
    ax.semilogy(radial_coordinate, relative_attenuation3, '--', linewidth=1.5, label='NRL')
    ax.semilogy(radial_coordinate, relative_attenuation4, 'r--', linewidth=1, label='BEB+Tabata 3D')
    ax.legend()

    matplotlib.pyplot.xlim(0.6, 0.7399)
    matplotlib.pyplot.minorticks_on()
    matplotlib.pyplot.grid(which='major')
    matplotlib.pyplot.xlabel('$R$ [m]', labelpad=-10.5, loc='right')
    matplotlib.pyplot.ylabel('neutral beam attenuation')
    matplotlib.pyplot.title('COMPASS #' + shot_number + ' (' + time + ' ms)')

    matplotlib.pyplot.savefig('atomic_log.png')
    matplotlib.pyplot.savefig('atomic_log.pdf')
    matplotlib.pyplot.show()


def export_beamlet_profile(shot_number='17178', time='1097', species='Li', energy='80'):
    z = 0
    tor = 0
    beamlet_geometry = set_beamlet(z, tor)

    r = RenateODManager(beamlet_geometry, shot_number, time, species, energy)
    radial_coordinates, relative_attenuation_rod = r.get_attenuation_profile()
    relative_attenuation_rod = relative_attenuation_rod.fillna(0)

    p = Profiles(is_export=False)
    if not numpy.array_equal(radial_coordinates.to_numpy, beamlet_geometry.rad):
        assert ValueError
    temperatures = p.get_temperature()
    densities = p.get_density()
    relative_attenuation_1d = numpy.zeros_like(beamlet_geometry.rad)
    relative_attenuation_nrl = numpy.zeros_like(beamlet_geometry.rad)
    relative_attenuation_3d = numpy.zeros_like(beamlet_geometry.rad)

    rate = Rate(species=species, beam_energy=float(energy) * 1000.0)
    v_beam = rate.get_speed()
    timestep = (beamlet_geometry.rad[0] - beamlet_geometry.rad[1]) / v_beam

    for i in range(beamlet_geometry.rad.size):
        print(i)
        rate.set_profiles(electron_temperature=temperatures[i], electron_density=densities[i])
        if i == 0:
            relative_attenuation_1d[i] = 1.0 - rate.get_attenuation_beb() * timestep
            relative_attenuation_nrl[i] = 1.0 - rate.get_attenuation_nrl() * timestep
            relative_attenuation_3d[i] = 1.0 - rate.get_attenuation_beb(True, 3) * timestep
        else:
            relative_attenuation_1d[i] = (1.0 - rate.get_attenuation_beb() * timestep) * relative_attenuation_1d[i - 1]
            relative_attenuation_nrl[i] = (1.0 - rate.get_attenuation_nrl() * timestep) * \
                                          relative_attenuation_nrl[i - 1]
            relative_attenuation_3d[i] = (1.0 - rate.get_attenuation_beb(True, 3) * timestep) * \
                                          relative_attenuation_3d[i - 1]
        if relative_attenuation_nrl[i] <= 0:
            relative_attenuation_nrl[i] = 0
        if i==100:
            break

    plot_attenuation_profile(shot_number, time, species, energy, radial_coordinates,
                             relative_attenuation_rod, relative_attenuation_1d, relative_attenuation_nrl, relative_attenuation_3d)


if __name__ == "__main__":
    a_shot_number = '17178'
    a_time = '1097'
    export_beamlet_profile(shot_number=a_shot_number, time=a_time)
