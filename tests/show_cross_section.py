import numpy
import matplotlib.pyplot
from bebim.cross_section import CrossSection
from tests.read_reference import read_reference


def plot_profiles(reaction, y_unit, energy, cross_section, reference_energy, reference_cross_section,
                  reference_relative_error=0):
    matplotlib.pyplot.semilogx(energy, cross_section, '.', label='simulation')
    matplotlib.pyplot.errorbar(reference_energy, reference_cross_section,
                               yerr=reference_cross_section * reference_relative_error, label='reference')
    y_exponent = str(int(numpy.log10(y_unit) + 4))
    matplotlib.pyplot.xlabel('T [eV]')
    matplotlib.pyplot.ylabel(r'$\sigma_{' + reaction + '}$ [$10^{' + y_exponent + '}$ cm$^2$]')
    matplotlib.pyplot.legend()
    matplotlib.pyplot.xlim(3e3,3e4)
    matplotlib.pyplot.ylim(0.0,4)
    matplotlib.pyplot.show()


def plot_uddin_profile(energy, cross_section):
    energy_beli, cross_section_beli = read_reference('Li_Beli')
    energy_lotz, cross_section_lotz = read_reference('Li_Lotz')
    fig, ax = matplotlib.pyplot.subplots()
    fig.set_size_inches(4, 4.3)
    matplotlib.pyplot.semilogx(energy, cross_section, label='simulation')
    matplotlib.pyplot.semilogx(energy_beli, cross_section_beli*1e16, label='Beli')
    matplotlib.pyplot.semilogx(energy_lotz, cross_section_lotz*1e16, label='Lotz')
    matplotlib.pyplot.ylim(0, 5)
    matplotlib.pyplot.xlim(1, 1e4)
    matplotlib.pyplot.xlabel('T [eV]')
    matplotlib.pyplot.ylabel(r'$\sigma_{Li\rightarrow Li^+}$ [$10^{-16}$ cm$^2$]')
    matplotlib.pyplot.legend(bbox_to_anchor=(1.15, 0.35))
    ax.xaxis.set_label_coords(0.5, -0.09)
    ax_img = fig.add_subplot(111, label='ax_img')
    img = matplotlib.pyplot.imread('tests/figures/uddin2015.jpg')
    crop_img = img[10:930, 100:960]
    ax_img.set_zorder(999)
    ax.set_zorder(1000)
    ax.patch.set_alpha(0.1)
    ax_img.imshow(crop_img, aspect='1', alpha=1)
    ax_img.set_yticks([])
    ax_img.set_xticks([])
    matplotlib.pyplot.show()


def plot_wareing_profile(normalised_energy, cross_section):
    energy_beli, cross_section_beli = read_reference('Li2_Beli')
    fig, ax = matplotlib.pyplot.subplots()
    fig.set_size_inches(5, 2.25)
    matplotlib.pyplot.plot(normalised_energy, cross_section, label='simulation')
    matplotlib.pyplot.plot(energy_beli/75.6400970, cross_section_beli*1e18, label='Beli')
    matplotlib.pyplot.ylim(0, 5)
    matplotlib.pyplot.xlim(0, 14)
    matplotlib.pyplot.xlabel('T/U')
    matplotlib.pyplot.ylabel(r'$\sigma_{Li^+\rightarrow Li^{2+}}$ [$10^{-18}$ cm$^2$]')
    matplotlib.pyplot.legend(loc='lower right')
    ax.xaxis.set_label_coords(1.08, -0.05)
    ax_img = fig.add_subplot(111, label='ax_img')
    img = matplotlib.pyplot.imread('tests/figures/wareing1967.jpg')
    crop_img = img[32:465, 130:1100]
    ax_img.set_zorder(999)
    ax.set_zorder(1000)
    ax.patch.set_alpha(0.1)
    ax_img.imshow(crop_img, aspect='1', alpha=1)
    ax_img.set_yticks([])
    ax_img.set_xticks([])
    matplotlib.pyplot.show()


def show_hydrogen_primary_ionisation():
    reference_energy, reference_cross_section = read_reference('H')
    energy = reference_energy
    c = CrossSection('H', 0, with_Q=False)
    y_unit = 1e-20
    cross_section = c.calculate(energy) / y_unit
    plot_profiles('H', y_unit, energy, cross_section, reference_energy, reference_cross_section)


# MA Uddin et al., Phys. Rev. A 72, 032715, 2015
def show_lithium_primary_ionisation():
    energy = numpy.logspace(0.75, 5, 50)
    c = CrossSection('Li', 0, with_Q=False)
    cross_section = c.calculate(energy) / 1e-20
    plot_uddin_profile(energy, cross_section)


# JB Wareing et al., Proc. Phys. Soc. 91 887, 1967
def show_lithium_secondary_ionisation():
    energy = numpy.linspace(75, 1000, 50)
    c = CrossSection('Li', 1, with_Q=False)
    cross_section = c.calculate(energy) / 1e-22
    normalised_energy = c.get_t()
    plot_wareing_profile(normalised_energy, cross_section)


# B Peart et al., Journ. Phys. B 2(12), 1969
def show_lithium_secondary_ionisation_high_energies():
    energy = numpy.logspace(3.5, 4.5, 50)
    c = CrossSection('Li', 1, with_Q=False)
    y_unit = 1e-22
    cross_section = c.calculate(energy) / y_unit
    reference_energy, reference_cross_section, reference_relative_error = read_reference('Li2')
    plot_profiles('Li^+->Li^{2+}', y_unit, energy, cross_section, reference_energy, reference_cross_section, reference_relative_error)


    energy_beli, cross_section_beli = read_reference('Li_Beli')
    plot_profiles('Li^+->Li^{2+} BELI', y_unit, energy_beli, cross_section_beli*1e18, reference_energy, reference_cross_section, reference_relative_error)


if __name__ == '__main__':
    show_hydrogen_primary_ionisation()
    show_lithium_primary_ionisation()
    show_lithium_secondary_ionisation()
    show_lithium_secondary_ionisation_high_energies()
