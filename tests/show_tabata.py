import numpy
import matplotlib.pyplot
from tabata_ctf.cross_section import CrossSection


def plot_tabata_profile(energy, cross_section, species, degree, factor=1):
    if degree == 'single':
        ylabel_text = '$\\sigma_{\\mathrm{' + species + '+D}^+\\rightarrow ' \
                                                        '\\mathrm{' + species + '}^++\\mathrm{D}}$ [cm$^2$]'
    elif degree == 'double':
        ylabel_text = '$\\sigma_{\\mathrm{2\,' + species + '+D}^+\\rightarrow ' \
                                                        '\\mathrm{2\,' + species + '}^++\\mathrm{D}^-}$ [cm$^2$]'
    else:
        ylabel_text = '$\\sigma$ [cm$^2$]'

    fig, ax = matplotlib.pyplot.subplots()
    matplotlib.pyplot.subplots_adjust(left=0.2, right=0.96)
    fig.set_size_inches(4, 4.78)
    matplotlib.pyplot.loglog(energy, cross_section)
    matplotlib.pyplot.ylim(1e-19*factor, 1e-13*factor)
    matplotlib.pyplot.xlim(10, 1e6)
    matplotlib.pyplot.xlabel('E [eV]')
    matplotlib.pyplot.ylabel(ylabel_text)
    ax.xaxis.set_label_coords(0.5, -0.09)
    ax_img = fig.add_subplot(111, label='ax_img')
    img = matplotlib.pyplot.imread('tabata1988_' + species + '_' + degree + '.jpg')
    crop_img = img[62:3425, 467:3252]
    ax_img.set_zorder(999)
    ax.set_zorder(1000)
    ax.patch.set_alpha(0.1)
    ax_img.imshow(crop_img, aspect='1', alpha=1)
    ax_img.set_yticks([])
    ax_img.set_xticks([])
    matplotlib.pyplot.show()


def show_single_lithium_cross_exchange_ionisation():
    species = 'Li'
    degree = 'single'
    energy = numpy.logspace(1.3, 6, 50)
    c = CrossSection(energy, species, degree)
    cross_section = c.calculate() * 1e4
    plot_tabata_profile(energy, cross_section, species, degree)


def show_double_magnesium_cross_exchange_ionisation():
    species = 'Mg'
    degree = 'double'
    energy = numpy.logspace(2, 6, 50)
    c = CrossSection(energy, species, degree)
    cross_section = c.calculate() * 1e4
    plot_tabata_profile(energy, cross_section, species, degree, factor=0.1)


if __name__ == '__main__':
    show_single_lithium_cross_exchange_ionisation()
    show_double_magnesium_cross_exchange_ionisation()
