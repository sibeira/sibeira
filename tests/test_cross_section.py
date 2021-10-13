import numpy
import unittest
import matplotlib.pyplot
from cross_section import CrossSection


def plot_profiles(energy, cross_section, reference_energy, reference_cross_section,
                  reference_relative_error=0):
    matplotlib.pyplot.semilogx(energy, cross_section, '.')
    matplotlib.pyplot.errorbar(reference_energy, reference_cross_section,
                               yerr=reference_cross_section*reference_relative_error)
    matplotlib.pyplot.show()


def plot_uddin_profile(energy, cross_section):
    fig, ax = matplotlib.pyplot.subplots()
    fig.set_size_inches(4, 4.3)
    matplotlib.pyplot.semilogx(energy, cross_section)
    matplotlib.pyplot.ylim(0, 5)
    matplotlib.pyplot.xlim(1, 1e4)
    ax_img = fig.add_subplot(111, label="ax_img")
    img = matplotlib.pyplot.imread("uddin2015.jpg")
    crop_img = img[10:930, 100:960]
    ax_img.set_zorder(999)
    ax.set_zorder(1000)
    ax.patch.set_alpha(0.1)
    ax_img.imshow(crop_img, aspect='1', alpha=1)
    ax_img.set_yticks([])
    ax_img.set_xticks([])
    matplotlib.pyplot.show()


def plot_wareing_profile(normalised_energy, cross_section):
    fig, ax = matplotlib.pyplot.subplots()
    fig.set_size_inches(5, 2.25)
    matplotlib.pyplot.plot(normalised_energy, cross_section)
    matplotlib.pyplot.ylim(0, 5)
    matplotlib.pyplot.xlim(0, 14)
    matplotlib.pyplot.xlabel('t')
    matplotlib.pyplot.ylabel(r'$\sigma_{Li^+\rightarrow Li^{2+}}$ [$10^{-18}$ cm$^2$]')
    ax.xaxis.set_label_coords(1.08, -0.05)
    ax_img = fig.add_subplot(111, label="ax_img")
    img = matplotlib.pyplot.imread("wareing1967.jpg")
    crop_img = img[32:465, 130:1100]
    ax_img.set_zorder(999)
    ax.set_zorder(1000)
    ax.patch.set_alpha(0.1)
    ax_img.imshow(crop_img, aspect='1', alpha=1)
    ax_img.set_yticks([])
    ax_img.set_xticks([])
    matplotlib.pyplot.show()


def read_reference(species):
    reference_data = numpy.loadtxt('reference_cross_section_' + species + '.dat')
    energy = reference_data[:, 0]
    cross_section = reference_data[:, 1]
    try:
        relative_error = reference_data[:, 2]/100.0
        return energy, cross_section, relative_error
    except IndexError:
        return energy, cross_section


class TestCrossSection(unittest.TestCase):

    def test_hydrogen_with_Q(self):
        reference_energy, reference_cross_section = read_reference('H_with_Q')

        energy = reference_energy
        c = CrossSection(energy, 'H', 0)
        cross_section = c.calculate() / 1e-20

        numpy.testing.assert_array_almost_equal(cross_section, reference_cross_section, decimal=4,
                                                err_msg='Cross section test for hydrogen')

    def test_hydrogen(self):
        reference_energy, reference_cross_section = read_reference('H')

        energy = reference_energy
        c = CrossSection(energy, 'H', 0, with_Q=False)
        cross_section = c.calculate() / 1e-20

        numpy.testing.assert_array_almost_equal(cross_section, reference_cross_section, decimal=4,
                                                err_msg='Cross section test for hydrogen (Q=1)')

    def test_helium_with_Q(self):
        reference_energy, reference_cross_section = read_reference('He_with_Q')

        energy = reference_energy
        c = CrossSection(energy, 'He', 0)
        cross_section = c.calculate() / 1e-20

        numpy.testing.assert_array_almost_equal(cross_section, reference_cross_section, decimal=3,
                                                err_msg='Cross section test for helium')

    def test_helium(self):
        reference_energy, reference_cross_section = read_reference('He')

        an_energy = reference_energy
        c = CrossSection(an_energy, 'He', 0, with_Q=False)
        cross_section = c.calculate() / 1e-20

        numpy.testing.assert_array_almost_equal(cross_section, reference_cross_section, decimal=3,
                                                err_msg='Cross section test for helium (Q=1)')

    def test_plot_lithium(self):
        energy = numpy.logspace(0.75, 5, 50)
        c = CrossSection(energy, 'Li', 0, with_Q=False)
        cross_section = c.calculate() / 1e-20
        plot_uddin_profile(energy, cross_section)

    def test_lithium_secondary_ionisation(self):
        energy = numpy.logspace(3.5, 4.5, 50)
        c = CrossSection(energy, 'Li', 1, with_Q=False)
        cross_section = c.calculate() / 1e-22
        reference_energy, reference_cross_section, reference_relative_error = read_reference('Li2')
        plot_profiles(energy, cross_section, reference_energy, reference_cross_section, reference_relative_error)

    def test_plot_lithium_secondary_ionisation(self):
        energy = numpy.linspace(75, 1000, 50)
        c = CrossSection(energy, 'Li', 1, with_Q=False)
        cross_section = c.calculate() / 1e-22
        normalised_energy = c.get_t()
        plot_wareing_profile(normalised_energy, cross_section)


if __name__ == '__main__':
    unittest.main()
