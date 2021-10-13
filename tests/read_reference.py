import numpy


def read_reference(species):
    reference_data = numpy.loadtxt('reference_cross_section_' + species + '.dat')
    energy = reference_data[:, 0]
    cross_section = reference_data[:, 1]
    try:
        relative_error = reference_data[:, 2]/100.0
        return energy, cross_section, relative_error
    except IndexError:
        return energy, cross_section
