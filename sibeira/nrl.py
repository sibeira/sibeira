import numpy

from sibeira.cross_section import CrossSection


def get_nrl_rate(species, ionisation_level, electron_temperature):
    c = CrossSection(species, ionisation_level)
    t = c.get_t(electron_temperature)
    return 1e-11 * numpy.sqrt(t) / c.B ** 1.5 / (6.0 + t) * numpy.exp(-1.0 / t)
