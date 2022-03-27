import unittest

import numpy

from renate_od import plot_attenuation_profile


class TestRenateODPlot(unittest.TestCase):
    def test_plot_attenuation_profile(self):
        r = numpy.linspace(0.6, 0.74, 11)
        p1 = numpy.random.random_sample(r.shape)
        p2 = numpy.random.random_sample(r.shape)
        p3 = numpy.random.random_sample(r.shape)
        plot_attenuation_profile(12345, 678, 'test', 123, 3, r, [p1, p2, p3], ['p1', 'p2', 'p3'])


if __name__ == '__main__':
    unittest.main()
