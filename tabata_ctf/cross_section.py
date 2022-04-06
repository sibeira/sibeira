import os
import numpy
import scipy.interpolate
import pandas


# T Tabata et al., Nucl. Inst. Meth. Phys. Res. B, 31 (3), 1988
class CrossSection:
    def __init__(self, species, degree='single'):
        self.species = species
        self.degree = degree
        self.tabata_data = self.get_tabata_data()

    def get_tabata_data(self):
        tabata_dataframe = pandas.read_csv(os.path.dirname(__file__) + '/' + self.degree + '.dat', sep='\t')
        tabata_dataframe.replace('Infinity', numpy.inf, inplace=True)
        tabata_data = tabata_dataframe[tabata_dataframe['target'] == self.species]
        if tabata_data.empty:
            raise ValueError('Invalid species for Tabata database: ' + self.species)
        return tabata_data

    def get_f(self, E1):
        ER = 25.00
        if self.tabata_data['a3'].to_numpy() == numpy.inf:
            return self.tabata_data['a1'].to_numpy() * (E1 / ER) ** self.tabata_data['a2'].to_numpy() / \
                   (1.0 +
                    (E1 / self.tabata_data['a5'].to_numpy()) **
                    (self.tabata_data['a2'].to_numpy() + self.tabata_data['a6'].to_numpy()))
        return self.tabata_data['a1'].to_numpy() * (E1 / ER) ** self.tabata_data['a2'].to_numpy() / \
               (1.0 +
                (E1 / self.tabata_data['a3'].to_numpy()) **
                (self.tabata_data['a2'].to_numpy() + self.tabata_data['a4'].to_numpy()) +
                (E1 / self.tabata_data['a5'].to_numpy()) **
                (self.tabata_data['a2'].to_numpy() + self.tabata_data['a6'].to_numpy()))

    @staticmethod
    def replace_nan_to_zero(a):
        nan_indices = numpy.isnan(a.astype(numpy.float))
        a[nan_indices] = 0.0
        return a

    def calculate(self, energy):
        try:
            E1 = energy / 1000.0 - self.tabata_data['Et'].to_numpy()
            sigma0 = 1e-20
            cross_section = sigma0 * (self.get_f(E1) +
                                      self.tabata_data['a7'].to_numpy() *
                                      self.get_f(E1 / self.tabata_data['a8'].to_numpy()))

            return self.replace_nan_to_zero(cross_section)
        except KeyError:
            raise KeyError('Broken database, missing argument')

    def get_polynomial(self):
        energy = numpy.logspace(0.75, 5, 50)
        cross_section = self.calculate(energy)
        return scipy.interpolate.interp1d(energy, cross_section, fill_value='extrapolate')
