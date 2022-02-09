import numpy
from bebim.beam import Beam
from bebim.cross_section import CrossSection
from integrator import ElectronRateIntegrator1D, IonRateIntegrator2D, IonRateIntegrator3D
from tabata_ctf.cross_section import CrossSection as CXCrossSection


class Rate(Beam):
    def __init__(self, species, beam_energy, ionisation_level=0):
        super().__init__(species, beam_energy, ionisation_level)
        self.beb = CrossSection(1000, self.species, self.ionisation_level)
        self.beb.set_polynomial()
        self.tabata = CXCrossSection(1000, self.species)
        self.tabata.set_polynomial()
        self.tabata_double = CXCrossSection(1000, self.species, degree='double')
        self.tabata_double.set_polynomial()

    def set_profiles(self, electron_temperature=numpy.nan, electron_density=numpy.nan):
        if ~numpy.isnan(electron_temperature):
            self.electron_temperature = electron_temperature
        if ~numpy.isnan(electron_density):
            self.electron_density = electron_density

    def get_attenuation_nrl(self, is_with_tabata=False, tabata_integration_dimension=2):
        if self.electron_density == 0:
            return 0.0
        c = CrossSection(self.electron_temperature, self.species, self.ionisation_level)
        t = c.get_t()
        r = 1e-11 * numpy.sqrt(t) / c.B ** 1.5 / (6.0 + t) * numpy.exp(-1.0 / t)
        if is_with_tabata:
            if tabata_integration_dimension == 2:
                r += IonRateIntegrator2D(self.electron_temperature, self.speed, self.tabata.f, self.tabata_double.f).get_coefficient()
            elif tabata_integration_dimension == 3:
                r += IonRateIntegrator3D(self.electron_temperature, self.speed, self.tabata.f, self.tabata_double.f).get_coefficient()
            else:
                raise ValueError
        return r * self.electron_density

    def get_attenuation_beb(self, is_with_tabata=False, tabata_integration_dimension=2):
        if self.electron_density == 0:
            return 0.0
        r = ElectronRateIntegrator1D(self.electron_temperature, self.beb.f).get_coefficient()
        if is_with_tabata:
            if tabata_integration_dimension == 2:
                r += IonRateIntegrator2D(self.electron_temperature, self.speed, self.tabata.f, self.tabata_double.f).get_coefficient()
            elif tabata_integration_dimension == 3:
                r += IonRateIntegrator3D(self.electron_temperature, self.speed, self.tabata.f, self.tabata_double.f).get_coefficient()
            else:
                raise ValueError
        return r * self.electron_density
