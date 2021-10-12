import scipy.constants

RYDBERG = scipy.constants.Rydberg * scipy.constants.Planck * scipy.constants.c / scipy.constants.elementary_charge


class OrbitalConstants:
    def __init__(self, species):
        self.orbital_constants = self.__get_orbital_constants(species)
        self.fill()

    def fill(self):
        if 'N' not in self.orbital_constants:
            self.orbital_constants['N'] = 1
        if 'Q' not in self.orbital_constants:
            self.orbital_constants['Q'] = 1
        if 'n' not in self.orbital_constants:
            self.fill_n()
        if 'U' not in self.orbital_constants:
            self.fill_U()

    def fill_n(self):
        if 'pqn' not in self.orbital_constants:
            self.orbital_constants['pqn'] = 1
        pqn = self.orbital_constants['pqn']
        if pqn == 2:
            pqn = 1
        self.orbital_constants['n'] = pqn

    def fill_U(self):
        self.orbital_constants['U'] = self.orbital_constants['B']

    @staticmethod
    def __get_orbital_constants(species):
        return orbital_constants[species]

    def get(self, field):
        return self.orbital_constants[field]


orbital_constants = {
    'H':  dict(pqn=1, B=RYDBERG, Q=0.5668),
    'He': dict(pqn=1, B=24.587389011, U=3*RYDBERG, N=2, Q=0.8841),
    'Li': dict(pqn=2, B=5.391714996),
    'Na': dict(pqn=3, B=5.13907696),
    'K':  dict(pqn=4, B=4.34066373),
    'Rb': dict(pqn=5, B=4.1771281),
    'Cs': dict(pqn=6, B=3.89390572743)
}


class IonOrbitalConstants(OrbitalConstants):
    # noinspection PyMissingConstructor
    def __init__(self, species):
        self.orbital_constants = self.__get_orbital_constants(species)
        self.fill()

    @staticmethod
    def __get_orbital_constants(species):
        return ion_orbital_constants[species]

    def fill_U(self):
        self.orbital_constants['U'] = (self.orbital_constants['Z']-5.0/16.0)**2 * 2.0 * RYDBERG


ion_orbital_constants = {
    'He': dict(pqn=1, Z=2,  B=4*RYDBERG, U=4*RYDBERG, Q=0.8841),
    'Li': dict(pqn=2, Z=3,  B=75.6400970, U=14.56*RYDBERG, N=2),
    'Na': dict(pqn=3, Z=11, B=47.28636, N=6),
    'K':  dict(pqn=4, Z=19, B=31.62500, N=6),
    'Rb': dict(pqn=5, Z=37, B=27.28954, N=6),
    'Cs': dict(pqn=6, Z=55, B=23.15745, N=6)
}
