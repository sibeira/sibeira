import scipy.constants

RYDBERG = scipy.constants.Rydberg * scipy.constants.Planck * scipy.constants.c / scipy.constants.elementary_charge


class OrbitalConstants:
    def __init__(self, species):
        self.orbital_constants = orbital_constants[species]
        if 'N' not in self.orbital_constants:
            self.orbital_constants['N'] = 1
        if 'Q' not in self.orbital_constants:
            self.orbital_constants['Q'] = 1
        if 'n' not in self.orbital_constants:
            self.orbital_constants['n'] = 1
        if 'B' not in self.orbital_constants:
            self.orbital_constants['B'] = RYDBERG
        if 'U' not in self.orbital_constants:
            self.orbital_constants['U'] = RYDBERG

    def get(self, field):
        return self.orbital_constants[field]


orbital_constants = {
    'H': dict(B=RYDBERG, U=RYDBERG, N=1, Q=0.5668),
    'He': dict(B=24.587, U=3*RYDBERG, N=2, Q=0.8841)
}
