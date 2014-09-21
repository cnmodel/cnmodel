from .population import Population

class Bushy(Population):
    def __init__(self, species='mouse'):
        self._species = species

    