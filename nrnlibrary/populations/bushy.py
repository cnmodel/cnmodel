from .population import Population

class Bushy(Population):
    def __init__(self, species='mouse'):
        # Completely fabricated cell distribution: uniform from 4kHz to 90kHz.
        # Note that `cf` is the mean value used when selecting SGCs to connect;
        # it is NOT the measured CF of the cell (although it should be close).
        size = 5000
        super(Bushy, self).__init__(species, size, fields=[('cf', float)])
        self._cells['cf'] = 4000 * 2**np.linspace(0, 4.5, size)
    
    