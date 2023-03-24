"""a test doc."""
from pycp.model.lattice import Lattice

la = Lattice.from_lengths_and_angles([12.1988, 12.9868, 7.0757],
                                     [116.997, 106.317, 97.917])
print(la.matrix)
