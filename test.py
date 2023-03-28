"""a test doc."""
from pycp.model.sites import Sites
from pycp.model.lattice import Lattice
from pycp.model.structure import Structure
import numpy as np
import pathlib

s = Structure.from_POSCAR(pathlib.Path("POSCAR"))
print(s)