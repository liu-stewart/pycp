"""test doc."""

from pycp.model.structure import Structure
from pycp.vasp.inputs import Poscar

s = Structure.read("POSCAR")
print(s)