"""test doc."""

from pycp.model.structure import Structure
from pycp.cal.inputs import Poscar, Incar, Kpoints, Potcar
from pycp.cal.script import Script
from pycp.cal.run import RunVasp
from pycp.model.lattice import Lattice
from pycp.model.sites import Sites
import copy
import pathlib
import sys

incar = Incar.read("INCAR")
poscar = Poscar.read("POSCAR.vasp")
vasp = RunVasp(
    incar=incar,
    poscar=poscar,
    work_dir="s"
)
vasp.run()