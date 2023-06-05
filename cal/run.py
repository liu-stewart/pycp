"""This module is used to run the script on the inputs."""
from __future__ import annotations
import pathlib
from pycp.cal.inputs import Poscar, Incar, Kpoints, Potcar
from pycp.cal.script import Script
import os
import subprocess
import re
import sys


class RunVasp():
    """This class is used to run the VASP."""

    def __init__(self,
                 poscar: Poscar,
                 incar: Incar,
                 kpoints: Kpoints | None = None,
                 potcar: Potcar | None = None,
                 script: Script | None = None,
                 work_dir: pathlib.Path | str = "defaultRun",
                 continue_list=[],
                 slab: int = 3
                 ):
        """Init this class."""
        self.poscar = poscar
        self.incar = incar
        if kpoints is None:
            kpoints = Kpoints.from_structure(poscar, slab=slab)
        self.kpoints = kpoints
        if potcar is None:
            potcar = Potcar.from_structure(poscar)
        self.potcar = potcar
        if script is None:
            script = Script()
        self.script = script
        if work_dir == "defaultRun":
            work_dir = pathlib.Path(sys.argv[0][:-3])
        self.work_dir = pathlib.Path(work_dir)
        self.continue_list = continue_list
        self.init_dir = pathlib.Path(".").absolute()

    def write(self) -> None:
        """Write the inputs and script to the work_dir."""
        self.work_dir.mkdir(exist_ok=True, parents=True)
        os.chdir(self.work_dir)
        self.poscar.write()
        self.incar.write()
        if "KSPACING" not in self.incar.keys():
            self.kpoints.write()
        self.potcar.write()
        self.script.write()
        os.chdir(self.init_dir)

    def run(self) -> int:
        """Run the script."""
        self.write()
        os.chdir(self.work_dir)
        result = subprocess.getoutput("bsub < script.lsf")
        print(result)
        pid = re.match('.*<([0-9]+)>', result)[1]  # type: ignore
        os.chdir(self.init_dir)
        return pid  # type: ignore
