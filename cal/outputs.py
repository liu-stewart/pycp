"""This module contains the output classes for VASP."""
from __future__ import annotations
import pathlib
from pycp.pycp_typing import Coords, NDArray
from pycp.model.structure import Structure
import numpy as np
from pycp.method.params import setting
import re
import subprocess
import numpy as np


class Oszicar():
    """This class used to represent an OSZICAR file."""

    @classmethod
    def getEnergy(cls, from_file: str | pathlib.Path) -> float:
        """Read the energy from an OSZICAR file."""
        if isinstance(from_file, str):
            from_file = pathlib.Path(from_file)
        if from_file.is_dir():
            from_file = from_file / "OSZICAR"
        with open(from_file, "r") as file:
            lines = file.readlines()
        for line in lines[::-1]:
            if "F=" in line:
                energy = float(line.split()[4])
                return energy
        raise ValueError("No energy found in OSZICAR file.")
