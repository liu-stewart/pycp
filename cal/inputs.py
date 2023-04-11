"""This module contains the input classes for VASP."""
from __future__ import annotations
import pathlib
from pycp.pycp_typing import Coords, NDArray
from pycp.model.structure import Structure
import numpy as np
from pycp.method.params import setting
import re
import subprocess
import numpy as np


class Poscar(Structure):
    """This class used to represent a POSCAR file."""

    def __str__(self) -> str:
        """Return the string representation of the POSCAR."""
        comment = self.comment
        matrix = self.matrix
        elements = self.elements
        independent_elements = ['start']
        elements_count = []
        for element in elements:
            if element != independent_elements[-1]:
                independent_elements.append(element)
                elements_count.append(1)
            else:
                elements_count[-1] += 1
        independent_elements = independent_elements[1:]
        selective_dynamics = [["T" if i else "F" for i in j]
                              for j in self.selective_dynamics]
        string = f"{comment}\n" + \
            f"   1\n" + \
            f"    {matrix[0][0]:.16f}" + \
            f" {matrix[0][1]:.16f} {matrix[0][2]:.16f}\n" + \
            f"    {matrix[1][0]:.16f}" + \
            f" {matrix[1][1]:.16f} {matrix[1][2]:.16f}\n" + \
            f"    {matrix[2][0]:.16f}" + \
            f" {matrix[2][1]:.16f} {matrix[2][2]:.16f}\n" + \
            f"   {'    '.join(independent_elements)}\n" + \
            f"   {'    '.join([str(i) for i in elements_count])}\n" + \
            f"Selective Dynamics\n" + \
            f"Direct\n"
        for i in range(len(self)):
            string += f"  {self.fractional_coords[i][0]:.16f}" + \
                f" {self.fractional_coords[i][1]:.16f}" + \
                f" {self.fractional_coords[i][2]:.16f}"
            string += f"  {selective_dynamics[i][0]}" + \
                f" {selective_dynamics[i][1]}" + \
                f" {selective_dynamics[i][2]}"
            string += "\n"
        string += "\n"
        for i in range(len(self)):
            string += f"  {self.velocities[i][0]:.16f}" + \
                f" {self.velocities[i][1]:.16f}" + \
                f" {self.velocities[i][2]:.16f}\n"
        return string

    def write(self, to_file: str | pathlib.Path = "POSCAR") -> None:
        """Write the POSCAR to a file."""
        if isinstance(to_file, str):
            to_file = pathlib.Path(to_file)
        if not to_file.parent.exists():
            to_file.parent.mkdir(parents=True)
        with open(to_file, "w+") as file:
            file.write(str(self))


class Incar(dict):
    """This class used to represent an INCAR file."""

    def __init__(self, parameters: dict[str, str | int | float | bool]):
        """Initialize the INCAR."""
        super().__init__(parameters)

    def __str__(self) -> str:
        """Return the string representation of the INCAR."""
        string = ""
        for key, value in self.items():
            if isinstance(value, bool):
                value = ".TRUE." if value else ".FALSE."
            string += f"{key} = {value} \n"
        return string

    @classmethod
    def read(cls, from_file: str | pathlib.Path) -> Incar:
        """Read the INCAR from a file."""
        if isinstance(from_file, str):
            from_file = pathlib.Path(from_file)
        with open(from_file, "r") as file:
            lines = file.readlines()
        parameters = {}
        for line in lines:
            if "=" in line:
                key, value = line.split("=")
                key = key.strip()
                value = value.strip()
                if value == ".TRUE.":
                    value = True
                elif value == ".FALSE.":
                    value = False
                parameters[key] = value
        return cls(parameters)

    def write(self, to_file: str | pathlib.Path = "INCAR") -> None:
        """Write the INCAR to a file."""
        if isinstance(to_file, str):
            to_file = pathlib.Path(to_file)
        if not to_file.parent.exists():
            to_file.parent.mkdir(parents=True)
        with open(to_file, "w+") as file:
            file.write(str(self))


class Kpoints():
    """This class use to represent a Kpoints file."""

    def __init__(self,
                 k: list | list[list],
                 mode: str = "Gamma",
                 number: int = 0,
                 line_Cartesian: bool = False,
                 label: list[str] = [],
                 comment: str = "KPOINTS"
                 ) -> None:
        """Initialize the Kpoints."""
        if np.array(k).ndim == 1:
            k = [k]
        self._k = k
        self.mode = mode
        self.number = number
        self._line_Cartesian = line_Cartesian
        if self.mode.upper().startswith("L") and len(label) != len(k):
            raise ValueError("The length of label must be the same as k.")
        self.label = label
        self.comment = comment

    @property
    def k(self) -> list | list[list]:
        """Return the kpoints."""
        return self._k

    @k.setter
    def k(self, k: list | list[list]) -> None:
        """Set the kpoints."""
        if np.array(k).ndim == 1:
            k = [k]
        self._k = k

    @property
    def line_Cartesian(self) -> bool:
        """Return the line_Cartesian."""
        return self._line_Cartesian

    @line_Cartesian.setter
    def line_Cartesian(self, line_Cartesian: bool) -> None:
        """Set the line_Cartesian."""
        if not isinstance(line_Cartesian, bool):
            raise TypeError("line_Cartesian must be a bool.")
        self._line_Cartesian = line_Cartesian

    def __str__(self) -> str:
        """Return the string representation of the Kpoints."""
        string = f"{self.comment}\n"
        string += f" {self.number}\n"
        if self.mode.upper().startswith("L"):
            string += "Line Mode\n"
            if self.line_Cartesian:
                string += "Cartesian\n"
            else:
                string += "Reciprocal\n"
        elif self.mode.upper().startswith("G"):
            string += "Gamma\n"
        else:
            string += "Monkhorst\n"
        count = 0
        for ak in self.k:
            if self.mode.upper().startswith("L"):
                string += " ".join([f"{i:.5f}" for i in ak])
                if len(self.label) != len(self.k):
                    raise ValueError("The length of "
                                     "label must be the same as k.")
                string += f"  {self.label[self.k.index(ak)]}"
                string += "\n"
                if count < 1:
                    count += 1
                else:
                    count = 0
                    string += "\n"
            else:
                string += " ".join([f"{int(i)}" for i in ak])
                string += "\n"
        string = string.strip()
        return string

    @classmethod
    def read(cls, from_file: str | pathlib.Path) -> Kpoints:
        """Read the Kpoints from a file."""
        if isinstance(from_file, str):
            from_file = pathlib.Path(from_file)
        with open(from_file, "r") as file:
            lines = file.readlines()
        comment = lines[0].strip()
        number = int(lines[1].strip())
        mode = lines[2].strip()
        line_Cartesian = False
        label = []
        if mode.upper().startswith("L"):
            line_Cartesian = lines[3].strip().upper().startswith("C")
            k = []
            label = []
            for line in lines[4:]:
                if line.strip():
                    k.append([float(i) for i in line.split()[:-1]])
                    label.append(line.split()[-1])
        else:
            k = [[float(i) for i in line.split()] for line in lines[3:]]
        return cls(k, mode, number, line_Cartesian, label, comment)

    def write(self, to_file: str | pathlib.Path = "KPOINTS") -> None:
        """Write the Kpoints to a file."""
        if isinstance(to_file, str):
            to_file = pathlib.Path(to_file)
        if not to_file.parent.exists():
            to_file.parent.mkdir(parents=True)
        with open(to_file, "w+") as file:
            file.write(str(self))

    @classmethod
    def from_structure(cls,
                       structure: Structure,
                       ka=40,
                       slab: int = 3) -> Kpoints:
        """Create the Kpoints from a Structure."""
        k = np.ceil(40/structure.length)  # type: ignore
        if slab != 0:
            k[slab-1] = 1
        return cls(k.tolist())


class Potcar():
    """This class used to represent a POTCAR file."""

    def __init__(self,
                 elements: list[str] = [],
                 func: str = "PAW_PBE") -> None:
        """Initialize the Potcar."""
        self.elements = elements
        self.func = func

    def write(self, to_file: pathlib.Path | str = "POTCAR"):
        """Write the POTCAR to a file."""
        if isinstance(to_file, str):
            to_file = pathlib.Path(to_file)
        pot_path = str(setting["POTCAR_PATH"][self.func])
        if pot_path.startswith("~"):
            pot_path = pathlib.Path.home()/pot_path.lstrip("~/")
        pot = ""
        for element in self.elements:
            with open(pot_path/element/"POTCAR") as file:  # type: ignore
                pot += file.read()
        with open(to_file, "w") as file:
            file.write(pot)

    @classmethod
    def read(cls, from_file: pathlib.Path | str):
        """Read the elements from POTCAR file."""
        if isinstance(from_file, str):
            from_file = pathlib.Path(from_file)
        title = subprocess.run(f"grep TIT {from_file}",
                               stdout=subprocess.PIPE,
                               shell=True).stdout
        title = str(title, encoding="utf-8").splitlines()
        elements = [t.split()[-2] for t in title]
        return cls(elements)

    @classmethod
    def from_structure(cls, structure: Structure):
        """Create Potcar from structure."""
        elements = structure.elements
        independent_elements = ['start']
        for element in elements:
            if element != independent_elements[-1]:
                independent_elements.append(element)
        independent_elements = independent_elements[1:]

        return cls(independent_elements)