"""This module contains the class Structure.

which is used to represent the structure of a solid.
"""

from __future__ import annotations
import numpy as np
from monty.json import MSONable
from pycp.pycp_typing import Coords, NDArray
from pycp.model.sites import Sites
from pycp.model.lattice import Lattice
from pycp.pattern import pattern_element
import pathlib


class Structure(Sites):
    """This class is used to represent the structure of a solid.

    Properties:
        lattice: The lattice of the structure.
        cartesian_coords: The cartesian coordinates of the sites.
        fractional_coords: The fractional coordinates of the sites.

    Methods:
        periodic_boundary_conditions: Apply periodic boundary conditions to the
            fractional coordinates.
        rotate: Rotate the structure.
        translate: Translate the structure.
        axial_symmetry: Apply axial symmetry to the structure.
        read: Read a structure from a file.
    """

    def __init__(self,
                 coordinates: Coords,
                 elements: list[str] | str,
                 lattice: Lattice,
                 *args):
        """Initialize a Structure.

        Args:
            coordinates: The coordinates of the sites.
            elements: The elements of the sites.
            lattice: The lattice of the structure.
        """
        if not isinstance(lattice, Lattice):
            raise TypeError("lattice must be a Lattice object")
        self.__lattice = lattice
        super().__init__(coordinates, elements)
        self.periodic_boundary_conditions()

    @property
    def lattice(self) -> Lattice:
        """Return the lattice of the structure."""
        return self.__lattice

    @lattice.setter
    def lattice(self, lattice: Lattice) -> None:
        """Set the lattice of the structure."""
        if not isinstance(lattice, Lattice):
            raise TypeError("lattice must be a Lattice object")
        self.__lattice = lattice
        self.periodic_boundary_conditions()

    @property
    def cartesian_coords(self) -> NDArray:
        """Return the cartesian coordinates of the sites."""
        return self.coordinates

    @property
    def fractional_coords(self) -> NDArray:
        """Return the fractional coordinates of the sites."""
        return np.linalg.solve(self.lattice.matrix, self.cartesian_coords.T).T

    def periodic_boundary_conditions(self) -> None:
        """Apply periodic boundary conditions to the fractional coordinates."""
        coords = self.fractional_coords
        coords -= np.floor(coords)
        self.coordinates = coords @ self.__lattice.matrix

    def rotate(self, angle: float,
               axis: Coords = [0, 0, 1],
               anchor: Coords = [0, 0, 0]) -> None:
        """Rotate the structure."""
        super().rotate(angle, axis, anchor)
        self.periodic_boundary_conditions()

    def translate(self, vector: Coords) -> None:
        """Translate the structure."""
        super().translate(vector)
        self.periodic_boundary_conditions()

    def axial_symmetry(self,
                       axis: Coords = ...,
                       anchor: Coords = [0, 0, 0]) -> None:
        """Apply axial symmetry to the structure."""
        super().axial_symmetry(axis, anchor)
        self.periodic_boundary_conditions()

    def __repr__(self) -> str:
        """Return the representation of the structure."""
        return super().__repr__() + f"\n\n{self.lattice}"

    def __str__(self) -> str:
        """Return the string representation of the structure."""
        return super().__str__() + f"\n\n{self.lattice}"

    @classmethod
    def read(cls, file: pathlib.Path | str, fmt: str = ""):
        """Create a Structure object from a file.

        Args:
            file: The path of the file.
            fmt: The format of the file.

        Returns:
            A Structure or a parent class of Structure.

        Raises:
            ValueError: If the file format is unknown.

        Examples:
            >>> structure = Structure.read("POSCAR")
            >>> poscar = Poscar.read("POSCAR")
        """
        if isinstance(file, str):
            file = pathlib.Path(file)
        if "POSCAR" in file.name.upper() or "CONTCAR" in file.name.upper() \
                or fmt == "poscar":
            return cls(*cls._from_POSCAR(file))
        else:
            raise ValueError("Unknown file format")

    @classmethod
    def _from_POSCAR(cls, file: pathlib.Path | str = "POSCAR") -> tuple:
        """Create a Structure object from a POSCAR file.

        Args:
            file: The path of the POSCAR file.

        Returns:
            A tuple containing the comment, the scale factor, the lattice,
            the elements, the coordinates, the selective dynamics and the
            velocities.
        """
        with open(file, 'r') as f:
            lines = f.readlines()
            lines = [line.strip() for line in lines]
            if lines[-1] == "":
                lines = "*#*#*".join(lines).strip("*#*#*").split("*#*#*")

        comment = lines[0]
        velocities = None
        scale = float(lines[1])
        lattice = Lattice([line.split() for line in lines[2:5]])
        lattice.matrix *= scale
        elements = lines[5].split()
        num_elements = [int(num) for num in lines[6].split()]
        elements = [element for element, num in zip(elements, num_elements)
                    for _ in range(num)]
        if not lines[7].startswith("S") or lines[7].startswith("s"):
            coords = np.array([line.split()
                               for line in lines[8:8 + sum(num_elements)]])
            coords = coords.astype(float)
            selective_dynamics = None
            if lines[7].startswith("D") or lines[7].startswith("d"):
                coords = coords @ lattice.matrix
            elif lines[7].startswith("C") or lines[7].startswith("c"):
                pass
            else:
                raise ValueError("Unknown coordinate type")
        else:
            coords = list(map(lambda x:
                              np.fromstring(x, sep=" ", dtype=np.float64),
                              lines[9:9 + sum(num_elements)]))
            coords = np.array(coords)
            selective_dynamics = [[True if s.startswith("T") else False
                                  for s in line.split()[3:]]
                                  for line in lines[9:9 + sum(num_elements)]]
            if lines[8].startswith("D") or lines[8].startswith("d"):
                coords = coords @ lattice.matrix
            elif lines[8].startswith("C") or lines[8].startswith("c"):
                pass
            else:
                raise ValueError("Unknown coordinate type")
            if len(selective_dynamics) != len(coords):
                raise ValueError("The number of selective_dynamic "
                                 "must equal to the elements or coords.")
            if len(coords) != len(elements):
                raise ValueError("The number of coords "
                                 "must equal to the elements.")
        if len(lines) > 9 + sum(num_elements):
            if lines[7].startswith("S") or lines[7].startswith("s"):
                velocities = np.array(
                    [line.split() for line in lines[10 + sum(num_elements):]],
                    dtype=np.float64)
            else:
                velocities = np.array(
                    [line.split() for line in lines[9 + sum(num_elements):]],
                    dtype=np.float64)
            if len(velocities) != len(coords):
                raise ValueError("The number of velocities "
                                 "must equal to the elements or coords.")
        return coords, elements, lattice, comment, \
            selective_dynamics, velocities
