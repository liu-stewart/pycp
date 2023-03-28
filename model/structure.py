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
        sites: The sites of the structure.
    """

    def __init__(self,
                 coordinates: Coords,
                 elements: list[str] | str,
                 lattice: Lattice):
        """Initialize a Structure.

        Args:
            lattice: The lattice of the structure.
            sites: The sites of the structure.
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
        return super().__repr__() + f" with lattice {self.lattice}"

    def __str__(self) -> str:
        """Return the string representation of the structure."""
        return super().__str__() + f" with lattice {self.lattice}"

    @classmethod
    def from_POSCAR(cls, file: pathlib.Path | str = "POSCAR"):
        """Create a Structure object from a POSCAR file.

        Args:
            file: The path of the POSCAR file.
        """
        with open(file, 'r') as f:
            lines = f.readlines()
            lines = [line.strip() for line in lines]
        zoom = float(lines[1])
        lattice = Lattice([line.split() for line in lines[2:5]])
        lattice.matrix *= zoom
        elements = lines[5].split()
        num_elements = [int(num) for num in lines[6].split()]
        elements = [element for element, num in zip(elements, num_elements)
                    for _ in range(num)]
        coords = np.array([line.split()
                           for line in lines[8:8 + sum(num_elements)]])
        coords = coords.astype(float)
        if lines[7].startswith("D") or lines[7].startswith("d"):
            coords = coords @ lattice.matrix
        elif lines[7].startswith("C") or lines[7].startswith("c"):
            pass
        else:
            raise ValueError("Unknown coordinate type")
        return cls(coords, elements, lattice)
