"""This module contains the input classes for VASP."""

from __future__ import annotations
import pathlib
from pycp.pycp_typing import Coords, NDArray
from pycp.model.sites import Sites
from pycp.model.lattice import Lattice
from pycp.model.structure import Structure


class Poscar(Structure):
    """This class used to represent a POSCAR file."""

    def __init__(self,
                 coords: Coords,
                 elements: list[str] | str,
                 lattice: Lattice,
                 comment: str = "This is a POSCAR file.",
                 selective_dynamics: list[bool] = None,  # type: ignore
                 velocities: Coords = None):  # type: ignore
        """Initialize a POSCAR.

        Args:
            coords: The coordinates of the sites.
            elements: The elements of the sites.
            lattice: The lattice of the structure.
            comment: The comment of the POSCAR.
            selective_dynamics: The selective dynamics of the sites.
            velocities: The velocities of the sites.
        """
        super().__init__(coords, elements, lattice)
        self.__comment = comment
        self.__selective_dynamics = selective_dynamics
        self.__velocities = velocities

    @property
    def comment(self) -> str:
        """Return the comment of the POSCAR."""
        return self.__comment

    @comment.setter
    def comment(self, comment: str) -> None:
        """Set the comment of the POSCAR."""
        self.__comment = comment

    @property
    def selective_dynamics(self) -> list[bool]:
        """Return the selective dynamics of the sites."""
        return self.__selective_dynamics

    @selective_dynamics.setter
    def selective_dynamics(self, selective_dynamics: list[bool]) -> None:
        """Set the selective dynamics of the sites."""
        if selective_dynamics is not None:
            if len(selective_dynamics) != len(self):
                raise ValueError("selective dynamics must "
                                 "be the same length as the number of sites")

    @property
    def velocities(self) -> Coords:
        """Return the velocities of the sites."""
        return self.__velocities

    @velocities.setter
    def velocities(self, velocities: Coords) -> None:
        """Set the velocities of the sites."""
        if velocities is not None:
            if len(velocities) != len(self):
                raise ValueError("velocities must be the same "
                                 "length as the number of sites")
