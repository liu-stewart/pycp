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


class Structure(MSONable):
    """This class is used to represent the structure of a solid.

    Properties:
        lattice: The lattice of the structure.
        sites: The sites of the structure.
    """

    def __init__(self, lattice: Lattice, sites: Sites):
        """Initialize a Structure.

        Args:
            lattice: The lattice of the structure.
            sites: The sites of the structure.
        """
        if not isinstance(lattice, Lattice):
            raise TypeError("lattice must be a Lattice object")
        if not isinstance(sites, Sites):
            raise TypeError("sites must be a Sites object")
        self.__lattice = lattice
        self.__sites = sites

    def __repr__(self) -> str:
        """Return a string representation of the structure."""
        return f"Structure({self.__lattice}, {self.__sites})"

    @property
    def lattice(self) -> Lattice:
        """Return the lattice of the structure."""
        return self.__lattice

    @property
    def sites(self) -> Sites:
        """Return the sites of the structure."""
        return self.__sites

    @lattice.setter
    def lattice(self, lattice: Lattice) -> None:
        """Set the lattice of the structure."""
        if not isinstance(lattice, Lattice):
            raise TypeError("lattice must be a Lattice object")
        self.__lattice = lattice

    @sites.setter
    def sites(self, sites: Sites) -> None:
        """Set the sites of the structure."""
        if not isinstance(sites, Sites):
            raise TypeError("sites must be a Sites object")
        self.__sites = sites
