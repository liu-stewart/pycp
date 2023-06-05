"""This module supply the class Sites.

Sites class is used to store the coordinates and elements of the sites,
and the class has some useful methods to transform the coordinates and
elements.
"""
from __future__ import annotations
import numpy as np
import re
from pycp.pycp_typing import Coords, NDArray
from pycp.method.coords import translation, rotation
from pycp.method.coords import axial_symmetry, perturbation
from pycp.method.coords import temp_mirror
from pycp.pattern import pattern_element
import copy


class Sites():
    """This class is used to store the coordinates and elements of the sites.

    Properties:
        coordinates: The coordinates of the sites.
        elements: The elements of the sites.
    """

    def __init__(self, coordinates: Coords, elements: list[str] | str) -> None:
        """Init Sites.

        Args:
            coordinates:
                a list or ndarray with type (n, 3) or (3,).
            elements:
                a list contain element or a string.
                example:
                ['Ag', 'Ag'] | 'Ag2' | 'AgAg'
        """
        if isinstance(coordinates, (np.ndarray, list)):
            self._coordinates = np.array(coordinates, dtype=np.float64)
        else:
            raise TypeError("You must supply a list or ndarray which can "
                            "convert to float.")
        if self._coordinates.shape == (3,):
            self._coordinates = np.expand_dims(a=self._coordinates, axis=0)
        elif self._coordinates.ndim == 2 and self._coordinates.shape[1] == 3:
            pass
        else:
            raise ValueError("Your input must meet one of the following "
                             "two dimensions: (3,) or (n, 3)")

        if isinstance(elements, (np.ndarray, list)):
            self._elements = list(elements)
        elif isinstance(elements, str):
            elements = pattern_element.findall(elements)
            self._elements = []
            for element in elements:
                if element.isalpha():
                    self._elements.append(element)
                else:
                    element = re.split(r'(\d+)', element)
                    self._elements.extend([element[0]] * int(element[1]))
        else:
            raise ValueError(
                "You must enter elements in the following format:\n"
                "['H','H','O'] or 'HHO' or 'H2O'.")
        if len(self._elements) != self._coordinates.shape[0]:
            raise ValueError("The length of elements must be equal to the "
                             "number of coordinates.")

    @property
    def coordinates(self) -> NDArray:
        """The coordinates of the sites.

        Returns:
            a ndarray with shape (n, 3).
        """
        return self._coordinates

    @coordinates.setter
    def coordinates(self, coordinates: Coords) -> None:
        """Set the coordinates of the sites.

        Args:
            coordinates:
                a list or ndarray with type (n, 3) or (3,).
        """
        if isinstance(coordinates, (np.ndarray, list)):
            self._coordinates = np.array(coordinates, dtype=np.float64)
        else:
            raise TypeError("You must supply a list or ndarray which can "
                            "convert to float.")
        if self._coordinates.shape == (3,):
            self._coordinates = np.expand_dims(a=self._coordinates, axis=0)
        elif self._coordinates.ndim == 2 and self._coordinates.shape[1] == 3:
            pass
        else:
            raise ValueError("Your input must meet one of the following "
                             "two dimensions: (3,) or (n, 3)")
        if type(self) is not Sites:
            self.periodic_boundary_conditions()  # type: ignore

    @property
    def elements(self) -> list[str]:
        """The elements of the sites.

        Returns:
            a list contain element.
        """
        return self._elements

    @elements.setter
    def elements(self, elements: list[str] | str) -> None:
        """Set the elements of the sites.

        Args:
            elements:
                a list contain element or a string.
                example:
                ['Ag', 'Ag'] | 'Ag2' | 'AgAg'
        """
        if isinstance(elements, (np.ndarray, list)):
            self._elements = list(elements)
        elif isinstance(elements, str):
            elements = pattern_element.findall(elements)
            self._elements = []
            for element in elements:
                if element.isalpha():
                    self._elements.append(element)
                else:
                    element = re.split(r'(\d+)', element)
                    self._elements.extend([element[0]] * int(element[1]))
        else:
            raise ValueError(
                "You must enter elements in the following format:\n"
                "['H','H','O'] or 'HHO' or 'H2O'.")

    def translate(self, vector: Coords) -> Sites:
        """Translate the coordinates.

        Args:
            vector: The translation vector.
        """
        self.coordinates = translation(self.coordinates, vector)
        return self

    def rotate(self,
               angle: float,
               axis: Coords = [0, 0, 1],
               anchor: Coords = [0, 0, 0]) -> Sites:
        """Rotate the coordinates.

        Args:
            angle: The rotation angle.
            axis: The rotation axis.
            anchor: The rotation anchor.
        """
        self.coordinates = rotation(self.coordinates, angle, axis, anchor)
        return self

    def axial_symmetry(self,
                       anchor1: Coords,
                       anchor2: Coords = [0, 0, 0]) -> Sites:
        """Axial symmetry the coordinates.

        Args:
            anchor1: The first anchor.
            anchor2: The second anchor.
        """
        self.coordinates = axial_symmetry(self.coordinates, anchor1, anchor2)
        return self

    def perturbation(self, scope: float = 0.1) -> Sites:
        """Perturbation the coordinates.

        Args:
            scope: The standard deviation of the normal distribution.
        """
        self.coordinates = perturbation(self.coordinates, scope)
        return self

    def temp_mirror(self,
                    anchor1: Coords,
                    anchor2: Coords) -> Sites:
        """Temp mirror the coordinates.

        Args:
            anchor1: The first anchor.
            anchor2: The second anchor.
        """
        self.coordinates = temp_mirror(self.coordinates, anchor1, anchor2)
        return self

    def __str__(self) -> str:
        """Return a string representation of the object."""
        return f"Sites(coordinates=\n{self.coordinates}, " \
               f"\n\nelements=\n{self.elements})\n"

    def __len__(self) -> int:
        """Return the number of sites."""
        return len(self.coordinates)

    def __getitem__(self, index: int) -> Sites:
        """Return the coordinates and element of the site.

        Args:
            index: The index of the site.

        Returns:
            a tuple contain the coordinates and element of the site.
        """
        return Sites(self.coordinates[index], self.elements[index])

    def __setitem__(self, index: int, value: tuple[Coords, str]) -> None:
        """Set the coordinates and element of the site.

        Args:
            index: The index of the site.
            value: The coordinates and element of the site.
        """
        self._coordinates[index] = value[0]
        self._elements[index] = value[1]

    def __add__(self, other: Sites | tuple[Coords, str]) -> Sites:
        """Add a site to the Sites.

        Args:
            other: The site which will be added.

        Returns:
            a new Sites object.
        """
        if isinstance(other, Sites):
            coordinates = np.concatenate((self.coordinates,
                                          other.coordinates))
            elements = self.elements + other.elements
        elif isinstance(other, (tuple, list)):
            coordinates = np.concatenate((self.coordinates,
                                          np.expand_dims(other[0], axis=0)))
            elements = self.elements + [other[1]]
        else:
            raise TypeError("You must supply a Sites or a tuple.")
        return Sites(coordinates, elements)

    def remove(self, index: int | list[int]) -> None:
        """Remove a site from the Sites.

        Args:
            index: The index of the site.
        """
        if isinstance(index, int):
            index = [index]
        for i in sorted(index, reverse=True):
            self._coordinates = np.delete(self._coordinates, i, axis=0)
            self._elements.pop(i)

    def copy(self) -> Sites:
        """Return a copy of the Sites."""
        return copy.deepcopy(self)

    def __iter__(self) -> Iterator[Sites]:
        """Return an iterator of the Sites."""
        self.__iterindex = 0
        return self

    def __next__(self) -> Sites:
        """Return the next site."""
        if self.__iterindex < len(self):
            self.__iterindex += 1
            return self[self.__iterindex - 1]
        else:
            raise StopIteration
