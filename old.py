"""This module supply some useful physical model.

Sites class is used to store the coordinates and elements of the sites,
and the class has some useful methods to transform the coordinates and
elements.
Lattice class is used to store the lattice vectors.
Structure class is used to store the lattice and sites.It inherits the
Lattice and Sites class.
"""


from __future__ import annotations
from pycp_typing import Coords, NDArray
import numpy as np
import re
from monty.json import MSONable
from pattern import pattern_element
from method import translation, rotation


class Sites(MSONable):
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
            self.__coordinates = np.array(coordinates, dtype=np.float64)
        else:
            raise TypeError("You must supply a list or ndarray which can "
                            "convert to float.")
        if self.__coordinates.shape == (3,):
            self.__coordinates = np.expand_dims(a=self.__coordinates, axis=0)
        elif self.__coordinates.ndim == 2 and self.__coordinates.shape[1] == 3:
            pass
        else:
            raise ValueError("Your input must meet one of the following "
                             "two dimensions: (3,) or (n, 3)")

        if isinstance(elements, (np.ndarray, list)):
            self.__elements = list(elements)
        elif isinstance(elements, str):
            elements = pattern_element.findall(elements)
            self.__elements = []
            for element in elements:
                """element is a string which construct by some letters and
                numbers, such as 'Ag15', and we need to split it into two
                parts, 'Ag' and '15'. Then we should add 'Ag' to the list
                15 times."""
                if element.isalpha():
                    self.__elements.append(element)
                else:
                    element = re.split(r'(\d+)', element)
                    self.__elements.extend([element[0]] * int(element[1]))
        else:
            raise ValueError(
                "You must enter elements in the following format:\n"
                "['H','H','O'] or 'HHO' or 'H2O'.")
        if len(self.__elements) != self.__coordinates.shape[0]:
            raise ValueError("The length of elements must be equal to the "
                             "number of coordinates.")

    @property
    def coordinates(self) -> NDArray:
        """The coordinates of the sites.

        Returns:
            a ndarray with shape (n, 3).
        """
        return self.__coordinates

    @property
    def elements(self) -> list[str]:
        """The elements of the sites.

        Returns:
            a list contain element.
        """
        return self.__elements

    def translate(self, vector: Coords) -> None:
        """Translate the coordinates.

        Args:
            vector: The translation vector.
        """
        self.__coordinates = translation(self.__coordinates, vector)

    def rotate(self,
               angle: float,
               axis: Coords = [0, 0, 1],
               anchor: Coords = [0, 0, 0]) -> None:
        """Rotate the coordinates.

        Args:
            angle: The rotation angle.
            axis: The rotation axis.
            anchor: The rotation anchor.
        """
        self.__coordinates = rotation(self.__coordinates, angle, axis, anchor)

    def __repr__(self) -> str:
        """Return a string representation of the object."""
        return f"Sites(coordinates={self.__coordinates}, " \
            f"elements={self.__elements})"

    def __str__(self) -> str:
        """Return a string representation of the object."""
        return f"Sites(coordinates={self.__coordinates}, " \
               f"elements={self.__elements})"

    def __len__(self) -> int:
        """Return the number of sites."""
        return len(self.__elements)

    def __getitem__(self, index: int) -> tuple[NDArray, str]:
        """Return the coordinates and element of the site.

        Args:
            index: The index of the site.

        Returns:
            a tuple contain the coordinates and element of the site.
        """
        return self.__coordinates[index], self.__elements[index]


class Lattice(MSONable):
    """This class is used to store the lattice vectors.

    Properties:
        matrix: The lattice matrix.
        length: The length of the lattice vectors.
        angle: The angle of the lattice vectors.
        volume: The volume of the lattice.
    """

    def __init__(self, matrix: Coords) -> None:
        """Init Lattice.

        Args:
            matrix: The lattice matrix.
        """
        if isinstance(matrix, (np.ndarray, list)):
            self.__matrix = np.array(matrix, dtype=np.float64)
        else:
            raise TypeError("You must supply a list or ndarray which can "
                            "convert to float.")
        if self.__matrix.shape != (3, 3):
            raise ValueError("The lattice matrix must be a 3x3 matrix.")

    @property
    def matrix(self) -> Coords:
        """Return the matrix of the lattice."""
        return self.__matrix

    @property
    def length(self) -> Coords:
        """Return the length of the lattice."""
        return np.linalg.norm(self.__matrix, axis=1)

    @property
    def angle(self) -> np.ndarray:
        """Return the angle of the lattice."""
        angle = np.arccos(
            np.sum(self.__matrix * np.roll(self.__matrix, 1, axis=0),
                   axis=1) / (self.length * np.roll(self.length, 1)))
        angle = np.degrees(np.roll(angle, 1, axis=0))
        return angle

    @property
    def volume(self) -> float:
        """Return the volume of the lattice."""
        return np.linalg.det(self.__matrix)


site = Sites([[0, 0, 0], [0.5, 0.5, 0.5]], ['Ag', 'Ag'])
site.rotate(90)
lattice = Lattice([[1, 0, 0], [0.6, 0.6, 0], [0, 0, 1]])
