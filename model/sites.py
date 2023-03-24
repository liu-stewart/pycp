"""This module supply the class Sites.

Sites class is used to store the coordinates and elements of the sites,
and the class has some useful methods to transform the coordinates and
elements.
"""


from __future__ import annotations
import numpy as np
import re
from monty.json import MSONable
from pycp.pycp_typing import Coords, NDArray
from pycp.method import translation, rotation
from pycp.pattern import pattern_element


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
