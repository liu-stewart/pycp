"""This module is used to store the class Lattice."""

import numpy as np
from pycp.pycp_typing import Coords
from monty.json import MSONable


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
    def matrix(self) -> np.ndarray:
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

    @matrix.setter
    def matrix(self, matrix: Coords) -> None:
        """Set the matrix of the lattice."""
        if isinstance(matrix, (np.ndarray, list)):
            self.__matrix = np.array(matrix, dtype=np.float64)
        else:
            raise TypeError("You must supply a list or ndarray which can "
                            "convert to float.")
        if self.__matrix.shape != (3, 3):
            raise ValueError("The lattice matrix must be a 3x3 matrix.")

    @classmethod
    def from_lengths_and_angles(cls, lengths: Coords, angles: Coords
                                ) -> "Lattice":
        """Generate a lattice from lengths and angles.

        Args:
            lengths: The length of the lattice vectors.
            angles: The angle of the lattice vectors.

        Returns:
            A Lattice object.
        """
        if isinstance(lengths, (np.ndarray, list)):
            lengths = np.array(lengths, dtype=np.float64)
        else:
            raise TypeError("You must supply a list or ndarray which can "
                            "convert to float.")
        if isinstance(angles, (np.ndarray, list)):
            angles = np.array(angles, dtype=np.float64)
        else:
            raise TypeError("You must supply a list or ndarray which can "
                            "convert to float.")
        if lengths.shape != (3,):
            raise ValueError("The lengths must be a 3x1 vector.")
        if angles.shape != (3,):
            raise ValueError("The angles must be a 3x1 vector.")
        angles = np.radians(angles)
        matrix = np.zeros((3, 3), dtype=np.float64)
        matrix[0, 0] = lengths[0]
        matrix[1, 0] = lengths[1] * np.cos(angles[2])
        matrix[1, 1] = lengths[1] * np.sin(angles[2])
        matrix[2, 0] = lengths[2] * np.cos(angles[1])
        matrix[2, 1] = lengths[2] * (np.cos(angles[0]) - np.cos(angles[1])
                                     * np.cos(angles[2])) / np.sin(angles[2])
        matrix[2, 2] = np.sqrt(lengths[2] ** 2 - matrix[2, 0] ** 2 -
                               matrix[2, 1] ** 2)
        return cls(matrix)

    def __repr__(self) -> str:
        """Return the representation of the lattice."""
        return "Lattice(matrix=\n{})".format(self.__matrix)

    def __str__(self) -> str:
        """Return the string of the lattice."""
        return "Lattice(matrix=\n{})".format(self.__matrix)
