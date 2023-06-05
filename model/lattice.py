"""This module is used to store the class Lattice."""
import numpy as np
from pycp.pycp_typing import Coords


class Lattice():
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
            self._matrix = np.array(matrix, dtype=np.float64)
        else:
            raise TypeError("You must supply a list or ndarray which can "
                            "convert to float.")
        if self._matrix.shape != (3, 3):
            raise ValueError("The lattice matrix must be a 3x3 matrix.")

    @property
    def matrix(self) -> np.ndarray:
        """Return the matrix of the lattice."""
        return self._matrix

    @matrix.setter
    def matrix(self, matrix: Coords) -> None:
        """Set the matrix of the lattice."""
        if isinstance(matrix, (np.ndarray, list)):
            self._matrix = np.array(matrix, dtype=np.float64)
        else:
            raise TypeError("You must supply a list or ndarray which can "
                            "convert to float.")
        if self._matrix.shape != (3, 3):
            raise ValueError("The lattice matrix must be a 3x3 matrix.")
        if type(self) is not Lattice:
            self.periodic_boundary_conditions()  # type: ignore

    @property
    def length(self) -> Coords:
        """Return the length of the lattice."""
        return np.linalg.norm(self._matrix, axis=1)

    @property
    def angle(self) -> np.ndarray:
        """Return the angle of the lattice."""
        angle = np.arccos(
            np.sum(self._matrix * np.roll(self._matrix, 1, axis=0),
                   axis=1) / (self.length * np.roll(self.length, 1)))
        angle = np.degrees(np.roll(angle, 1, axis=0))
        return angle

    @property
    def volume(self) -> float:
        """Return the volume of the lattice."""
        return np.linalg.det(self._matrix)

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

    def __str__(self) -> str:
        """Return the string of the lattice."""
        return "Lattice(matrix=\n{})".format(self._matrix)

    def __getitem__(self, index) -> Coords:
        """Return the vector of the lattice."""
        return self.matrix[index]

    def __setitem__(self, index, vector: Coords) -> None:
        """Set the vector of the lattice."""
        self.matrix[index] = vector

    def distance(self, coord1: Coords, coord2: Coords) -> float:
        """Calculate the distance between two sites.

        Consider the periodic boundary conditions.

        Args:
            site1: The first site.
            site2: The second site.

        Returns:
            The distance between two sites.
        """
        coord1 = np.array(coord1)
        coord2 = np.array(coord2)
        lattice = self.matrix
        dist = np.linalg.norm(coord1 - coord2)  # type: ignore
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    a = np.linalg.norm(coord1 - coord2 +  # type: ignore
                                       i * lattice[0] +
                                       j * lattice[1] +
                                       k * lattice[2])
                    dist = min(dist, a)  # type: ignore
        return dist

    def is_same_coord(self, coord1: Coords, coord2: Coords, tol=1e-2) -> bool:
        """Check if two coordinates are the same.

        Consider the periodic boundary conditions.

        Args:
            coord1: The first coordinate.
            coord2: The second coordinate.

        Returns:
            True if two coordinates are the same.
        """
        coord1 = np.array(coord1)
        coord2 = np.array(coord2)
        matrix = self.matrix
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    if np.allclose(coord1, coord2 +  # type: ignore
                                   i * matrix[0] +
                                   j * matrix[1] +
                                   k * matrix[2], atol=tol):
                        return True
        return False
