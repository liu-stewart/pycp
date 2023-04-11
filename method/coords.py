"""This module supply some useful operation method."""
import numpy as np
from pycp.pycp_typing import Coords


def translation(coords: Coords, vector: Coords) -> np.ndarray:
    """Translate the coordinates.

    Args:
        coords: The coordinates to be translated.
        vector: The translation vector.

    Returns:
        The translated coordinates.
    """
    coords = np.array(coords, dtype=np.float64)
    vector = np.array(vector, dtype=np.float64)
    return coords + vector


def rotation(coords: Coords,
             angle: float,
             axis: Coords = [0, 0, 1],
             anchor: Coords = [0, 0, 0],
             ) -> np.ndarray:
    """Rotate the coordinates.

    Args:
        coords: The coordinates to be rotated.
        angle: The rotation angle in degree.
        axis: The rotation axis.
        anchor: The rotation anchor.

    Returns:
        The rotated coordinates.
    """
    coords = np.array(coords, dtype=np.float64)
    anchor = np.array(anchor, dtype=np.float64)
    axis = np.array(axis, dtype=np.float64)
    angle = np.deg2rad(angle)
    coords = coords - anchor
    c = np.cos(angle)
    s = np.sin(angle)
    t = 1 - c
    x, y, z = axis
    R = np.array([[t * x * x + c, t * x * y - s * z, t * x * z + s * y],
                  [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
                  [t * x * z - s * y, t * y * z + s * x, t * z * z + c]])
    return np.dot(R, coords.T).T + anchor


def axial_symmetry(coords: Coords,
                   anchor1: Coords,
                   anchor2: Coords = [0, 0, 0],
                   ) -> np.ndarray:
    """Axial symmetry the coordinates.

    Args:
        coords: The coordinates to be axial symmetry.
        anchor1: The first anchor.
        anchor2: The second anchor.

    Returns:
        The axial symmetry coordinates.
    """
    coords = np.array(coords, dtype=np.float64)
    coords = coords - anchor2
    anchor1 = np.array(anchor1, dtype=np.float64)
    anchor2 = np.array(anchor2, dtype=np.float64)
    axial_vector = (anchor1 - anchor2) / np.linalg.norm(anchor1 - anchor2)
    axial_vector = np.expand_dims(axial_vector, axis=1)
    vector = np.dot(coords, axial_vector) * axial_vector.T
    coords = coords + 2 * (vector - coords)
    return coords + anchor2


def temp_mirror(coords: Coords,
                anchor1,
                anchor2):
    """Temp mirror the coordinates.

    Args:
        coords: The coordinates to be temp mirror.
        anchor1: The first anchor.
        anchor2: The second anchor.

    Returns:
        The temp mirror coordinates.
    """
    ascoords = axial_symmetry(coords, anchor1, anchor2)
    for index in range(len(ascoords)):
        ascoords[index][2] = coords[index][2]
    return ascoords


def perturbation(coords, scope=0.1):
    """Perturb the coordinates.

    Args:
        coords: The coordinates to be perturbed.
        scope: The perturbation scope.

    Returns:
        The perturbed coordinates.
    """
    coords = np.array(coords, dtype=np.float64)
    return coords + np.random.uniform(-scope, scope, coords.shape)
