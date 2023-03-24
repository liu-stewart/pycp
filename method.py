"""This module supply some useful operation method."""

import numpy as np
from pycp_typing import Coords


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
