"""This module contains some geometry functions."""
import numpy as np


def create_regular_polygon(n, length: float = 1):
    """
    Create a regular polygon.

    Args:
        n: The number of vertices.
        length: The length of the edge.

    Returns:
        the coordinates of the polygon.
    """
    r = length / (2 * np.sin(np.pi / n))
    angle = 2 * np.pi / n
    angles = np.arange(0, 2 * np.pi, angle)

    x = r * np.cos(angles)
    y = r * np.sin(angles)

    center_x = np.mean(x)
    center_y = np.mean(y)

    x -= center_x
    y -= center_y
    z = np.zeros_like(x)
    return np.array([x, y, z]).T
