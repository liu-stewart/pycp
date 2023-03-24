"""This module supply some typing."""

from __future__ import annotations
from numpy.typing import NDArray
from typing import Union
from pymatgen.core.structure import Structure

Coords = Union[NDArray, list]
