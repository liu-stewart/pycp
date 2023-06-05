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
import pathlib
import copy


class Structure(Sites, Lattice):
    """This class is used to represent the structure of a solid.

    Properties:
        lattice: The lattice of the structure.
        cartesian_coords: The cartesian coordinates of the sites.
        fractional_coords: The fractional coordinates of the sites.

    Methods:
        periodic_boundary_conditions: Apply periodic boundary conditions to the
            fractional coordinates.
        rotate: Rotate the structure.
        translate: Translate the structure.
        axial_symmetry: Apply axial symmetry to the structure.
        read: Read a structure from a file.
    """

    def __init__(self,
                 coordinates: Coords,
                 elements: list[str] | str,
                 matrix: Coords,
                 comment: str = "Default comment",
                 selective_dynamics: list[list[bool]] = None,  # type: ignore
                 velocities: Coords = None):  # type: ignore
        """Initialize a Structure.

        Args:
            coordinates: The coordinates of the sites.
            elements: The elements of the sites.
            matrix: The lattice of the structure.
            comment: The comment of the POSCAR.
            selective_dynamics: The selective dynamics of the sites.
            velocities: The velocities of the sites.
        """
        Sites.__init__(self, coordinates, elements)
        Lattice.__init__(self, matrix)
        self._comment = comment
        if selective_dynamics is not None:
            self._selective_dynamics = selective_dynamics
        else:
            self._selective_dynamics = [[True, True, True]
                                        for i in range(len(self))]
        if velocities is not None:
            self._velocities = velocities
        else:
            self._velocities = np.array([[0, 0, 0] for i in range(len(self))])
        self.periodic_boundary_conditions()
        self.init_md = False

    @property
    def fractional_coords(self) -> NDArray:
        """Return the fractional coordinates of the sites."""
        return self.coordinates @ np.linalg.inv(self.matrix)

    @property
    def selective_dynamics(self) -> list[list[bool]]:
        """Return the selective dynamics of the sites."""
        return self._selective_dynamics

    @property
    def comment(self) -> str:
        """Return the comment of the POSCAR."""
        return self._comment

    @comment.setter
    def comment(self, comment: str) -> None:
        """Set the comment of the POSCAR."""
        self._comment = comment

    @selective_dynamics.setter
    def selective_dynamics(self, selective_dynamics: list[list[bool]]) -> None:
        """Set the selective dynamics of the sites."""
        if selective_dynamics is not None:
            if len(selective_dynamics) != len(self):
                raise ValueError("selective dynamics must "
                                 "be the same length as the number of sites")
            self._selective_dynamics = selective_dynamics
        else:
            self._selective_dynamics = [[True, True, True]
                                        for i in range(len(self))]

    @property
    def velocities(self) -> Coords:
        """Return the velocities of the sites."""
        return self._velocities

    @velocities.setter
    def velocities(self, velocities: Coords) -> None:
        """Set the velocities of the sites."""
        if velocities is not None:
            if len(velocities) != len(self):
                raise ValueError("velocities must be the same "
                                 "length as the number of sites")
            self._velocities = velocities
        else:
            self._velocities = np.array([[0, 0, 0] for i in range(len(self))])

    def periodic_boundary_conditions(self) -> None:
        """Apply periodic boundary conditions to the fractional coordinates."""
        coords = self.fractional_coords
        coords -= np.floor(coords)
        self._coordinates = coords @ self.matrix

    def __repr__(self) -> str:
        """Return the representation of the structure."""
        sites = super(Structure, self).__str__()
        lattice = super(Sites, self).__str__()
        return sites + f"\n\n{lattice}"

    def __str__(self) -> str:
        """Return the string representation of the structure."""
        return self.string_POSCAR()

    @classmethod
    def read(cls, from_file: pathlib.Path | str, fmt: str = "poscar"):
        """Create a Structure object from a file.

        Args:
            file: The path of the file.
            fmt: The format of the file.

        Returns:
            A Structure or a parent class of Structure.

        Raises:
            ValueError: If the file format is unknown.

        Examples:
            >>> structure = Structure.read("POSCAR")
            >>> poscar = Poscar.read("POSCAR")
        """
        if isinstance(from_file, str):
            from_file = pathlib.Path(from_file)
        if "POSCAR" in from_file.name.upper() or \
            "CONTCAR" in from_file.name.upper() or \
                fmt == "poscar":
            return cls(*cls._from_POSCAR(from_file))
        else:
            raise ValueError("Unknown file format")

    @classmethod
    def _from_POSCAR(cls, from_file: pathlib.Path | str = "POSCAR") -> tuple:
        """Create a Structure object from a POSCAR file.

        Args:
            file: The path of the POSCAR file.

        Returns:
            A tuple containing the comment, the scale factor, the lattice,
            the elements, the coordinates, the selective dynamics and the
            velocities.
        """
        with open(from_file, 'r') as f:
            lines = f.readlines()
            lines = [line.strip() for line in lines]
            if lines[-1] == "":
                lines = "*#*#*".join(lines).strip("*#*#*").split("*#*#*")

        comment = lines[0]
        velocities = None
        scale = float(lines[1])
        matrix = np.array([line.split() for line in lines[2:5]],
                          dtype=np.float64)
        matrix *= scale
        elements = lines[5].split()
        num_elements = [int(num) for num in lines[6].split()]
        elements = [element for element, num in zip(elements, num_elements)
                    for _ in range(num)]
        if not lines[7].startswith("S") or lines[7].startswith("s"):
            coords = np.array([line.split()[:3]
                               for line in lines[8:8 + sum(num_elements)]])
            coords = coords.astype(float)
            selective_dynamics = None
            if lines[7].startswith("D") or lines[7].startswith("d"):
                coords = coords @ matrix
            elif lines[7].startswith("C") or lines[7].startswith("c"):
                pass
            else:
                raise ValueError("Unknown coordinate type")
        else:
            coords = list(map(lambda x:
                              np.fromstring(x, sep=" ", dtype=np.float64),
                              lines[9:9 + sum(num_elements)]))
            coords = np.array(coords)
            selective_dynamics = [[True if s.startswith("T") else False
                                  for s in line.split()[3:]]
                                  if len(line.split()) == 6
                                  else
                                  [True if s.startswith("T") else False
                                   for s in line.split()[4:]]
                                  for line in lines[9:9 + sum(num_elements)]]
            if lines[8].startswith("D") or lines[8].startswith("d"):
                coords = coords @ matrix
            elif lines[8].startswith("C") or lines[8].startswith("c"):
                pass
            else:
                raise ValueError("Unknown coordinate type")
            if len(selective_dynamics) != len(coords):
                raise ValueError("The number of selective_dynamic "
                                 "must equal to the elements or coords.")
            if len(coords) != len(elements):
                raise ValueError("The number of coords "
                                 "must equal to the elements.")
        if len(lines) > 9 + sum(num_elements):
            if lines[7].startswith("S") or lines[7].startswith("s"):
                velocities = np.array(
                    [line.split() for line in lines[10 + sum(num_elements):]],
                    dtype=np.float64)
            else:
                velocities = np.array(
                    [line.split() for line in lines[9 + sum(num_elements):]],
                    dtype=np.float64)
            if len(velocities) != len(coords):
                raise ValueError("The number of velocities "
                                 "must equal to the elements or coords.")
        return coords, elements, matrix, comment, \
            selective_dynamics, velocities

    def remove_duplicates(self, tol: float = 1e-2) -> None:
        """Remove duplicate sites.

        Args:
            tol: The tolerance of the distance between two sites.
        """
        self.periodic_boundary_conditions()
        coords = self.coordinates
        elements = self.elements
        unique_coords = []
        unique_elements = []
        velocities = self.velocities
        selective_dynamics = self.selective_dynamics
        unique_velocities = []
        unique_selective_dynamics = []
        for coord, element, velocity, selective_dynamic \
                in zip(coords, elements, velocities, selective_dynamics):
            if not any(self.is_same_coord(coord, unique_coord, tol=tol)
                       for unique_coord in unique_coords):
                unique_coords.append(coord)
                unique_elements.append(element)
                unique_velocities.append(velocity)
                unique_selective_dynamics.append(selective_dynamic)
        self.coordinates = np.array(unique_coords)
        self.elements = unique_elements
        self.velocities = np.array(unique_velocities)
        self.selective_dynamics = unique_selective_dynamics

    def sort(self) -> None:
        """Sort the sites by the elements."""
        all = list(zip(self.coordinates, self.elements,
                       self.velocities, self.selective_dynamics))
        all.sort(key=lambda x: x[1])
        all = list(zip(*all))
        self.coordinates = np.array(all[0])
        self.elements = list(all[1])  # type: ignore
        self.velocities = np.array(all[2])
        self.selective_dynamics = list(all[3])  # type: ignore

    def select(self, elements: list[str] | str) -> Structure:
        """Select sites by the elements.

        Args:
            elements: The elements to select.
        """
        if isinstance(elements, str):
            elements = [elements]
        all = list(zip(self.coordinates, self.elements,
                       self.velocities, self.selective_dynamics))
        all = [x for x in all if x[1] in elements]
        all = list(zip(*all))
        coordinates = np.array(all[0])
        elements = list(all[1])  # type: ignore
        velocities = np.array(all[2])
        selective_dynamics = list(all[3])

        return self.__class__(coordinates, elements, self.matrix, self.comment,
                              selective_dynamics, velocities)  # type: ignore

    def __add__(self,
                other: Sites | tuple[Coords, str] | Structure):
        """Add sites to the structure.

        Args:
            other: The sites to add.
        """
        if isinstance(other, (Sites, tuple, list)):
            if isinstance(other, (tuple, list)):
                other = Sites(other[0], other[1])
            coords = other.coordinates
            elements = other.elements
            other = Structure(coords, elements, self.matrix)
        coords = np.concatenate((self.coordinates, other.coordinates))
        elements = self.elements + other.elements
        velocities = np.concatenate((self.velocities, other.velocities))
        selective_dynamics = self.selective_dynamics + other.selective_dynamics
        return self.__class__(coords, elements, self.matrix,
                              self.comment, selective_dynamics, velocities)

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
            self._velocities = np.delete(self._velocities, i, axis=0)
            self._selective_dynamics.pop(i)

    def copy(self) -> Structure:
        """Return a copy of the structure."""
        return copy.deepcopy(self)

    def supercell(self, n: int | tuple[int, int, int] | list[int]
                  ) -> Structure:
        """Make a supercell of the structure.

        Args:
            n: The number of supercells in each direction.
        """
        if isinstance(n, int):
            n = (n, n, n)
        coords = []
        elements = []
        velocities = []
        selective_dynamics = []
        matrix = self.matrix * n
        for i in range(n[0]):
            for j in range(n[1]):
                for k in range(n[2]):
                    for coord, element, velocity, selective_dynamic \
                            in zip(self.coordinates, self.elements,
                                   self.velocities, self.selective_dynamics):
                        coords.append(coord + i * self.matrix[0] +
                                      j * self.matrix[1] +
                                      k * self.matrix[2])
                        elements.append(element)
                        velocities.append(velocity)
                        selective_dynamics.append(selective_dynamic)
        coords = np.array(coords)
        velocities = np.array(velocities)
        return self.__class__(coords, elements, matrix,
                              self.comment, selective_dynamics, velocities)

    def base_transformation(self, vectors):
        """Transform the structure to a new base.

        Args:
            base: The new base.
        """
        vectors = np.array(vectors[:2])
        vectors = np.row_stack((vectors, vectors[0] + vectors[1], np.zeros(3)))
        fvectors = vectors @ np.linalg.inv(self.matrix)
        ran = np.array([max(component)-min(component)
                        for component in fvectors.T])
        ran = np.ceil(ran).astype(int)
        ran[2] = 1
        structure = self.supercell(ran)  # type: ignore
        vectors = np.row_stack((vectors[:2], self.matrix[2]))
        structure.matrix = vectors
        structure.remove_duplicates()
        return structure

    def write(self, to_file: pathlib.Path | str = "POSCAR", fmt="poscar"):
        """Write the POSCAR to a file."""
        if isinstance(to_file, str):
            to_file = pathlib.Path(to_file)
        if not to_file.parent.exists():
            to_file.parent.mkdir(parents=True)
        if fmt == "poscar":
            with open(to_file, "w+") as file:
                file.write(self.string_POSCAR())

    def string_POSCAR(self) -> str:
        """Return the string of the POSCAR file."""
        comment = self.comment
        matrix = self.matrix
        elements = self.elements
        independent_elements = ['start']
        elements_count = []
        for element in elements:
            if element != independent_elements[-1]:
                independent_elements.append(element)
                elements_count.append(1)
            else:
                elements_count[-1] += 1
        independent_elements = independent_elements[1:]
        selective_dynamics = [["T" if i else "F" for i in j]
                              for j in self.selective_dynamics]
        string = f"{comment}\n" + \
            f"   1\n" + \
            f"    {matrix[0][0]:.16f}" + \
            f" {matrix[0][1]:.16f} {matrix[0][2]:.16f}\n" + \
            f"    {matrix[1][0]:.16f}" + \
            f" {matrix[1][1]:.16f} {matrix[1][2]:.16f}\n" + \
            f"    {matrix[2][0]:.16f}" + \
            f" {matrix[2][1]:.16f} {matrix[2][2]:.16f}\n" + \
            f"   {'    '.join(independent_elements)}\n" + \
            f"   {'    '.join([str(i) for i in elements_count])}\n" + \
            f"Selective Dynamics\n" + \
            f"Direct\n"
        for i in range(len(self)):
            string += f"  {self.fractional_coords[i][0]:.16f}" + \
                f" {self.fractional_coords[i][1]:.16f}" + \
                f" {self.fractional_coords[i][2]:.16f}"
            string += f"  {selective_dynamics[i][0]}" + \
                f" {selective_dynamics[i][1]}" + \
                f" {selective_dynamics[i][2]}"
            string += "\n"
        if self.init_md is False:
            string += "\n"
            for i in range(len(self)):
                string += f"  {self.velocities[i][0]:.16f}" + \
                    f" {self.velocities[i][1]:.16f}" + \
                    f" {self.velocities[i][2]:.16f}\n"
        return string
