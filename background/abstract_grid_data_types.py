from dataclasses import dataclass
from enum import Enum
import math
import numpy as np
import warnings

from mashumaro.mixins.json import DataClassJSONMixin
from pandera import DataFrameModel
from pandera.typing import Index, Series
from pandera.typing.geopandas import GeoSeries
from pyproj import CRS


class GridGeometry(str, Enum):
    SQUARE = "square"
    HEXAGON = "hexagon"

    def __str__(self):
        return self.name


@dataclass
class AbstractGridGDF(DataFrameModel):
    """
    Pandera model for a GeoDataFrame containing geometries for an arbitrary set of AbstractGrid
    grid-cells.
    """

    grid_x: Index[int]
    grid_y: Index[int]
    geometry: GeoSeries


@dataclass
class GridIDAbstractGridGDF(DataFrameModel):
    """
    Pandera model for a GeoDataFrame containing geometries for an arbitrary set of
    BoundedAbstractGrid grid-cells.
    """

    grid_x: Index[int]
    grid_y: Index[int]
    grid_id: Series[int]
    geometry: GeoSeries


@dataclass
class AbstractGrid(DataClassJSONMixin):
    """
    Specification of an AbstractGrid.
    """

    # EPSG code of the UTM projection being used
    epsg_code: int

    # South-west corner (origin) of the grid, units are the CRS of the UTM projection
    # (easting/northing in m)
    origin_x: float
    origin_y: float

    # Hexagon or square grid
    geometry_type: GridGeometry

    # Area of each cell of the grid
    # Unit: m^2 for UTM CRS's. Only UTM CRS's are intended to be used with ObservationGrid!
    area: float

    def __post_init__(self):
        # Max area: 100km x 100km = 10_000 km^2
        AREA_LIMIT = 100e3**2

        crs = CRS.from_epsg(self.epsg_code)

        if crs.utm_zone is None:
            raise ValueError(f"Must provide a UTM CRS. {self.epsg_code} is non-UTM.")

        if crs.ellipsoid.name != "WGS 84":
            raise ValueError(f"Your CRS must use the WGS 84 ellipse ({crs.ellipsoid.name} invalid)")

        # Modify the grid geometry to allow loading from json.
        # [DP-4788] TODO: Int support should be removed after transition.
        if not isinstance(self.geometry_type, GridGeometry):
            if isinstance(self.geometry_type, int):
                self.geometry_type = (
                    GridGeometry.SQUARE if self.geometry_type == 1 else GridGeometry.HEXAGON
                )
            elif isinstance(self.geometry_type, str):
                self.geometry_type = GridGeometry(self.geometry_type)
            else:
                raise TypeError("geometry_type must be a GridGeometry")

        if self.area <= 0:
            raise ValueError(f"area must be positive, but was: {self.area}")

        if self.area > AREA_LIMIT:
            raise ValueError("grid-cell area is way too big, greater than 100km * 100km")


@dataclass
class BoundedAbstractGrid(AbstractGrid):
    """
    Bounded variant of an AbstractGrid. Any AbstractGrid parameterization is acceptable - hexagons,
    squares, area, origin, etc.
    """

    # North-east corner (extent) of the grid
    #
    # Units: [meters] (meters east/north, in the CRS of the UTM projection)
    extent_x: float
    extent_y: float

    def __post_init__(self):
        super().__post_init__()  # call `AbstractGrid`'s __post_init__()

        if self.extent_x <= self.origin_x:
            raise ValueError("extent_x was West of origin_x, or equal to origin_x")

        if self.extent_y <= self.origin_y:
            raise ValueError("extent_y was South of origin_y, or equal to origin_y")

        if self.extent_x <= self.origin_x + math.sqrt(self.area):
            warnings.warn(
                "extent_x was only about 1 grid-cell east of origin_x. Consider expanding?"
            )

        if self.extent_y <= self.origin_y + math.sqrt(self.area):
            warnings.warn(
                "extent_y was only about 1 grid-cell north of origin_y. Consider expanding?"
            )


@dataclass
class AbstractGridCoordinates():
    """
    AbstractGridCoordinates links a grid and a set of coordinates on the grid.
    It helps handling a specific set of coordinates on a grid, like the reported
    emitters on the emitter grid or the observations on the observation grid.
    TODO: Apply this dataclass thouroughy to the codebase.

    Contains:
        - grid: The grid on which the coordinates are defined
        - coords: The coordinates on the grid, shape (... , 2) of (x,y) pairs.
    """
    grid: AbstractGrid
    coords: np.ndarray
