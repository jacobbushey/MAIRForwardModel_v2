from typing import Optional, Union

import geopandas as gpd
import numpy as np
import pandas as pd
import pyproj
from shapely import box, polygons

#from msat_level4_gim.abstract_grid.coords import generate_tiled_grid_id_coord, convert_grid_to_utm
#from msat_level4_gim.abstract_grid.data_types import (
#    AbstractGrid,
#    AbstractGridGDF,
#    BoundedAbstractGrid,
#    GridGeometry,
#    GridIDAbstractGridGDF,
#)
from coords import generate_tiled_grid_id_coord, convert_grid_to_utm
from data_types import (
    AbstractGrid,
    AbstractGridGDF,
    BoundedAbstractGrid,
    GridGeometry,
    GridIDAbstractGridGDF,
)

def _make_hexagonal_abstract_grid_gdf(
    grid: AbstractGrid,
    max_utm_x: int,
    max_utm_y: int,
) -> AbstractGridGDF:
    """
    Generate vector geometry for each location in a hexagonal AbstractGrid, up until some maximum
    x and y UTM coordinate. This can be used for visualization purposes.

    Inputs:
    grid: An AbstractGrid
    max_utm_x, max_utm_y: Maximum X and Y coordinate that the grid must cover, in UTM

    Output: An AbstractGridGDF (aka a GeoDataFrame, cols: ['grid_x', 'grid_y', 'geometry'[Polygon]])
    """
    min_x, min_y, max_x, max_y = (grid.origin_x, grid.origin_y, max_utm_x, max_utm_y)

    # Hexagon: A = 3/2 * sqrt(3) * R^2 where R is "size" (distance from center to any point)
    # Solving for R when A is `grid.area`
    side_length_m = np.sqrt(2 / (3 * np.sqrt(3)) * grid.area)

    # Number of rows and columns needed to completely cover the requested extent
    nx = np.ceil((max_x - min_x) / (np.sqrt(3) * side_length_m)).astype(int)
    ny = np.ceil((max_y - min_y) / (3 / 2 * side_length_m)).astype(int)

    # --- <Build the hexagons> ---

    # X and Y indices (these are *not* UTM coordinates)
    xs, ys = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    xs, ys = np.reshape(xs, (-1,)), np.reshape(ys, (-1,))

    # X and Y coordinates of hexagon centers (these are UTM coordinates).
    #
    # Note: mod() is used to account for the x-offset of odd-numbered columns
    cx = min_x + np.where(
        np.mod(ys, 2),
        side_length_m * (np.sqrt(3) * xs + np.sqrt(3) / 2),
        side_length_m * np.sqrt(3) * xs,
    )
    cy = min_y + side_length_m * 3 / 2 * ys

    # X and Y coordinates of a hexagon centered on (0,0)
    #
    # We need a 7th point to close the polygon - so there will be two points corresponding to angle
    # pi/6
    angles_for_points_on_a_pointy_topped_hexagon = 2 * np.pi * np.linspace(0, 1, 7) + np.pi / 6
    polygon_xs = side_length_m * np.cos(angles_for_points_on_a_pointy_topped_hexagon)
    polygon_ys = side_length_m * np.sin(angles_for_points_on_a_pointy_topped_hexagon)

    # X and Y coordinates for each hexagon.
    # Dims: hexagon, coordinate, x/y
    # Shape: nx*ny, 7, 2
    polygons_coords = np.stack(
        (
            # [7 repeats of the center-coordinates] + [the "7" points on the hexagon]
            np.reshape(np.repeat(cx, 7) + np.tile(polygon_xs, nx * ny), (-1, 7)),
            np.reshape(np.repeat(cy, 7) + np.tile(polygon_ys, nx * ny), (-1, 7)),
        ),
        axis=-1,
    )

    # --- </Build the hexagons> --

    return gpd.GeoDataFrame(
        {
            "grid_x": xs,
            "grid_y": ys,
            "geometry": polygons(polygons_coords),
        },
        crs=pyproj.CRS.from_epsg(grid.epsg_code),
    ).set_index(["grid_x", "grid_y"])


def _make_square_abstract_grid_gdf(
    grid: AbstractGrid,
    max_utm_x: int,
    max_utm_y: int,
) -> AbstractGridGDF:
    """
    Generate vector geometry for each location in a square AbstractGrid, up until some maximum
    x and y UTM coordinate. This can be used for visualization purposes.

    Inputs:
    grid: An AbstractGrid
    max_utm_x, max_utm_y: Maximum X and Y coordinate that the grid must cover, in UTM

    Output: An AbstractGridGDF (aka a GeoDataFrame, cols: ['grid_x', 'grid_y', 'geometry'[Polygon]])
    """
    min_x, min_y, max_x, max_y = (grid.origin_x, grid.origin_y, max_utm_x, max_utm_y)

    side_length_m = np.sqrt(grid.area)

    xs, ys = np.arange(min_x, max_x, side_length_m), np.arange(min_y, max_y, side_length_m)
    xx, yy = np.meshgrid(xs, ys)

    xind, yind = np.meshgrid(np.arange(xs.shape[0]), np.arange(ys.shape[0]))

    boxes = box(
        xmin=np.reshape(xx, -1),
        ymin=np.reshape(yy, -1),
        xmax=np.reshape(xx, -1) + side_length_m,
        ymax=np.reshape(yy, -1) + side_length_m,
    )

    return gpd.GeoDataFrame(
        {
            "grid_x": np.reshape(xind, -1),
            "grid_y": np.reshape(yind, -1),
            "geometry": boxes,
        },
        crs=pyproj.CRS.from_epsg(grid.epsg_code),
    ).set_index(["grid_x", "grid_y"])


def make_abstract_grid_gdf(
    grid: Union[AbstractGrid, BoundedAbstractGrid],
    max_utm_x: Optional[int] = None,
    max_utm_y: Optional[int] = None,
    tile_size: Optional[tuple[int, int]] = None,
) -> Union[AbstractGridGDF, GridIDAbstractGridGDF]:
    """
    Generate vector geometry for each location in a AbstractGrid, up until some maximum
    x and y UTM coordinate. This can be used for visualization purposes.

    Inputs:
    grid: An AbstractGrid

    max_utm_x, max_utm_y: Maximum X and Y coordinate that the grid must cover, in UTM.
        For `AbstractGrid`: You must pass this.
        For `BoundedAbstractGrid`: You may pass this. Otherwise, `extent_x` and `extent_y` will be
            used, and a GeoDataFrame of the entire grid will be returned.

    tile_size:
        BoundedAbstractGrid only, optional: generate a `grid_id` column using
        `coords.generate_tiled_grid_id_coord()`. This makes it possible to sort data using
        `grid_id`, which will improve the spatial locality. But, because the `grid_id` depends on
        the `tile_size`, it should not be used as the primary key of a grid-cell.

    Output:
        An AbstractGridGDF (aka a GeoDataFrame, cols: ['grid_x', 'grid_y', 'geometry'[Polygon]])
        Or, if `tile_size` is not None:
            A GridIDAbstractGridGDF (with extra col 'grid_id')
    """
    if not isinstance(grid, BoundedAbstractGrid) and (max_utm_x is None or max_utm_y is None):
        raise ValueError("Must pass `max_utm_x`, `max_utm_y` if not using a BoundedAbstractGrid")

    # If the user passed a BoundedAbstractGrid, but they didn't specify what the maximum
    # grid-coordinates are to be, give them all grid-cells in the bounding box specified by
    # [origin_x, origin_y, extent_x, extent_y].
    if isinstance(grid, BoundedAbstractGrid) and (max_utm_x is None or max_utm_y is None):
        max_utm_x, max_utm_y = grid.extent_x, grid.extent_y

    if grid.geometry_type == GridGeometry.HEXAGON:
        result = _make_hexagonal_abstract_grid_gdf(grid, max_utm_x, max_utm_y)
    elif grid.geometry_type == GridGeometry.SQUARE:
        result = _make_square_abstract_grid_gdf(grid, max_utm_x, max_utm_y)
    else:
        raise ValueError(f"Unsupported grid geometry type {grid.geometry_type}")

    if isinstance(grid, BoundedAbstractGrid) and tile_size is not None:
        grid_ids = generate_tiled_grid_id_coord(
            result["grid_x"], result["grid_y"], tile_size=tile_size, grid=grid
        )

        # Build the new index, using "grid_id" alone.
        result["grid_id"] = pd.Series(grid_ids)

    return result


def make_gdf_from_coords(coords: np.ndarray, grid: AbstractGrid) -> AbstractGridGDF:
    """Create an geodataframe from the coordinates of the grid cells. The geodataframe
    will contain the geometries of the tiles at the given coordinates. The geodataframe
    will have a multiindex with the grid coordinates in the order of the input array.

    Args:
        coords (np.ndarray): (n,2) AbstractGrid-coordinates of the grid-cells
        grid (AbstractGrid): AbstractGrid that the coords are defined on

    Returns:
        AbstractGridGDF (aka a GeoDataFrame, cols: ['grid_x', 'grid_y', 'geometry'[Polygon]])
    """
    if isinstance(grid, BoundedAbstractGrid):
        raise AssertionError(
            "Cannot create a GeoDataFrame from coordinates for a BoundedAbstractGrid"
        )
    max_grid_x, max_grid_y = np.max(coords[:, 0]), np.max(coords[:, 1])
    max_utm_x, max_utm_y = convert_grid_to_utm(
        np.array([max_grid_x + 1, max_grid_y + 1]), grid=grid
    )
    gdf = make_abstract_grid_gdf(
        grid, max_utm_x=max_utm_x, max_utm_y=max_utm_y, tile_size=grid.area
    )
    index = pd.MultiIndex.from_arrays([coords[:, 0], coords[:, 1]], names=["grid_x", "grid_y"])
    return gdf.loc[index]
