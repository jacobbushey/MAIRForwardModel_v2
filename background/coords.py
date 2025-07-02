from typing import Union

import numpy as np
import rioxarray  # noqa: F401
from xarray import DataArray

#from msat_level4_gim.abstract_grid.data_types import (
#    AbstractGrid,
#    BoundedAbstractGrid,
#    GridGeometry,
#)
#from msat_level4_gim.abstract_grid.hexagon_utils import (
#    cartesian_to_pointy_hex_offset_coord,
#    hexagon_size_for_area,
#)
from abstract_grid_data_types import (
    AbstractGrid,
    BoundedAbstractGrid,
    GridGeometry,
)
from hexagon_utils import (
    cartesian_to_pointy_hex_offset_coord,
    hexagon_size_for_area,
)


def _crs_x_to_square_grid_x_coords(crs_x_coords: np.ndarray, grid: AbstractGrid) -> np.ndarray:
    """
    Computes the following:
      "What are the grid-x-coordinates in a square grid `grid` for a collection of x-coordinates
       that reference the same CRS as `grid`?"

    Inputs:
    crs_x_coords: N-length array of x coordinates, must be in the same CRS as `grid`
    grid: An AbstractGrid

    Outputs:
    N-length array of grid-x coordinates
    """
    if grid.geometry_type != GridGeometry.SQUARE:
        raise ValueError(f"grid geometry must be square, but was: {grid.geometry_type}")

    grid_x_res = np.sqrt(grid.area)
    dx = crs_x_coords - grid.origin_x
    return np.floor(dx / grid_x_res).astype(int)


def _crs_y_to_square_grid_y_coords(crs_y_coords: np.ndarray, grid: AbstractGrid) -> np.ndarray:
    """
    Computes the following:
      "What are the grid-y-coordinates in a square grid `grid` for a collection of y-coordinates
       that reference the same CRS as `grid`?"

    Inputs:
    crs_y_coords: N-length array of y coordinates, must be in the same CRS as `grid`
    grid: An AbstractGrid

    Outputs:
    N-length array of grid-y coordinates
    """
    if grid.geometry_type != GridGeometry.SQUARE:
        raise ValueError(f"grid geometry must be square, but was: {grid.geometry_type}")

    grid_y_res = np.sqrt(grid.area)
    dy = crs_y_coords - grid.origin_y
    return np.floor(dy / grid_y_res).astype(int)


def _crs_xy_to_pointy_hexagon_grid_coords(
    crs_xy_coords: np.ndarray, grid: AbstractGrid
) -> np.ndarray:
    """
    Computes the following:
      "What are the grid-xy-coordinates in a hexagonal grid `grid` for a collection of
       xy-coordinates that reference the same CRS as `grid`?"

    Inputs:
    crs_xy_coords: N * 2[x,y] array of x-y coordinates, must be in the same CRS as `grid`
    grid: An AbstractGrid

    Outputs:
    N * 2[grid_x,grid_y] array
    """
    if grid.geometry_type != GridGeometry.HEXAGON:
        raise ValueError(f"grid geometry must be hexagon, but was: {grid.geometry_type}")

    crs_x_coords, crs_y_coords = crs_xy_coords[:, 0], crs_xy_coords[:, 1]
    return cartesian_to_pointy_hex_offset_coord(
        points=np.stack((crs_x_coords - grid.origin_x, crs_y_coords - grid.origin_y), axis=-1),
        size=hexagon_size_for_area(grid.area),
    )


def convert_utm_to_grid(utm_coords: np.ndarray, grid: AbstractGrid) -> np.ndarray:
    """
    Create `grid_x' and `grid_y` coordinates for a UTM-projected set of coordiantes on a given
    AbstractGrid.

    Inputs:
    utm_coords: Shape (..., 2) coordiantes in grid.epsg_code projected UTM coordinates.
    grid: AbstractGrid for which to calculate indexes.

    Output: Same shape as `utm_coords` but integer dtype, represents what grid cell each of the
            `utm_coords` falls in on the `grid`.
    """
    if grid.geometry_type == GridGeometry.HEXAGON:
        # Convert xy-CRS-coordinates to grid-x, grid-y coordinates in a hexagonal grid
        # Shape: ... * 2[grid_x, grid_y]
        grid_coords = _crs_xy_to_pointy_hexagon_grid_coords(
            crs_xy_coords=utm_coords.reshape((-1, 2)), grid=grid
        )

        return grid_coords.reshape(utm_coords.shape)

    elif grid.geometry_type == GridGeometry.SQUARE:
        xs, ys = utm_coords[..., 0], utm_coords[..., 1]
        grid_x = _crs_x_to_square_grid_x_coords(xs, grid)
        grid_y = _crs_y_to_square_grid_y_coords(ys, grid)

        return np.stack((grid_x, grid_y), axis=-1)

    else:
        raise ValueError(f"Unsupported geometry type {grid.geometry_type}")


def generate_grid_xy_coords(
    crs_x_coords: DataArray, crs_y_coords: DataArray, grid: AbstractGrid
) -> tuple[
    # xarray coordinates are assigned by passing tuples that look like this:
    #   > ("orchard_location", ["apple" "orange" "mango"])
    #
    # or like this, when the coordinate you're assigning is a function of more than one dimension:
    #   > (("orchard_row", "orchard_column"), [["apple" "apple"] ["orange" "mango"]])
    tuple[Union[str, tuple[str, str]], np.ndarray],  # <- the grid_x coordinate
    tuple[Union[str, tuple[str, str]], np.ndarray],  # <- the grid_y coordinate
]:
    """
    Create `grid_x' and `grid_y` coordinates for a UTM-projected level3 product (segment or mosaic).
    If a square grid, `grid_x` and `grid_y` will be a function of `x` and `y`, respectively.
    If a hex grid, `grid_x` and `grid_y` will *both* be a function of `x` *and* `y`.

    Inputs:
    x_coords: "x" coords DataArray from a level3 product. Must have already been projected into the
        same UTM zone as is specified by `grid.epsg_code`.
    y_coords: "y" coords DataArray from a level3 product. Must have already been projected into the
        same UTM zone as is specified by `grid.epsg_code`.
    grid: AbstractGrid

    Outputs:
    Coordinate tuples for grid_x and grid_x. These can be passed to `xr.assign_coords()`, or
    assigned to a Dataset using the `ds.coords["coord_name"] = coords` syntax.
    """
    crs_x_coords = crs_x_coords.to_numpy()  # Convert DataArray -> ndarray
    crs_y_coords = crs_y_coords.to_numpy()
    # Convert to grid coordinates
    grid_coords = convert_utm_to_grid(np.stack((crs_x_coords, crs_y_coords), axis=-1), grid)
    grid_x, grid_y = grid_coords[..., 0], grid_coords[..., 1]

    # If BoundedAbstractGrid, we check here to make sure the grid_x, grid_y values being produced
    # are valid (within the specified extent of the BoundedAbstractGrid).
    #
    # A note on validity: A CRS coordinate can be considered valid even if it lies to the east of
    # `extent_x` and/or to the north of `extent_y`, iff its enclosing BoundedAbstractGrid grid-cell
    # is a member of the easternmost and/or northernmost grid_x/grid_y column.
    #
    # => In other words, geospatial extent is evaluated at the grid-cell level, not at the CRS
    #    level.
    if isinstance(grid, BoundedAbstractGrid):
        extent_coords = convert_utm_to_grid(np.array([[grid.extent_x, grid.extent_y]]), grid)
        max_grid_x, max_grid_y = extent_coords[0, 0], extent_coords[0, 1]

        if grid.geometry_type == GridGeometry.HEXAGON:
            # The maximum grid values of UTM coords on the borders can be higher
            # than the grid coords of the UTM extent values in a hex grid
            max_grid_x = max_grid_x + 1
            max_grid_y = max_grid_y + 1

        if np.any(grid_x > max_grid_x):
            offending_grid_cells = np.where(grid_x > max_grid_x)
            raise ValueError(
                f"The max grid-x for this grid is {max_grid_x}, but the following grid-cells had"
                f"higher grid-x values:\n{offending_grid_cells}. The value of `extent_x` was "
                f"{grid.extent_x}."
            )

        if np.any(grid_y > max_grid_y):
            offending_grid_cells = np.where(grid_y > max_grid_y)
            raise ValueError(
                f"The max grid-y for this grid is {max_grid_y}, but the following grid-cells had"
                f"higher grid-y values:\n{offending_grid_cells}. The value of `extent_y` was "
                f"{grid.extent_y}."
            )

    # Assign coordinates to Dataset, mapping combinations of the "x" and "y" dimensions "grid_x"
    # and "grid_y"
    if grid.geometry_type == GridGeometry.HEXAGON:
        grid_x_xarray_coords = (("x", "y"), grid_x)
        grid_y_xarray_coords = (("x", "y"), grid_y)

        return grid_x_xarray_coords, grid_y_xarray_coords

    elif grid.geometry_type == GridGeometry.SQUARE:
        grid_x_xarray_coords = ("x", grid_x)
        grid_y_xarray_coords = ("y", grid_y)

        return grid_x_xarray_coords, grid_y_xarray_coords

    else:
        raise ValueError(f"Unsupported geometry type {grid.geometry_type}")


def convert_grid_to_utm(
    grid_coords: np.ndarray, grid: AbstractGrid
) -> tuple[np.ndarray, np.ndarray]:
    """
    Given a set of indices for a given `grid`, return the centers of the corresponding grid cells.

    Inputs:
    grid_coords: Shape (..., 2) x/y grid indices.
    grid: AbstractGrid on which the coordinates are defined.

    Output: Tuple of UTM x's and y's. Each of the two arrays has the same shape as `grid_coords`,
        the coordinates of those grid locations in UTM defined by `grid.epsg_code`.
    """
    # X and Y indices (these are *not* UTM coordinates)
    xs, ys = grid_coords[..., 0], grid_coords[..., 1]
    if grid.geometry_type == GridGeometry.HEXAGON:
        min_x, min_y = grid.origin_x, grid.origin_y

        # Hexagon: A = 3/2 * sqrt(3) * R^2 where R is "size" (distance from center to any point)
        # Solving for R when A is `grid.area`
        size = np.sqrt(2 / (3 * np.sqrt(3)) * grid.area)

        # X and Y coordinates of hexagon centers (these are UTM coordinates).
        #
        # Note: mod() is used to account for the x-offset of odd-numbered columns
        cx = min_x + np.where(
            np.mod(ys, 2), size * (np.sqrt(3) * xs + np.sqrt(3) / 2), size * np.sqrt(3) * xs
        )
        cy = min_y + size * 3 / 2 * ys
        return cx, cy

    elif grid.geometry_type == GridGeometry.SQUARE:
        square_dimension = np.sqrt(grid.area)

        utm_x_center = grid.origin_x + square_dimension * xs + 0.5 * square_dimension
        utm_y_center = grid.origin_y + square_dimension * ys + 0.5 * square_dimension

        return utm_x_center, utm_y_center

    else:
        raise ValueError(f"Unsupported geometry type {grid.geometry_type}")


def generate_utm_xy_coords(
    grid_x_coords: DataArray, grid_y_coords: DataArray, grid: AbstractGrid
) -> tuple[
    # xarray coordinates are assigned by passing tuples that look like this:
    #   > ("orchard_location", ["apple" "orange" "mango"])
    #
    # or like this, when the coordinate you're assigning is a function of more than one dimension:
    #   > (("orchard_row", "orchard_column"), [["apple" "apple"] ["orange" "mango"]])
    tuple[Union[str, tuple[str, str]], np.ndarray],  # <- the grid_x coordinate
    tuple[Union[str, tuple[str, str]], np.ndarray],  # <- the grid_y coordinate
]:
    """
    For a given list of grid_x and grid_y coords, return the UTM x-and-y-coordinates of the centers
    of the grid-cells in the cartesian product of the x-and-y-coordinates.

    Inputs:
    grid_x_coords: An integer DataArray
    grid_y_coords: An integer DataArray

    Outputs:
    Coordinate tuples for utm_x and utm_y. These can be passed to `xr.assign_coords()`, or
    assigned to a Dataset using the `ds.coords["coord_name"] = coords` syntax.
    """
    grid_x_coords = grid_x_coords.to_numpy()
    grid_y_coords = grid_y_coords.to_numpy()

    utm_x_center, utm_y_center = convert_grid_to_utm(
        np.stack((grid_x_coords, grid_y_coords), axis=-1), grid
    )

    if grid.geometry_type == GridGeometry.HEXAGON:
        utm_x_xarray_coords = (("grid_x", "grid_y"), utm_x_center)
        utm_y_xarray_coords = (("grid_x", "grid_y"), utm_y_center)
        return utm_x_xarray_coords, utm_y_xarray_coords

    elif grid.geometry_type == GridGeometry.SQUARE:
        utm_x_xarray_coords = ("grid_x", utm_x_center)
        utm_y_xarray_coords = ("grid_y", utm_y_center)
        return utm_x_xarray_coords, utm_y_xarray_coords

    else:
        raise ValueError(f"Unsupported geometry type {grid.geometry_type}")


def generate_tiled_grid_id_coord(
    grid_xs: np.ndarray,
    grid_ys: np.ndarray,
    tile_size: tuple[int, int],
    grid: BoundedAbstractGrid,
) -> np.ndarray:
    """
    For a given list of grid_x and grid_y coords, return unique integer IDs (`grid_id`s) for each
    grid-cell. A tiling approach is used to assign IDs, which increases the average spatial locality
    of `grid_id`s with similar values.

    Example with `tile_size == (3, 2)`:

      ┌─────────────────────────────────────┐      ┌─────────────────────────────────────┐
      │  ┌───┐  ┌───┐  ┌───┐  ┌───┐  ┌─22┐  │      │  ┌─────────────────┐  ┌─────────────│
     3┼  │15 │  │16 │  │17 │  │21 │ ▲extent │     3┼  │                 │  │     ▲extent │
      │  └───┘  └───┘  └───┘  └───┘  └───┘  │      │  │    tile 2       │  │     tile 3  │
      │  ┌───┐  ┌───┐  ┌───┐  ┌───┐  ┌───┐  │      │  │                 │  │             │
     2┼  │12 │  │13 │  │14 │  │18 │  │19 │  │     2┼  │                 │  │             │
      │  └───┘  └───┘  └───┘  └───┘  └───┘  │      │  └─────────────────┘  └─────────────│
      │  ┌───┐  ┌───┐  ┌───┐  ┌───┐  ┌───┐  │      │  ┌─────────────────┐  ┌─────────────│
     1┼  │ 3 │  │ 4 │  │ 5 │  │ 9 │  │ 10│  │     1┼  │                 │  │             │
      │  └───┘  └───┘  └───┘  └───┘  └───┘  │      │  │    tile 0       │  │      tile 1 │
      │  ┌───┐  ┌───┐  ┌───┐  ┌───┐  ┌───┐  │      │  │                 │  │             │
     0┼  │ 0 │  │ 1 │  │ 2 │  │ 6 │  │ 7 │  │     0┼  │                 │  │             │
      │  └───┘  └───┘  └───┘  └───┘  └───┘  │      │  └─────────────────┘  └─────────────│
      └────┼──────┼──────┼──────┼──────┼────┘      └────┼──────┼──────┼──────┼──────┼────┘
           0      1      2      3      4                0      1      2      3      4

    ```python
    grid = BoundedAbstractGrid(...)

    generate_tiled_grid_id_coord(
        np.array[[1, 1, 1]],
        np.array([[1, 2, 3]]),
        tile_size=(3, 2),
        grid=grid,
    )  # => returns np.array([[4, 13, 16]])
    ```

    How `grid_id`s are assigned, starting with grid-cell (0, 0):

    1. Within a tile, increase while going east until you reach the end of the tile.
       While there is still a row to the north, go to the westernmost cell in the row, repeat.

    2. After completing a tile, travel east to the next tile and repeat step 1, unless the
       current tile is the easternmost tile in that column of tiles. In that case, while there
       is still a row of tiles to the north, go to the westernmost tile of that next row, and
       repeat step 1.

    The easternmost column of tiles is the tile-column that contains the UTM coordinate
    (extent_x, extent_y).

    Inputs:
    grid_xs: An 2D integer ndarray, shape (n, m) or (n,)
    grid_ys: An 2D integer ndarray, shape (n, m) or (n,)
    tile_size: (int, int) tuple
        Size of a tile. This determines how `grid_id` is generated for each grid-cell.
        Must be positive integers.

        units: [grid_x, grid_y] (number of grid-cells in the x- and y- directions).

    grid: BoundedAbstractGrid instance.

    Outputs:
    ndarray of `grid_id`s, same shape as `grid_xs`, `grid_ys`
    """
    if not isinstance(grid, BoundedAbstractGrid):
        raise ValueError(
            "`grid` was not a BoundedAbstractGrid. AbstractGrids can't have `grid_id`s."
        )

    # Check for 2-tuple
    if len(tile_size) != 2:
        raise ValueError(f"`tile_size` must be a 2-tuple, but a {len(tile_size)}-tuple was passed")

    # Check for int-tuple
    if type(tile_size[0]) is not int or type(tile_size[1]) is not int:
        raise ValueError("tile_size must be an int-tuple, denominated in grid-cells")

    if tile_size[0] < 1 or tile_size[1] < 1:
        raise ValueError(f"Tile XY dimensions must be positive, but were: {tile_size}")

    tile_size_x, tile_size_y = tile_size

    # tile-x-index of the furthest-east tile-column. `extent` is contained within this column
    extent_grid_coords = convert_utm_to_grid(np.array([[grid.extent_x, grid.extent_y]]), grid)

    max_tile_x = np.floor_divide(extent_grid_coords[0, 0], tile_size_x)

    # x- and y- indices of the tile enclosing each AbstractGrid cell
    (tile_x, rem_x), (tile_y, rem_y) = (
        np.divmod(grid_xs, tile_size_x),
        np.divmod(grid_ys, tile_size_y),
    )

    tile_local_idx = rem_y * tile_size_x + rem_x

    # The tile id of each AbstractGrid cell
    tile_id = tile_y * (max_tile_x + 1) + tile_x

    # The "grid_id" of the "first" cell in each tile.
    # The "first" cell is the cell in the southwestern corner of the tile.
    tile_start = (tile_size_x * tile_size_y) * tile_id

    grid_id = tile_start + tile_local_idx

    # "grid_id" is just the "grid_id" of the first tile in your cell, plus your "tile-local index"
    # within that cell.
    return grid_id
