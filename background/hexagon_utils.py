import numpy as np


def cube_round(cube: np.ndarray) -> np.ndarray:
    """
    This is a vectorized implementation of this function:
    https://www.redblobgames.com/grids/hexagons/#rounding

    n * 3(q,r,s) -> n * 3(q,r,s)
    """
    q_orig, r_orig, s_orig = cube[:, 0], cube[:, 1], cube[:, 2]
    q, r, s = np.round(cube[:, 0]), np.round(cube[:, 1]), np.round(cube[:, 2])

    d_q = np.abs(q - q_orig)
    d_r = np.abs(r - r_orig)
    d_s = np.abs(s - s_orig)

    cond1 = np.logical_and(d_q > d_r, d_q > d_s)

    q = np.where(cond1, -1 * r - s, q)

    cond2 = d_r > d_s

    r = np.where(
        np.logical_and(np.logical_not(cond1), cond2),
        -1 * q - s,
        r
    )

    s = np.where(
        np.logical_and(np.logical_not(cond1), np.logical_not(cond2)),
        -1 * q - r,
        s
    )

    return np.stack((q, r, s), axis=-1).astype(int)


def axial_to_cube(axial: np.ndarray) -> np.ndarray:
    """
    This is a vectorized implementation of this function:
    https://www.redblobgames.com/grids/hexagons/#conversions-axial

    n * 2(q,r) -> n * 3(q,r,s)
    """
    q, r = axial[:, 0], axial[:, 1]
    s = -1 * q - r
    return np.stack((q, r, s), axis=-1)


def cube_to_axial(cube: np.ndarray) -> np.ndarray:
    """
    This is a vectorized implementation of this function:
    https://www.redblobgames.com/grids/hexagons/#conversions-axial

    n * 3(q,r,s) -> n * 2(q,r)
    """
    return cube[:, 0:2]


def axial_round(axial: np.ndarray) -> np.ndarray:
    """
    This is a vectorized implementation of this function:
    https://www.redblobgames.com/grids/hexagons/#rounding

    n * 2(q,r) -> n * 2(q,r)
    """
    return cube_to_axial(cube_round(axial_to_cube(axial)))


def axial_to_oddr(axial: np.ndarray) -> np.ndarray:
    """
    This is a vectorized implementation of this function:
    https://www.redblobgames.com/grids/hexagons/#conversions-offset

    n * 2(q,r) -> n * (r,c)
    """
    q, r = axial[:, 0], axial[:, 1]
    row = q + (r - np.mod(r, 2)) / 2
    col = r
    return np.stack((row, col), axis=-1).astype(int)


def cartesian_to_pointy_hex(points: np.ndarray, size: float) -> np.ndarray:
    """
    This is a vectorized implementation of this function:
    https://www.redblobgames.com/grids/hexagons/#pixel-to-hex

    n * 2(x,y) -> n * 2(q,r)
    """
    x, y = points[:, 0], points[:, 1]
    q = (np.sqrt(3) / 3 * x - 1 / 3 * y) / size
    r = (2 / 3 * y) / size
    return axial_round(np.stack((q, r), axis=-1))


def cartesian_to_pointy_hex_offset_coord(points: np.ndarray, size: float = 1) -> np.ndarray:
    """
    Given an array of points in a cartesian coordinate system, calculates which hexagon each point
    falls into, given a pointy-topped "odd-r" hexagonal grid with a hexagon centered on (0,0).

    For a definition of "odd-r", see:
    https://www.redblobgames.com/grids/hexagons/#coordinates-offset

    This is a vectorized implementation of this function:
    https://www.redblobgames.com/grids/hexagons/#pixel-to-hex

    n * 2(x,y) -> n * 2(r,c)
    """
    return axial_to_oddr(cartesian_to_pointy_hex(points, size))


def hexagon_size_for_area(area: float) -> float:
    """
    Return hexagon "size" for a given area (size is distance from center to a point on the edge)
    See https://www.redblobgames.com/grids/hexagons/#basics for more info.
    """
    return np.sqrt(2 / (3 * np.sqrt(3)) * area)
