from numba import jit
import numpy as np


@jit(nopython=True)  # Set "nopython" mode for best performance, equivalent to @njit
def distance_matrix(x_coords: np.ndarray, y_coords: np.ndarray) -> np.ndarray:
    """
    Calculates the Euclidian distances between all pairs of points [x_i, y_i], [x_j, y_j] for
    all i in 0...n-1 and all j in 0...n-1.

    Inputs:
    x_coords: x-coordinates of the 2D points
    y_coords: y-coordinates of the 2D points

    Output:
    np.ndarray, (n, n), symmetric, 0s on the diagonal
    """
    if x_coords.shape != y_coords.shape:
        raise ValueError(
            "x_coords and y_coords must have same shape, but were "
            f"{x_coords.shape}, {y_coords.shape}"
        )

    # Allocate the distance matrix
    distances = np.zeros((x_coords.shape[0], x_coords.shape[0]))

    # Calculate the entries for each column
    # (only calculate the entries that are on or above the diagonal)
    for i in range(x_coords.shape[0]):
        # Euclidian distance
        distances[: (i + 1), i] = np.sqrt(
            np.power(x_coords[i] - x_coords[: (i + 1)], 2)
            + np.power(y_coords[i] - y_coords[: (i + 1)], 2)
        )

    # add the transpose, forming a symmetric distance matrix
    return distances + distances.T
