from typing import Union

import numpy as np

#from msat_level4_gim.errors import ShapeMismatchException
from errors import ShapeMismatchException

def negative_reflection(
    full_data: np.ndarray,
    averaging_kernel: Union[np.ndarray, float],
    background: Union[np.ndarray, float],
) -> np.ndarray:
    """
    Given a dataset, an averaging kernel, and a proposed background, return the reflected negative
    noise distribution. See README for more details.

    Inputs:
        full_data: shape S in n-dimensions. Must be in same units as background.
        averaging_kernel: shape S or single value.
        background: Single value or shape that broadcasts to S. Must be in same units as full_data.
            Note: This is equivalent to *c* in the documentation.
    Returns:
        1-d array of same size as full_data, containing the distribution.
            Note: Will be nan-free if the inputs are.
    """
    if isinstance(averaging_kernel, np.ndarray) and averaging_kernel.shape != full_data.shape:
        raise ShapeMismatchException("Averaging Kernel and Data must have same shape.")
    if isinstance(background, np.ndarray) and background.shape != full_data.shape:
        raise ShapeMismatchException("Background and Data must have same shape.")
    data_negative = (full_data - background * averaging_kernel).flatten()
    data_0 = data_negative[data_negative == 0]
    data_negative = data_negative[data_negative < 0]
    data_positive = -1 * data_negative
    combined = np.concatenate([data_negative, data_0, data_positive])
    return combined
