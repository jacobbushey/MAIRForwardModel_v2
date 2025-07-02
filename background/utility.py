from typing import Callable, Optional

import numpy as np


def whole_scene_test(
    l3_data: np.ndarray,
    mean_ak: np.ndarray,
    test_function: Callable[[np.ndarray, np.ndarray, float], np.ndarray],
    start: Optional[float] = None,
    stop: Optional[float] = None,
    steps: Optional[int] = 200,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Calls the test function on all possible backgrounds and returns the results.
    Inputs:
        l3_data: Shape S raw L3 xch4 values, likely observed XCH4 minus prior. Can contain NaN.
        mean_ak: Shape S column averaged averaging kernel.
        test_function: Function that will be applied to scan over possible backgrounds.
            It should take l3_data, mean_ak, and a float candidate background. NaNs in arrays okay.
            It should return array size 3 of statistics for that combination. NaNs okay.
                We expect the return to be a statistic, p_value, and calculated std, but that is
                not strictly necessary.
        start: Predetermined start of background search. Defaults to 25% percentile of l3_data.
        stop: Predetermined stop of background search. Defaults to 60th percentile of l3_data.
        steps: The number of linearly spaced backgrounds to try.
    Returns:
        Return two np.ndarrays. The first is a size {steps} list of backgrounds.
        The second is of  size {steps}*3, and is the output from your test_function.
    """
    # Assume the correct background is between the 25% and 60%. We actually expect it to be <50%.
    minimum_guess = start if start else np.nanpercentile(l3_data, 25)
    maximum_guess = stop if stop else np.nanpercentile(l3_data, 60)

    bgs_to_test = np.linspace(minimum_guess, maximum_guess, steps)
    results = np.full((len(bgs_to_test), 3), np.nan)

    for ib, bg in enumerate(bgs_to_test):
        results[ib] = test_function(l3_data, mean_ak, bg)

    return bgs_to_test, np.stack(results, axis=0)


def convolutional_test(
    l3_data: np.ndarray,
    mean_ak: np.ndarray,
    test_function: Callable[[np.ndarray, np.ndarray, float], np.ndarray],
    start: Optional[float] = None,
    stop: Optional[float] = None,
    steps: Optional[int] = 50,
    region_size: Optional[int] = 50,
    region_step: Optional[int] = 25,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Calls the test function on convolved windows of the data for all possible backgrounds and
    returns the results.
    Inputs:
        l3_data: Shape S (2D) raw L3 xch4 values, likely observed XCH4 minus prior. Can contain NaN.
        mean_ak: Shape S (2D) column averaged averaging kernel.
        test_function: Function that will be applied to scan over possible backgrounds.
            It should take l3_data, mean_ak, and a float candidate background. NaNs in arrays okay.
            It should return array size 3 of statistics for that combination. NaNs okay.
                We expect the return to be a statistic, p_value, and calculated std, but that is
                not strictly necessary.
        start: Predetermined start of background search. Defaults to 10% percentile of l3_data.
        stop: Predetermined stop of background search. Defaults to 60th percentile of l3_data.
        steps: The number of linearly spaced backgrounds to try.
        region_size: How large is a convolutional region? Manhattan distance radius.
        region_step: How many pixels between region centers?
    Returns:
        Return two np.ndarrays. The first is a size {steps} list of backgrounds.
        The second is of size {NX, Ny, steps}*3, and is the output from your test_function.
            NX and NY refer to the number of convolution windows in the x and y dimensions of S.
    """
    # Assume the correct background is between the 10% and 60%. We actually expect it to be <50%.
    # l3_data = np.copy(l3_data)
    minimum_guess = start if start else np.nanpercentile(l3_data, 10)
    maximum_guess = stop if stop else np.nanpercentile(l3_data, 60)

    bgs_to_test = np.linspace(minimum_guess, maximum_guess, steps)

    x_max, y_max = l3_data.shape
    num_x, num_y = (x_max + region_step) // region_step, (y_max + region_step) // region_step

    # Shape is (x, y, bg, 3). 3 is an arbitrary size chosen to simplify this function.
    results = np.full((num_x, num_y, len(bgs_to_test), 3), np.nan)

    # For some efficiency, we take a subregion first, then do all the analysis we need on it.
    offset = region_size // 2
    for ix, x in enumerate(range(offset, x_max, region_step)):
        for iy, y in enumerate(range(offset, y_max, region_step)):
            # Take the subregion in xch4 and ak. Be careful not to go negative on indices.
            x_start, y_start = max(0, x - region_size), max(0, y - region_size)
            region = l3_data[x_start : x + region_size, y_start : y + region_size]
            region_ak = mean_ak[x_start : x + region_size, y_start : y + region_size]
            # Now loop over all backgrounds
            for ib, bg in enumerate(bgs_to_test):
                results[ix, iy, ib, :] = test_function(region, region_ak, bg)

    return bgs_to_test, results
