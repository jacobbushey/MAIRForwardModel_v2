# NOTE: THIS IS A DRAFT METHOD THAT DOES NOT WORK (YET)
# ZSIGMA TEST IS THE PROPOSEd METHOD
from numba import njit

import numpy as np

#from msat_level4_gim.background.data_types import (
#    BackgroundCalculationAnalytics,
#    BackgroundCalculationOptions,
#    BackgroundCalculationResults,
#)
#from msat_level4_gim.background.negative_reflection import negative_reflection
#from msat_level4_gim.background.utility import convolutional_test

from data_types import (
    BackgroundCalculationAnalytics,
    BackgroundCalculationOptions,
    BackgroundCalculationResults,
)
from negative_reflection import negative_reflection
from utility import convolutional_test

@njit
def _local_deviation_in_std(data, mean_ak, bg, correlation_factor=np.nan) -> np.ndarray:
    """Returns the mean local deviation of the measuremnts to the background in units of the
    standard deviation of the region. Since it is the mean deviation, this method requires
    and estimate of the local information centent. Since this is unknown at the moment, this
    method is **NOT OPERATIONAL**.

    Inputs:
        data: xch4 offset from prior. Any shape
        mean_ak: The column averaged averaging kernel. Same shape as l3_data.
        bg: The background concentration to test. Same shape as l3_data.
        correlation_factor: The autocorrelation factor between [0,1]. Unknown at the moment.
    Output:
        [deviation_in_sigma, local_std, NAN]. Three returned values are required.
    """
    if data.size < 20:
        return np.array([np.nan, np.nan, np.nan])

    # robust against NANs
    data = data[np.isfinite(data)]
    N = data.sum()
    local_std = np.std(data)
    deviation_in_sigma = np.mean(data - bg) / local_std / np.sqrt(N * correlation_factor)
    return np.array([deviation_in_sigma, local_std, np.nan])


def global_std_fit(
    l3_data: np.ndarray, mean_ak: np.ndarray, options: BackgroundCalculationOptions
) -> tuple[BackgroundCalculationResults, BackgroundCalculationAnalytics]:
    raise NotImplementedError("Global std test not implemented.")


def local_std_fit(
    l3_data: np.ndarray, mean_ak: np.ndarray, options: BackgroundCalculationOptions
) -> tuple[BackgroundCalculationResults, BackgroundCalculationAnalytics]:
    """
    Skeleton for a local standard deviation test for the background.
    Computes the standard deviation of a region of the data, and then
    uses the negative reflected noise, normalized by the region standard deviation,
    to determine the background concentration.
    That works up to an unknown factor for the effective sample size within the region.
    Since this is unknown, this tes is ***NOT OPERATIONAL*** at the moment.

    This correction factor could be (pixel-varying) sqrt(total sample weight), plus a constant
    derived from the L2 sampling autocorrelation properties plus some other pixel-varying constan
    for atmospheric autocorrelation.

    Inputs:
        l3_data: The raw L3 xch4 offset from prior. Any shape
        mean_ak: The column averaged averaging kernel. Same shape as l3_data.
        start: Optional starting point for search
        stop: Optional end point for search
        region_size: How large a region to do estimation on. Manhattan distance radius pixels.
        region_step: How far between regions? Manhattan distance pixels.
    Returns:
        BackgroundCalculationResults: The best background concentration and standard deviation
        BackgroundCalculationAnalytics: Analytics
    """
    # Handle masked arrays in order to use numba.
    if isinstance(l3_data, np.ma.MaskedArray):
        l3_data = l3_data.filled(np.nan)
    if isinstance(mean_ak, np.ma.MaskedArray):
        mean_ak = mean_ak.filled(np.nan)

    bgs_to_test, results = convolutional_test(
        l3_data,
        mean_ak,
        _local_deviation_in_std,
        options.search_start,
        options.search_stop,
        options.search_steps,
        options.region_size,
        options.region_spacing,
    )
    # results has shape (x, y, bg, 3)
    x, y, bg, r = results.shape
    assert bg == bgs_to_test.size
    # flatten xy-coords
    results = results.reshape(x * y, bg, r)  # shape (x*y, bg, 3)

    deviation_in_sigma, local_std = (
        results[..., 0],
        results[..., 1],
    )  # shape (x*y, bg)

    res = np.full(bgs_to_test.shape, np.nan)
    for i, bg in enumerate(bgs_to_test):
        bg_dev = deviation_in_sigma[..., i]
        bg_std = local_std[..., i]

        # mask out very large std-areas (plumes, topo-edges, cloud-effects, stripes, you name it)
        # this mask proved reasonable when I investigated the real data
        mask = bg_std < 3 * np.median(bg_std)

        # find negative deviations of valid regions
        bg_dev = bg_dev[np.logical_and(mask, bg_dev <= 0)]

        if bg_dev.size > 10:
            # reflected negative noise
            res[i] = np.nanstd(negative_reflection(bg_dev, 0, 0))
        else:
            # no valid regions
            res[i] = np.nan

    # Find the best background, which is the one with the std closest to 1
    best_bg = bgs_to_test[np.nanargmin(np.abs(res - 1))]
    bg_std = 0

    return BackgroundCalculationResults(best_bg, bg_std), BackgroundCalculationAnalytics(
        background_candidates=bgs_to_test.tolist(),
        background_candidate_calculated_standard_deviation=res.tolist(),
    )
