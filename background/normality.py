import warnings

import numpy as np
from scipy import stats

#from msat_level4_gim.background.data_types import (
#    BackgroundCalculationAnalytics,
#    BackgroundCalculationOptions,
#    BackgroundCalculationResults,
#)
#from msat_level4_gim.background.negative_reflection import negative_reflection
#from msat_level4_gim.background.utility import convolutional_test, whole_scene_test

from data_types import (
    BackgroundCalculationAnalytics,
    BackgroundCalculationOptions,
    BackgroundCalculationResults,
)
from negative_reflection import negative_reflection
from utility import convolutional_test, whole_scene_test

def _normality_test(data: np.ndarray) -> np.ndarray:
    """
    The actual test, given a dataset centered around 0, test for normality.
    Matches form expected by `whole_scene_test` and `convolutional_test` in utility.py.
    Inputs:
        data: any shape, should be centered on 0, likely from calling negative_reflection
    Output:
        [statistic, p_value, implied sigma].
        Statistic is 0 for perfectly normal and higher for less normal.
    """
    if data.size < 8:
        return np.array([np.nan, np.nan, np.nan])
    # This has some deprecated numpy usage.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        stat, p = stats.normaltest(data)
    sigma = np.std(data)
    return np.array([stat, p, sigma])


def global_normality_test(
    l3_data: np.ndarray, mean_ak: np.ndarray, options: BackgroundCalculationOptions
) -> tuple[BackgroundCalculationResults, BackgroundCalculationAnalytics]:
    """
    Given some L3 xCH4 observations and the mean averaging kernel at those points, estimate the
    background using a global normality test.

    Inputs:
        l3_data: The raw L3 xch4 offset from prior. Any shape
        mean_ak: The column averaged averaging kernel. Same shape as l3_data.
        options: Configuration for the test. Note that we ignore convolution info.
    Returns:
        Background estimate results, and analytics, showing the normalcy metric (lower is better)
        of tested backgrounds.
    """
    bgs_to_test, results = whole_scene_test(
        l3_data,
        mean_ak,
        lambda l3, ak, bg: _normality_test(negative_reflection(l3, ak, bg)),
        options.search_start,
        options.search_stop,
        options.search_steps,
    )
    metrics, ps, sigmas = results[..., 0], results[..., 1], results[..., 2]
    # Normaltest statistics is best when as small as possible.
    best_result = np.nanmin(metrics)
    # This stat can get very small near the right answer (<1e-3). This line basically considers
    # any bg viable if its metric is within 10% or 10 (absolute value) of the best metric.
    viable_bgs = bgs_to_test[(metrics - best_result) / max(best_result, 100) < 0.1]
    best_bg = bgs_to_test[np.nanargmin(metrics)]
    spread = np.nanmax(viable_bgs) - np.nanmin(viable_bgs)

    return BackgroundCalculationResults(best_bg, spread / 2), BackgroundCalculationAnalytics(
        bgs_to_test.tolist(), metrics.tolist(), ps.tolist(), sigmas.tolist()
    )


def local_normality_test(
    l3_data: np.ndarray, mean_ak: np.ndarray, options: BackgroundCalculationOptions
) -> tuple[BackgroundCalculationResults, BackgroundCalculationAnalytics]:
    """
    Given some L3 xCH4 observations and the mean averaging kernel at those points, estimate the
    background using a local normality test.

    Inputs:
        l3_data: The raw L3 xch4 offset from prior. Any shape
        mean_ak: The column averaged averaging kernel. Same shape as l3_data.
        options: Configuration for the test. Note that we assume convolution info is present.
    Returns:
        Background estimate results, and analytics, showing the normalcy metric (lower is better)
        of tested backgrounds.
    """
    bgs_to_test, results = convolutional_test(
        l3_data,
        mean_ak,
        lambda l3, ak, bg: _normality_test(negative_reflection(l3, ak, bg)),
        options.search_start,
        options.search_stop,
        options.search_steps,
        options.region_size,
        options.region_spacing,
    )

    # At the moment, the convolutional version doesn't care about the p-values or std's.
    metrics, _, _ = results[..., 0], results[..., 1], results[..., 2]

    # Normaltest values that are low indicate more normal data, so for each region, figure out what
    # background resulted in the smallest statistic.
    min_bg_index = np.full(metrics.shape[:2], np.nan)
    for x in range(metrics.shape[0]):
        for y in range(metrics.shape[1]):
            # Get non-nan results for this region and take the best bg index.
            if np.any(np.isfinite(metrics[x, y])):
                min_bg_index[x, y] = np.nanargmin(metrics[x, y])
            else:
                min_bg_index[x, y] = np.nan
    # Now take the 25th percentile (somewhat arbitrary) to exclude low lying outliers.
    best_index = int(np.nanpercentile(min_bg_index, 25))
    best_bg = bgs_to_test[best_index]
    # The STD of this result is from 5%-45% (also somewhat arbitrary).
    min_index, max_index = np.nanpercentile(min_bg_index, 5), np.nanpercentile(min_bg_index, 45)
    possible_bgs = bgs_to_test[int(min_index) : int(max_index)]
    bg_std = (np.nanmax(possible_bgs) - np.nanmin(possible_bgs)) / 2

    # Calculate the vote share for every background.
    vote_share = [np.count_nonzero(min_bg_index == i) for i in range(len(bgs_to_test))]
    vote_share = np.array(vote_share, dtype=float)
    vote_share /= np.sum(vote_share)

    return BackgroundCalculationResults(best_bg, bg_std), BackgroundCalculationAnalytics(
        bgs_to_test.tolist(), vote_share.tolist()
    )
