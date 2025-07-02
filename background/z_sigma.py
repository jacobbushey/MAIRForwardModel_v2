from scipy.ndimage import uniform_filter, binary_fill_holes
from sklearn.cluster import DBSCAN
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

#from msat_level4_gim.background.data_types import (
#    BackgroundCalculationAnalytics,
#    BackgroundCalculationOptions,
#    BackgroundCalculationResults,
#)
#from msat_level4_gim.background.negative_reflection import negative_reflection
#from msat_level4_gim.background.plotting import plot_zsigma_overview

from data_types import (
    BackgroundCalculationAnalytics,
    BackgroundCalculationOptions,
    BackgroundCalculationResults,
)
from negative_reflection import negative_reflection
from plotting import plot_zsigma_overview

def nan_aware_local_std(arr, size):
    """Compute the local standard deviation of a 2D array using a uniform filter.

    Args:
        arr (np.ndarray): (n,m) L3 data
        size (int): Size of the local window is (size x size) with pixel count size*size

    Returns:
        local_std (np.ndarray): (n,m) Local standard deviation of the input array
        local_count (np.ndarray): (n,m) Count of valid pixels in the local window
        local_mean (np.ndarray): (n,m) Local mean of the input array
    """
    # NaN mask: 1 where valid, 0 where NaN
    valid_mask = np.isfinite(arr).astype(float)

    # Replace NaNs with 0 for sum computation (will be masked later)
    arr_filled = np.nan_to_num(arr, nan=0.0)

    # Sum and count in the window
    local_mean = uniform_filter(arr_filled, size=size, mode="reflect")
    local_valid_fraction = uniform_filter(valid_mask, size=size, mode="reflect")
    local_sq_mean = uniform_filter(arr_filled**2, size=size, mode="reflect")

    # Remove regions with (nearly) no valid pixels
    local_mean[local_valid_fraction < 1e-5] = np.nan
    local_sq_mean[local_valid_fraction < 1e-5] = np.nan
    local_valid_fraction[local_valid_fraction < 1e-5] = np.nan

    # adjust for local count values (since we artificially added 0s in the mean calculation)
    local_mean /= local_valid_fraction
    local_sq_mean /= local_valid_fraction

    # Std = sqrt(E[x^2] - (E[x])^2)
    local_std = np.sqrt(local_sq_mean - local_mean**2)
    return local_std, local_valid_fraction, local_mean


def get_largest_clusters(binary_mask, eps=1.5, min_samples=6, top_n=5):
    """
    Takes a binary mask of 1s (neg. residual) and 0s (pos. residual) and uses DBSCAN to find
    the largest connected clusters of negative residuals.

    Returns the cluster IDs, sizes, and densities of the top N clusters. The density is defined as
    the number of pixels in the cluster divided by the area of the cluster (in pixels), where the
    area is is the cluster size after filling holes. It remains unused for now, but could be
    useful for future analysis.

    Inputs:
        binary_mask: 2D numpy array of 1s (foreground) and 0s (background)
        eps: DBSCAN eps parameter (in pixels)
        min_samples: DBSCAN min_samples parameter
        top_n: number of top clusters to analyze by size

    Returns:
        Dict with cluster IDs, sizes, and densities fot the top N clusters.
        2D numpy int array with 1's denoting the largest cluster.
    """

    # Step 1: Get coordinates of 1s (being negative residuals)
    coords = np.column_stack(np.nonzero(binary_mask))

    if len(coords) == 0:
        print("No foreground pixels found.")
        return {"cluster_id": [], "size": [], "density": []}

    # Step 2: Run DBSCAN
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(coords)
    labels = db.labels_

    labels_on_map = np.full(binary_mask.shape, np.nan)
    labels_on_map[coords[:, 0], coords[:, 1]] = labels

    # Step 3: Count cluster sizes (ignore noise = label -1)
    counts = Counter(labels)
    if -1 in counts:
        del counts[-1]

    # Step 4: Get top N clusters by size
    top_clusters = counts.most_common(top_n)
    results = {"cluster_id": [], "size": [], "density": []}
    # Step 5: Calculate size and density for each cluster
    for cluster_id, size in top_clusters:
        cluster_coords = (labels_on_map == cluster_id).astype(int)
        size = np.sum(cluster_coords)
        area = binary_fill_holes(cluster_coords).sum()

        density = size / area if area and not np.isnan(area) else np.nan

        results["cluster_id"].append(cluster_id)
        results["size"].append(float(size))
        results["density"].append(density)

    return results, np.where(labels_on_map == top_clusters[0][0], 1, 0)


def _get_binary_sign_mask(ch4_offset, mean_ak, bg, mask):
    """
    Get the binary mask of the sign of the enhancement. -1 is 1, +1 is 0.
    """
    sign_of_enhancement = np.sign(np.where(mask, ch4_offset - bg * mean_ak, np.nan))
    # replace -1 with 1 since we look out for those and 1 with 0
    return np.where(sign_of_enhancement == -1, 1, 0)


def z_sigma_fit(
    ch4_offset: np.ndarray, mean_ak: np.ndarray, options: BackgroundCalculationOptions
) -> tuple[BackgroundCalculationResults, BackgroundCalculationAnalytics]:
    """
    Given some L3 xCH4 observations and the mean averaging kernel at those points, estimate the
    background using the Z-Sigma method.
        1.  Calculate the local std and mean for the whole scene using a uniform filter.
        2.  For each background guess, take the negative residuals, divide them by the local std,
            and reflect them around 0.
        3.  Calculate the std of the reflected distribution. The best background is the one with
            the std closest to 1.
    The background may, depending on the data, have an obvious bias. Connected patches of negative
    residuals are a quality criterion for the background estimate. If the largest cluster is very
    small, the background is probably underestimated. If the largest cluster is very large, the
    background is probably overestimated. We can try to salvage the background estimate by
    iteratively optimizing the background guess:
        1.  Large clusters: Probably due to enhanced regions contributing to negative residuals.
            --> Iteratively remove data in regions with high local means. An area with an enhanced
            local mean is less likely to be a background region, but can have negative residuals
            due to noise. Since we remove data from high-mean regions, this should be fine.
        2.  Small clusters: Probably due to low-lying outliers in the scene.
            --> Iteratively remove the smallest values from the enhancements. This is very sketchy
            and I have much less confidence in this approach than in the large cluster approach.

    Inputs:
        ch4_offset: The raw L3 xch4 offset from prior. Any shape
        mean_ak: The column averaged averaging kernel. Same shape as ch4_offset.
        options: BackgroundCalculationOptions object with parameters for the search.
            - search_start: Optional starting point for search
            - search_stop: Optional end point for search
            - search_steps: Number of steps to take between start and stop.
            - region_size: Edge size of the square convolution kernel, defining the neighborhood.
            - dbscan_epsilon: DBSCAN eps parameter (in pixels)
            - dbscan_min_samples: DBSCAN min_samples parameter
            - dbscan_top_n: number of top clusters to analyze by size
            - max_iterations: maximum number of iterations for finding a suitable background
                estimate. Set to 1 to turn iterations off.
            - dbscan_max_cluster_size: maximum size for a valid background estimate.
            - dbscan_min_cluster_size: minimum size for a valid background estimate.
            - plot_dir: Directory to save plots to. No plots are saved if None.
            - plot_name: Name of the plot. If None, the default name is used.

    Returns:
        BackgroundCalculationResults: The best background estimate and its standard deviation.
            The standard deviation is not available for this method and is set to 0.
        BackgroundCalculationAnalytics: Additional information about the background calculation
            process.
    """

    # calculate the local std and mean for the whole scene
    local_std, local_valid_fraction, local_mean = nan_aware_local_std(
        ch4_offset, options.region_size
    )

    # Create a mask for valid pixels // TODO: some more investigations needed here
    mask = np.logical_and(
        local_std < np.nanmedian(local_std) * 3,
        local_valid_fraction > 0.75,
    )

    # Assume the correct background is between the 1% and 60%. We actually expect it to be <50%.
    minimum_guess = (
        options.search_start if options.search_start else np.nanpercentile(ch4_offset, 1)
    )
    maximum_guess = options.search_stop if options.search_stop else np.nanpercentile(ch4_offset, 60)
    bgs_to_test = np.linspace(minimum_guess, maximum_guess, options.search_steps)

    # have some lists to store the best background and the max cluster size during iterations
    best_bgs = []
    largest_cluster_sizes = []
    # initialize parameters to enter loop at least once
    largest_cluster = options.dbscan_max_cluster_size + 1
    iteration = 0
    # While loop continues until the largest_cluster size is within the limits
    # [options.dbscan_min_cluster_size, options.dbscan_max_cluster_size] or the maximum number
    # of iterations is reached.
    while (
        largest_cluster > options.dbscan_max_cluster_size
        or largest_cluster < options.dbscan_min_cluster_size
    ) and iteration < options.max_iterations:

        # IF we enter this part again, start adjusting the data mask
        if largest_cluster > options.dbscan_max_cluster_size:
            # iteratively remove data in regions with high local means
            percentile = 100 - iteration * 5
            if percentile <= 0:
                break
            mask = np.logical_and(mask, local_mean < np.nanpercentile(local_mean, percentile))
        elif largest_cluster < options.dbscan_min_cluster_size:
            # iteratively remove the smallest values from the enhancements in 0.05% increments
            mask = np.logical_and(mask, ch4_offset > np.nanpercentile(ch4_offset, iteration * 0.05))

        #
        std_below_zero = np.full(bgs_to_test.shape, np.nan)
        for i, bg in enumerate(bgs_to_test):
            dev_in_std = (ch4_offset[mask] - bg * mean_ak[mask]) / local_std[mask]
            neg_reflected_distr = negative_reflection(dev_in_std[dev_in_std <= 0], 0, 0)
            if neg_reflected_distr.size > 1000:  # arbitrary threshold, usually much larger on L3
                std_below_zero[i] = np.nanstd(neg_reflected_distr)
            else:
                std_below_zero[i] = np.nan

        # Find the best background, which is the one with the std closest to 1
        best_bg = bgs_to_test[np.nanargmin(np.abs(std_below_zero - 1))]

        # A placeholder estimate for the standard deviation for this method. Depends on search_steps
        top_bgs = bgs_to_test[
            np.abs(std_below_zero - 1)
            <= np.quantile(np.abs(std_below_zero - 1), options.top_bgs_quantile)
        ]
        bg_std = np.std(top_bgs)

        # find the connected negative residual clusters
        binary_mask = _get_binary_sign_mask(ch4_offset, mean_ak, best_bg, mask)
        cluster_results, largest_cluster_map = get_largest_clusters(
            binary_mask,
            eps=options.dbscan_epsilon,
            min_samples=options.dbscan_min_samples,
            top_n=options.dbscan_top_n,
        )

        largest_cluster = max(cluster_results["size"])
        largest_cluster_sizes.append(largest_cluster)
        best_bgs.append(float(best_bg))
        iteration += 1

    if options.plot_dir:
        try:
            fig = plot_zsigma_overview(
                ch4_offset,
                mean_ak,
                best_bg,
                binary_mask,
                local_std,
                largest_cluster_map,
                bgs_to_test,
                std_below_zero,
                best_bgs,
                largest_cluster_sizes,
                mask,
            )
            if options.plot_name is None:
                options.plot_name = "z_sigma"
            fig.savefig(
                Path(options.plot_dir) / f"bg_overview_{options.plot_name}.png",
                dpi=300,
                bbox_inches="tight",
            )
            plt.close(fig)
        except Exception as e:
            print(f"Error while plotting: {e}")

    return BackgroundCalculationResults(best_bg, bg_std), BackgroundCalculationAnalytics(
        background_candidates=bgs_to_test.tolist(),
        background_candidate_calculated_standard_deviation=std_below_zero.tolist(),
        cluster_sizes=cluster_results["size"],
        largest_cluster_sizes=largest_cluster_sizes,
        best_bgs=best_bgs,
    )
