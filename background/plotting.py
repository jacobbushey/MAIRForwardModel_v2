import matplotlib.pyplot as plt
import numpy as np


def plot_zsigma_overview(
    ch4_ofs,
    mean_ak,
    best_bg,
    binary_mask,
    local_std,
    largest_cluster_map,
    bgs_to_test,
    std_below_zero,
    best_bgs,
    max_cluster_sizes,
    mask,
):
    """
    plot_zsigma_overview _summary_

    Args:
        ch4_ofs: [n x m] The raw L3 xch4 offset from prior for L3 dimensions  n x m
        mean_ak: [n x m] The column averaged averaging kernel. Same shape as ch4_offset
        best_bg: scalar, the background with the std closest to 1
        binary_mask: [n x m] int, 1's denote background pixels
        local_std: [n x m] local std for the whole scene
        largest_cluster_map: [n x m] int array with 1's denoting the largest cluster
        bgs_to_test: [k] backgrounds to test
        std_below_zero: [k] for each background, the std of the negative residuals (RNN)
        best_bgs: list, length >= 1 of the best backgrounds for each (optional) iteration.
        max_cluster_sizes: list, length >= 1 of the size of the largest cluster for each iteration.
        mask: [n x m]: Outlier mask, denotes candidate pixels for background selection.

    Returns:
        6-panel plot
    """
    # compute stuff
    predicted_enhancements = ch4_ofs + best_bg * mean_ak
    background_evolution = [bg - best_bg for bg in best_bgs]

    fig, axs = plt.subplots(3, 2, figsize=(12, 12), dpi=300)

    x = predicted_enhancements[mask]
    vmin, vmax = np.quantile(x[np.isfinite(x)], [0.05, 0.95])
    if not np.isfinite(vmin) or not np.isfinite(vmax):
        vmin, vmax = None, None
    im = axs[0, 0].imshow(
        predicted_enhancements, cmap="viridis", origin="lower", vmin=vmin, vmax=vmax
    )
    axs[0, 0].set_title("Predicted Enhancements")
    fig.colorbar(im, ax=axs[0, 0], label="Enhancement (ppb)")

    im = axs[0, 1].imshow(np.where(mask, local_std, np.nan), cmap="cividis", origin="lower")
    axs[0, 1].set_title("Local Standard Deviation")
    fig.colorbar(im, ax=axs[0, 1], label="Standard Deviation (ppb)")

    axs[1, 0].imshow(np.where(mask, binary_mask, np.nan), cmap="cividis", origin="lower")
    axs[1, 0].set_title("Binary Mask")

    axs[1, 1].imshow(np.where(mask, largest_cluster_map, np.nan), cmap="cividis", origin="lower")
    axs[1, 1].set_title(f"Largest Cluster: {largest_cluster_map.sum():.0f} px")

    axs[2, 0].plot(bgs_to_test, std_below_zero)
    axs[2, 0].axvline(best_bg, color="red", label=f"Best Background: {best_bg:.2f} ppb")
    axs[2, 0].axhline(1, color="k", linewidth=0.5)
    axs[2, 0].set_xlabel("Background Offset Estimate (ppb)")
    axs[2, 0].set_ylabel("Normalized Standard Deviation of NRR")
    axs[2, 0].legend()

    n, bins, _ = axs[2, 1].hist(predicted_enhancements.flatten(), bins=100, log=True, alpha=0.5)
    axs[2, 1].hist(predicted_enhancements[mask], bins=bins)
    axs[2, 1].axvline(0, color="black", label="Final Estimate")
    for bg, cl in zip(background_evolution, max_cluster_sizes):
        axs[2, 1].axvline(
            bg, color="red", label=f"Background: {bg:.2f} ppb / {cl:.0f} px", linestyle="--"
        )
    axs[2, 1].set_xlabel("Enhancement (ppb)")
    axs[2, 1].set_ylabel("Pixel count")
    axs[2, 1].legend(loc="lower right")
    plt.tight_layout()
    return fig
