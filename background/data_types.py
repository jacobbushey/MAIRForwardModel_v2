from dataclasses import dataclass
from enum import auto
from typing import Optional

from mashumaro.mixins.json import DataClassJSONMixin
from strenum import StrEnum


class BackgroundCalculationMethod(StrEnum):
    NORMALITY = auto()
    FIT_STD = auto()
    ZSIGMA = auto()


@dataclass
class BackgroundCalculationOptions(DataClassJSONMixin):
    # Technique to use.
    method: BackgroundCalculationMethod = BackgroundCalculationMethod.NORMALITY

    # Filtering L3 data parameters
    num_samples_threshold: float = 0.5
    # Optional parameters that can fix the bg search space.
    search_start: Optional[float] = None
    search_stop: Optional[float] = None
    search_steps: int = 200

    # If doing a Convolution, these should be filled.
    region_size: Optional[int] = None  # How large in Manhattan distance radius is the region?
    region_spacing: Optional[int] = None  # How many pixels should we leave between region centers?

    # Options for Z-Sigma test
    # These are default parameters for the DBSCAN clustering algorithm.
    dbscan_epsilon: Optional[float] = 1.5
    dbscan_min_samples: Optional[int] = 6
    dbscan_top_n: Optional[int] = 3
    # These are parameters for iterative corrections if unreasonable clusters sizes
    max_iterations: int = 1  # set to 1 to turn iterations off // recommend 10-20
    # For 45m L3 pixels
    dbscan_max_cluster_size: Optional[int] = 5e5  # Typical largest cluster size for desert scenes
    dbscan_min_cluster_size: Optional[int] = 1e4  # Obvious background underestimations cluster size
    # Upper quantile limit to use when choosing backgrounds with normalized sigma
    # close to 1 to compute std for the method
    top_bgs_quantile: Optional[int] = 0.1

    plot_dir: Optional[str] = None  # Directory to save plots to
    plot_name: Optional[str] = None  # Name of the plot


@dataclass
class BackgroundCalculationResults(DataClassJSONMixin):
    # At the moment, we assume that all background calculation methods estimate a normal
    # distribution and return its mean and standard deviation.
    mean: float
    std: float


@dataclass
class BackgroundCalculationAnalytics(DataClassJSONMixin):
    """ "
    This class is not used for follow up work, but can be nice for plotting. Should be extended
    and every object should have a default value, in case some methods don't provide it."
    """

    # If other background values were tried, report their metrics. These lists should be the same
    # length if they are not None.
    background_candidates: Optional[list[float]] = None
    background_candidate_statistic: Optional[list[float]] = None
    background_candidate_p_values: Optional[list[float]] = None
    # This is our estimate of the noise model std.
    background_candidate_calculated_standard_deviation: Optional[list[float]] = None

    # These are diagnostic fields filled in automatically for every scene
    mean_albedo: Optional[float] = None
    mean_sza: Optional[float] = None

    # Additions for Z-Sigma method
    # Cluster sizes and densities of last iteration
    cluster_densities: Optional[list[float]] = None
    cluster_sizes: Optional[list[float]] = None
    # Largest cluster sizes and best background candidates during iterations
    largest_cluster_sizes: Optional[list[float]] = None
    # Best background for each iteration
    best_bgs: Optional[list[float]] = None
