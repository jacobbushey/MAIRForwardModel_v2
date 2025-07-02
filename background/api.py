from importlib import resources
from pathlib import Path
from typing import Optional

import numpy as np
import yaml
#from msat_netcdf.schema_validator import validate_against_schema
from schema_validator import validate_against_schema
from netCDF4 import Dataset
#from shapely import MultiPolygon, Polygon, contains_xy, prepare
from shapely.geometry import MultiPolygon, Polygon
from shapely import contains_xy, prepare

#from msat_level4_gim.background.data_types import (
#    BackgroundCalculationAnalytics,
#    BackgroundCalculationMethod,
#    BackgroundCalculationOptions,
#    BackgroundCalculationResults,
#)
#from msat_level4_gim.background.fit_std import global_std_fit, local_std_fit
#from msat_level4_gim.background.normality import global_normality_test, local_normality_test
#from msat_level4_gim.background.z_sigma import z_sigma_fit
#from msat_level4_gim.errors import ConfigurationException, SchemaValidationException

from data_types import (
    BackgroundCalculationAnalytics,
    BackgroundCalculationMethod,
    BackgroundCalculationOptions,
    BackgroundCalculationResults,
)
from fit_std import global_std_fit, local_std_fit
from normality import global_normality_test, local_normality_test
from z_sigma import z_sigma_fit
from errors import ConfigurationException, SchemaValidationException

def estimate_background_for_l3(
    l3_path: Path, options: BackgroundCalculationOptions, regions: Optional[list[Polygon]] = None
) -> tuple[BackgroundCalculationResults, BackgroundCalculationAnalytics]:
    """
    This function encompasses many methods for estimating background conentration and is the main
    entrypoint to this codebase from the broader L4 GIM perspective.
    Inputs:
        l3_path: Path to a local L3 Regrid file.
        options: Options that select a background estimation method and control it.
        regions: Optional list of regions to filter. Expected in WGS84 Lat/Lon.
    Output:
        tuple of result and analytics objects. The result is expected to get used in later stages,
        while the analytics is likely to be written to a file for manual analysis.
    """
    # Verify region parameters are either both set or both unset.
    if (options.region_size is None) != (options.region_spacing is None):
        raise ConfigurationException("Region Size and Region Spacing must both be None or valued.")
    is_convolution = not (options.region_size is None)

    # First verify the file matches expected L3 schema and extract necessary data from it.
    #schema = yaml.safe_load(
    #    resources.files("msat_level4_gim").joinpath("schemas/l3_for_background.yaml").read_text()
    #)
    with open("l3_for_background.yaml", "r") as f:
        schema = yaml.safe_load(f)
    with Dataset(str(l3_path), "r") as ds:
        valid, errors = validate_against_schema(ds, schema)
        if not valid:
            raise SchemaValidationException(f"L3 does not match background schema: {errors}")

        raw_xch4 = np.where(
            ds["num_samples"][:] > options.num_samples_threshold,
            ds["xch4"][:].filled(np.nan),
            np.nan,
        )
        averaging_kernel = ds["co2proxy_fit_diagnostics/ch4_averaging_kernel_mean"][:].filled(
            np.nan
        )
        prior = ds["apriori_data/xch4"][:].filled(np.nan)
        # Calculate the diagnostics immediately.
        albedo = ds["apriori_data/albedo_ch4band"][:].filled(np.nan)
        sza = ds["geolocation/sza"][:].filled(np.nan)
        # Don't bother reading the lat/long unless they are going to be used (for region filtering).
        if regions:
            lon = ds["lon"][:].filled(np.nan)
            lat = ds["lat"][:].filled(np.nan)
            # L3 is lat/lon order.
            lat, lon = np.meshgrid(lat, lon, indexing="ij")

    # See README for details, but we essentially want to model the offset between the observed
    # XCH4 and the prior. Copy is necessary to avoid read-only arrays.
    xch4_offset = np.copy(raw_xch4 - prior)

    # Filter by regions if they are provided.
    if regions:
        if len(regions) > 1:
            # For speed, convert to a multipolygon, if possible
            regions = MultiPolygon(regions)
        elif len(regions) == 1:
            regions = regions[0]
        else:
            raise ValueError("Regions must be a non-empty list of polygons.")

        prepare(regions)
        in_region_mask = contains_xy(regions, lon, lat)
        xch4_offset = np.where(in_region_mask, xch4_offset, np.nan)

    # Switch on the method selected.
    if options.method == BackgroundCalculationMethod.NORMALITY:
        if is_convolution:
            result, analytics = local_normality_test(xch4_offset, averaging_kernel, options)
        else:
            result, analytics = global_normality_test(xch4_offset, averaging_kernel, options)
    elif options.method == BackgroundCalculationMethod.FIT_STD:
        if is_convolution:
            result, analytics = local_std_fit(xch4_offset, averaging_kernel, options)
        else:
            result, analytics = global_std_fit(xch4_offset, averaging_kernel, options)
    elif options.method == BackgroundCalculationMethod.ZSIGMA:
        if is_convolution:
            result, analytics = z_sigma_fit(xch4_offset, averaging_kernel, options)
        else:
            raise ConfigurationException("ZSIGMA method requires region size be set.")
    # elif options.method == NEXT METHOD
    else:
        raise NotImplementedError(f"Background method {options.method} is not supported.")

    # Do some common operations to add to analytics. Force floats for yaml conversion.
    analytics.mean_albedo = float(np.nanmean(albedo))
    analytics.mean_sza = float(np.nanmean(sza))
    return result, analytics
