from scipy.ndimage import uniform_filter, binary_fill_holes
from sklearn.cluster import DBSCAN
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from dataclasses import dataclass
from enum import auto
from typing import Optional
from mashumaro.mixins.json import DataClassJSONMixin
from strenum import StrEnum
from netCDF4 import Dataset
import os
import csv
from numba import njit
import geojson
import shapely
import xarray as xr
import fnmatch

# Needed to install mashumaro version 3.14 to fix a dependency issue
# Need to install shapely 2.0 to fix a dependency issue
  # but this wasn't working with some "pinning" errors with python v 3.9?

# Import custom functions from scripts -----------------------

from data_types import (
    BackgroundCalculationAnalytics,
    BackgroundCalculationOptions,
    BackgroundCalculationResults,
    BackgroundCalculationMethod
)
from negative_reflection import negative_reflection
from plotting import plot_zsigma_overview

from z_sigma import nan_aware_local_std
from z_sigma import get_largest_clusters
from z_sigma import _get_binary_sign_mask
from z_sigma import z_sigma_fit

from numba import njit
from fit_std import global_std_fit, local_std_fit

from api import estimate_background_for_l3

#def find_all(name, path):
#    result = []
#    for root, dirs, files in os.walk(path):
#        if name in files:
#            result.append(os.path.join(root, name))
#    return result

def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result
# from https://stackoverflow.com/questions/1724693/find-a-file-in-python

# Set the path to the file you're interested in reading -----------------------------

# Loop through all of these
# make the output figure title the filename
#path_to_file = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model_v2/Inputs/L3/RF06_Permian/MethaneAIR_L3_segment_20210806T161742_20210806T162243_dpp_ak.nc"

#path_to_file = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model_v2/Inputs/L3/RF06_Permian/MethaneAIR_L3_segment_20210806T161742_20210806T162243_dpp.nc"

filepath = '/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model_v2/Inputs/L3/RF06_Permian/'
name = 'RF06_Permian'

file_list = find("*dpp_ak.nc", filepath)
file_list.sort()

#file_list = ['MethaneAIR_L3_segment_20210806T161742_20210806T162243_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T162243_20210806T162745_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T162745_20210806T163247_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T163247_20210806T163748_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T163748_20210806T164250_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T164250_20210806T164751_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T164751_20210806T165253_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T165253_20210806T165754_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T165754_20210806T170256_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T170256_20210806T170758_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T170758_20210806T171259_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T171259_20210806T171801_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T171801_20210806T172302_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T172302_20210806T172804_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T172804_20210806T173306_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T173306_20210806T173807_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T173807_20210806T174309_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T174309_20210806T174810_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T174811_20210806T175312_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T175312_20210806T175814_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T175814_20210806T180315_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T180315_20210806T180747_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T180747_20210806T181218_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T181218_20210806T181650_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T181650_20210806T182121_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T182121_20210806T182553_dpp_ak.nc',
#    'MethaneAIR_L3_segment_20210806T182553_20210806T183024_dpp_ak.nc']
#    # find a way to retrieve this list according to 'ak.nc'

results = []
results.append(['segment_name', 'time_coverage_start', 'time_coverage_end', 'mean', 'sd'])
    # set the column headers so you don't lose the first row of data as headers

for i in range(len(file_list)):
    # remember, it's zero indexing!
    
    #file_tick = file_list[i]
    #filepath_tick = filepath + file_tick
    filepath_tick = file_list[i]

    #x = file_tick.split(".")
    ##plot_name_tick = 'TEST_' + x[0]
    #plot_name_tick = x[0]
    x = filepath_tick.split("/")
    plot_name_tick = x[-1]

    print(i)
    print(plot_name_tick)

    options = BackgroundCalculationOptions(
        method=BackgroundCalculationMethod.ZSIGMA,
        num_samples_threshold=0.50,
        region_size=20,
        region_spacing=20,  # I thought this parameter didn't get used? Still has to be set I guess
    #    dbscan_max_cluster_size=5e5,
    #    dbscan_min_cluster_size=1e4,
    #    dbscan_max_cluster_size=5e3,  # this doesn't seem to make any difference
    #    dbscan_min_cluster_size=1e2,
        plot_dir="/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model_v2/MAIRForwardModel_v2/background",
        plot_name=plot_name_tick
    )

    # changing the local_valid_fraction threshold didn't change anything
    # changing num_samples_threshold didn't change anything about this (why?)

    # WHY AM I GETTING MORE ERRORS AND A DIFFERENT OFFSET TODAY? THE FILE IS THE SAME RIGHT?
    # DIFFERENCE IN WHETHER I CHANGED SZA OR NOT?

    result, analytics = estimate_background_for_l3(filepath_tick, options)
    # result, analytics = estimate_background_for_l3(path_to_file, options)

    #result_mean = result.mean
    #result_std = result.std

    result_name_tick = plot_name_tick + '_result.csv'
    analytics_name_tick = plot_name_tick + '_analytics.csv'

    result_new = [result.mean, result.std]
    analytics_new = [
        [analytics.background_candidates],
        [analytics.background_candidate_statistic],
        [analytics.background_candidate_p_values],
        [analytics.background_candidate_calculated_standard_deviation],
        [analytics.mean_albedo],
        [analytics.mean_sza],
        [analytics.cluster_densities],
        [analytics.cluster_sizes],
        [analytics.largest_cluster_sizes],
        [analytics.best_bgs]
    ]

    nc_file = Dataset(filepath_tick, 'r')
    time_coverage_start = nc_file.time_coverage_start
    time_coverage_end = nc_file.time_coverage_end
    nc_file.close()

    with open(result_name_tick, mode = 'w', newline = '') as f:
        writer = csv.writer(f)
        writer.writerow(result_new)

    with open(analytics_name_tick, mode = 'w', newline = '') as f:
        writer = csv.writer(f)
        writer.writerow(analytics_new)

    # Adds the output to a growing list for all segments that will be saved at the end
    #results.append([file_tick, time_coverage_start, time_coverage_end, result.mean, result.std])
    results.append([plot_name_tick, time_coverage_start, time_coverage_end, result.mean, result.std])

filename = name + '_background_results.csv'
with open(filename, mode = 'w', newline = '') as f:
    writer = csv.writer(f)
    writer.writerows(results)

    #with open('data.csv', mode='w', newline='') as file:
    #writer = csv.writer(file)
    #writer.writerows(data)


# CURRENTLY NO CROPPING WITH GEOJSONS IMPLEMENTED FOR MAIR
# comes from msat_level4_gim/flyte/task_background.py
#with open(path_to_geojson, "r") as f:
#    geometry_collection = shapely.from_geojson(f.read())
#    regions = list(geometry_collection.geoms)

# Create the options object for the method you're using to solve -----------------

#options = BackgroundCalculationOptions(
#    method=BackgroundCalculationMethod.ZSIGMA,
#    num_samples_threshold=0.50,
#    search_start=None,
#    search_stop=None,
#    search_steps=200,
#    region_size=20,
#    region_spacing=10,
#    dbscan_epsilon=1.5,
#    dbscan_min_samples=6,
#    dbscan_top_n=3,
#    max_iterations=10,
#    dbscan_max_cluster_size=5e5,
#    dbscan_min_cluster_size=1e4,
#    top_bgs_quantile=0.1,
#    plot_dir="/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/MSATForwardModel/background",
#    plot_name="TEST_MAIR_RF06_Permian"
#)

#options = BackgroundCalculationOptions(
#    method=BackgroundCalculationMethod.ZSIGMA,
#    num_samples_threshold=0.50,
#    region_size=20,
#    region_spacing=20,
##    dbscan_max_cluster_size=5e5,
##    dbscan_min_cluster_size=1e4,
#    plot_dir="/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model_v2/MAIRForwardModel_v2/background",
#    plot_name="TEST_MethaneAIR_L3_segment_20210806T161742_20210806T162243_dpp.nc"
#)

# Run the function ----------------------------------------------------

##result, analytics = estimate_background_for_l3(path_to_file, options, regions)
#result, analytics = estimate_background_for_l3(path_to_file, options)
  # currently no cropping with geojsons implemented for MAIR

# Save the output as a csv to be read into other files ----------------------

#with open('testfile.csv', 'rb') as f:
#    data = list(csv.reader(f))





