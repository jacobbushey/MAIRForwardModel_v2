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


# Set the path to the file you're interested in reading -----------------------------

# MAIR_RF06_Permian +++++++++++++++++++++++++++++++++++++

path_to_file = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/Inputs/L3/MAIR_RF06_Permian/MethaneAIR_L3_mosaic_20210806T161742_20210806T183024_dpp_sza_ak.nc"
  # took longer than I would've expected?
  # maybe just taking a little long to converge?
  
  # should be a similar file size. Why is it taking so long? Weird shapes? But NZ has weird shapes too.

# MSAT_005_Permian ++++++++++++++++++++++++++++++++++++++

#path_to_file = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/Inputs/L3/MSAT_005_01430050/regrid_t5_20240911_c01430050_p738_45m_MSAT_L3_45m_c01430050_p738_v01010000_20240911T203025Z_203057Z.nc"
  # fully sized MethaneSAT scenes (no water) take longer (minutes)

path_to_file = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/Inputs/L3/MSAT_005_01430050/MSAT_L3_45m_c01430050_p738_v01010000_20240911T203025Z_203057Z.nc"

path_to_geojson = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/MSATForwardModel/background/geojsons/MSAT_L3_45m_c01430050_p738_v01010000_20240911T203025Z_203057Z.geojson"

# MSAT_196_Canterbury_1 ++++++++++++++++++++++++++++++++++++

#path_to_file = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/Inputs/L3/MSAT_196_07EE0C40/regrid_t196_20250301_c07EE0C40_p1339_45m_MSAT_L3_45m_c07EE0C40_p1339_v01010001_20250301T015137Z_015209Z.nc"

path_to_file = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/Inputs/L3/MSAT_196_07EE0C40/MSAT_L3_45m_c07EE0C40_p1339_v01010001_20250301T015137Z_015209Z.nc"

path_to_geojson = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/MSATForwardModel/background/geojsons/MSAT_L3_45m_c07EE0C40_p1339_v01010001_20250301T015137Z_015209Z.geojson"

# MSAT_196_Canterbury_2 ++++++++++++++++++++++++++++++++++++

#path_to_file = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/Inputs/L3/MSAT_196_09030C40/regrid_t196_20250320_c09030C40_p1399_45m_MSAT_L3_45m_c09030C40_p1399_v01011000_20250320T020530Z_020602Z.nc"

path_to_file = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/Inputs/L3/MSAT_196_09030C40/MSAT_L3_45m_c09030C40_p1399_v01011000_20250320T020530Z_020602Z.nc"

path_to_geojson = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/MSATForwardModel/background/geojsons/MSAT_L3_45m_c09030C40_p1399_v01011000_20250320T020530Z_020602Z.geojson"

# MSAT_196_Canterbury_3 ++++++++++++++++++++++++++++++++++++

#path_to_file = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/Inputs/L3/MSAT_196_09230C40/regrid_t196_20250321_c09230C40_p1388_45m_MSAT_L3_45m_c09230C40_p1388_v01011000_20250321T021106Z_021138Z.nc"

path_to_file = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/Inputs/L3/MSAT_196_09230C40/MSAT_L3_45m_c09230C40_p1388_v01011000_20250321T021106Z_021138Z.nc"

path_to_geojson = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/MSATForwardModel/background/geojsons/MSAT_L3_45m_c09230C40_p1388_v01011000_20250321T021106Z_021138Z.geojson"

# MSAT_065_Turkmenistan
#path_to_file = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/Inputs/L3/MSAT_065_01A30410/regrid_t65_20240923_c01A30410_p853_45m_MSAT_L3_45m_c01A30410_p853_v01010000_20240923T104207Z_104239Z.nc"

path_to_file = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/Inputs/L3/MSAT_065_01A30410/MSAT_L3_45m_c01A30410_p853_v01010000_20240923T104207Z_104239Z.nc"

path_to_geojson = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/MSATForwardModel/background/geojsons/MSAT_L3_45m_c01A30410_p853_v01010000_20240923T104207Z_104239Z.geojson"



# Loop through all of these
# make the output figure title the filename
path_to_file = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model_v2/Inputs/L3/RF06_Permian/MethaneAIR_L3_segment_20210806T161742_20210806T162243_dpp_ak.nc"

filepath = '/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model_v2/Inputs/L3/RF06_Permian/'

file_list = ['MethaneAIR_L3_segment_20210806T162243_20210806T162745_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T162745_20210806T163247_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T163247_20210806T163748_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T163748_20210806T164250_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T164250_20210806T164751_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T164751_20210806T165253_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T165253_20210806T165754_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T165754_20210806T170256_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T170256_20210806T170758_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T170758_20210806T171259_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T171259_20210806T171801_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T171801_20210806T172302_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T172302_20210806T172804_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T172804_20210806T173306_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T173306_20210806T173807_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T173807_20210806T174309_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T174309_20210806T174810_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T174811_20210806T175312_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T175312_20210806T175814_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T175814_20210806T180315_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T180315_20210806T180747_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T180747_20210806T181218_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T181218_20210806T181650_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T181650_20210806T182121_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T182121_20210806T182553_dpp_ak.nc',
    'MethaneAIR_L3_segment_20210806T182553_20210806T183024_dpp_ak.nc']
    # find a way to retrieve this list according to 'ak.nc'

for i in range(len(file_list)):
    file_tick = file_list[i]
    filepath_tick = filepath + file_tick

    plot_name_tick = 'TEST' + file_tick

    options = BackgroundCalculationOptions(
        method=BackgroundCalculationMethod.ZSIGMA,
        num_samples_threshold=0.50,
        region_size=20,
        region_spacing=20,
    #    dbscan_max_cluster_size=5e5,
    #    dbscan_min_cluster_size=1e4,
        plot_dir="/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model_v2/MAIRForwardModel_v2/background",
        plot_name=plot_name_tick
    )

    result, analytics = estimate_background_for_l3(path_to_file, options)


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

options = BackgroundCalculationOptions(
    method=BackgroundCalculationMethod.ZSIGMA,
    num_samples_threshold=0.50,
    region_size=20,
    region_spacing=20,
    dbscan_max_cluster_size=5e5,
    dbscan_min_cluster_size=1e4,
    plot_dir="/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model_v2/MAIRForwardModel_v2/background",
    plot_name="TEST_MethaneAIR_L3_segment_20210806T161742_20210806T162243_dpp.nc"
)

# Run the function ----------------------------------------------------

#result, analytics = estimate_background_for_l3(path_to_file, options, regions)
result, analytics = estimate_background_for_l3(path_to_file, options)
  # currently no cropping with geojsons implemented for MAIR

# Save the output as a csv to be read into other files ----------------------

with open('testfile.csv', 'rb') as f:
    data = list(csv.reader(f))








Traceback (most recent call last):
  File "<python-input-154>", line 1, in <module>
    result, analytics = estimate_background_for_l3(path_to_file, options)
                        ~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^
  File "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/MSATForwardModel/background/api.py", line 116, in estimate_background_for_l3
    result, analytics = z_sigma_fit(xch4_offset, averaging_kernel, options)
                        ~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/MSATForwardModel/background/z_sigma.py", line 242, in z_sigma_fit
    cluster_results, largest_cluster_map = get_largest_clusters(
                                           ~~~~~~~~~~~~~~~~~~~~^
        binary_mask,
        ^^^^^^^^^^^^
    ...<2 lines>...
        top_n=options.dbscan_top_n,
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^
    )
    ^
  File "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/MSATForwardModel/background/z_sigma.py", line 90, in get_largest_clusters
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(coords)
  File "/n/home03/jbushey/.conda/envs/jacob_env/lib/python3.13/site-packages/sklearn/base.py", line 1363, in wrapper
    return fit_method(estimator, *args, **kwargs)
  File "/n/home03/jbushey/.conda/envs/jacob_env/lib/python3.13/site-packages/sklearn/cluster/_dbscan.py", line 421, in fit
    neighborhoods = neighbors_model.radius_neighbors(X, return_distance=False)
  File "/n/home03/jbushey/.conda/envs/jacob_env/lib/python3.13/site-packages/sklearn/neighbors/_base.py", line 1277, in radius_neighbors
    chunked_results = Parallel(n_jobs, prefer="threads")(
        delayed_query(X[s], radius, return_distance, sort_results=sort_results)
        for s in gen_even_slices(X.shape[0], n_jobs)
    )
  File "/n/home03/jbushey/.conda/envs/jacob_env/lib/python3.13/site-packages/sklearn/utils/parallel.py", line 82, in __call__
    return super().__call__(iterable_with_config_and_warning_filters)
           ~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/n/home03/jbushey/.conda/envs/jacob_env/lib/python3.13/site-packages/joblib/parallel.py", line 1986, in __call__
    return output if self.return_generator else list(output)
                                                ~~~~^^^^^^^^
  File "/n/home03/jbushey/.conda/envs/jacob_env/lib/python3.13/site-packages/joblib/parallel.py", line 1914, in _get_sequential_output
    res = func(*args, **kwargs)
  File "/n/home03/jbushey/.conda/envs/jacob_env/lib/python3.13/site-packages/sklearn/utils/parallel.py", line 147, in __call__
    return self.function(*args, **kwargs)
           ~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^
KeyboardInterrupt










