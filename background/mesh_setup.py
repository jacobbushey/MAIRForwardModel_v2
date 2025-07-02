import math

import numpy as np
import numpy.typing as npt

#from msat_level4_gim.abstract_grid.data_types import AbstractGrid, GridGeometry
#from msat_level4_gim.abstract_grid.geometry import make_abstract_grid_gdf, make_gdf_from_coords
#from msat_level4_gim.inference.interpolate import simplex_info_for_points, triangulate
#from msat_level4_gim.inference.utils import distance_matrix

from abstract_grid_data_types import AbstractGrid, GridGeometry
from geometry import make_abstract_grid_gdf, make_gdf_from_coords
from interpolate import simplex_info_for_points, triangulate
from utils import distance_matrix

def create_bg_grid_for_observations(
    observation_grid: AbstractGrid, bg_grid_area: float
) -> AbstractGrid:
    """
    The background mesh is defined by a "bg_grid", which is a hexagonal grid of a user-defined
    resolution (`bg_grid_area`). The centroids of a hexagonal grid become the vertices in the
    triangular background mesh. This function defines the origin of the `bg_grid` such that all of
    these vertices will be in the first quadrant of the `bg_grid`.

    Inputs:
        observation_grid: AbstractGrid that the observations are on
        bg_grid_area: Area of each desired grid-cell (unit: m^2). A smaller value ultimately means a
            denser mesh.
    """

    # Make sure that the origin of the grid (the SW corner) is south-west of any of the hexagons
    # that are going to be used to build the mesh structure.
    buffer = 5 * math.sqrt(observation_grid.area)

    return AbstractGrid(
        epsg_code=observation_grid.epsg_code,
        origin_x=observation_grid.origin_x - buffer,
        origin_y=observation_grid.origin_y - buffer,
        geometry_type=GridGeometry.HEXAGON,
        area=bg_grid_area,
    )


def define_mesh(
    observation_coords: npt.NDArray[np.int_],
    observation_grid: AbstractGrid,
    bg_grid: AbstractGrid,
) -> tuple[
    npt.NDArray[np.float_],
    npt.NDArray[np.float_],
    npt.NDArray[np.int_],
]:
    """
    This function creates the geometry of an approximately-equal-area triangular mesh so that it
    completely covers the centroids of some user-specified observations (`observation_coords`) that
    exist on some `observation_grid`.

    - The area of each simplex in the triangular mesh is given by `bg_grid.area`. Use
      `create_bg_grid_for_observations()` to create `bg_grid`.

    - Each observation-grid-cell-centroid will be within one simplex, and it will have a barycentric
      coordinate within that simplex.

    - Each simplex will contain at least one observation-grid-cell-centroid

    Inputs:
        observation_grid_coords: (o,2) [int] AbstractGrid coordinates of each observation
        observation_grid: the AbstractGrid that the observations are on
        bg_grid: the AbstractGrid that defines the density of the mesh

    Outputs: 3 different sets of coordinates in different formats:
        bg_barycentric_coords: (o, v)
            A matrix, each row of the matrix is a set of barycentric coordinates for some
            observation `o`. (Therefore, all row-sums are equal to 1). 1-3 entries in each row
            will be nonzero, because an observation-grid-cell-centroid is in one and only one
            simplex (usually, 3 entries will be nonzero because the centroid will not be on the
            edges of the simplex). All column-sums are strictly positive, because each simplex in
            the mesh must contain at least one observation-grid-cell-centroid.

            Multiplying this matrix against a vector of vertex values performs a piecewise linear
            interpolation of the triangular mesh onto the observations.

        bg_utm_coords: (v, 2) UTM x/y coords for each of the `v` vertices in the mesh (unit: m)
            This is used for setting the distance-based correlation decay of the mesh prior.

        bg_grid_coords: (v, 2) AbstractGrid-coords for each of the `v` vertices, on `bg_grid`.
    """

    observation_domain_gdf = make_gdf_from_coords(coords=observation_coords, grid=observation_grid)

    # We don't yet know the set of `bg_grid` grid-cells that covers the observations.
    #
    # Our mesh-building strategy is to build geometries for `bg_grid` grid-cells in excess of what
    # we actually need (the "unclipped" geometries"), and then triangulate, and then figure out
    # which simplexes we need from that triangulation, and then "clip" out the vertices that weren't
    # involved in those simplices.
    buffer = 10 * math.sqrt(bg_grid.area)  # Much bigger buffer than we actually need here
    unclipped_bg_domain_gdf = make_abstract_grid_gdf(
        grid=bg_grid,
        max_utm_x=observation_domain_gdf.total_bounds[2] + buffer,
        max_utm_y=observation_domain_gdf.total_bounds[3] + buffer,
    )

    observation_domain_centroid_coords = np.stack(
        (
            observation_domain_gdf.centroid.get_coordinates()["x"].to_numpy(),
            observation_domain_gdf.centroid.get_coordinates()["y"].to_numpy(),
        ),
        axis=-1
    )

    # Triangulate the background domain. Many of the generated triangles are not needed because they
    # won't contain any observations; we will cut them out later.
    #
    # HACK: Add ~1m of jitter so that the triangulation succeeds. Triangulation will fail when
    # several of each point's nearest neighbors are equidistant from one another (up to precision).
    #
    # TODO: Improve this.
    unclipped_bg_domain_centroid_coords = np.stack(
        (
            unclipped_bg_domain_gdf.centroid.get_coordinates()["x"].to_numpy(),
            unclipped_bg_domain_gdf.centroid.get_coordinates()["y"].to_numpy(),
        ),
        axis=-1
    )
    unclipped_bg_domain_centroid_coords += np.random.normal(
        0, 1, size=unclipped_bg_domain_centroid_coords.shape
    )

    tri = triangulate(unclipped_bg_domain_centroid_coords)

    (
        # (o,) the index of the simplex containing each observation
        unclipped_simplex_indices_for_each_observation,
        # (o, 3) the 3 barycentric coords for each observation
        barycentric_coords_for_each_observation
    ) = simplex_info_for_points(point_coords=observation_domain_centroid_coords, tri=tri)

    # shape: (o, 3)
    #
    # The indices in `unclipped_bg_domain_gdf` of the 3 vertices involved the simplex containing
    # each `o`.
    unclipped_bg_vertex_indices = np.array([
        tri.vertices[i] for i in unclipped_simplex_indices_for_each_observation
    ])

    # shape: (o, 3, 2[x,y])
    #
    # For each observation, the grid-x and grid-y coordinates of each of the 3 vertices in the
    # simplex that contains the observation.
    bg_grid_coords_of_vertices_for_each_observation = np.stack(
        (
            unclipped_bg_domain_gdf.reset_index().grid_x.to_numpy()[unclipped_bg_vertex_indices],
            unclipped_bg_domain_gdf.reset_index().grid_y.to_numpy()[unclipped_bg_vertex_indices],
        ),
        axis=-1
    )

    # Figure out which vertices are involved in a simplex that actually contains an observation.
    #
    # This is how we "clip" out the vertices that wouldn't do anything in an interpolation on the
    # observations.
    #
    # shape: (v, 2[x,y])   v = # clipped/final vertices
    #   A list of bg_grid coordinates. Each grid-cell listed is involved in the interpolation of
    #   some observation grid-cell
    #
    # shape: (o, 3)
    #   Each element is an index into `clipped_bg_grid_coords`. This maps each observation to the
    #   indices of its 3 vertices as they appear in the "final list"
    clipped_bg_grid_coords, clipped_bg_vertex_indices = np.unique(
        bg_grid_coords_of_vertices_for_each_observation.reshape((-1, 2)),
        axis=0,
        return_inverse=True,
    )
    clipped_bg_vertex_indices = clipped_bg_vertex_indices.reshape((-1, 3))

    o = observation_coords.shape[0]
    v = clipped_bg_grid_coords.shape[0]  # At last, we know the final number of vertices

    # Build the interpolation matrix
    # shape: (o, v)
    bg_barycentric_coords = np.zeros((o, v))
    for o in range(bg_barycentric_coords.shape[0]):
        bg_indices = clipped_bg_vertex_indices[o]
        bg_barycentric_coords[o, bg_indices] = barycentric_coords_for_each_observation[o]

    # gdf: grid_x, grid_y, geometry:Polygon
    clipped_bg_domain_gdf = unclipped_bg_domain_gdf.loc[clipped_bg_grid_coords.tolist()]

    # UTM coordinates of each vertex. Used for setting distance-based priors
    bg_utm_coords = np.stack(
        (
            clipped_bg_domain_gdf.centroid.get_coordinates().x.to_numpy(),
            clipped_bg_domain_gdf.centroid.get_coordinates().y.to_numpy(),
        ),
        axis=-1
    )

    # Grid-coordinates of each vertex.
    bg_grid_coords = np.stack(
        (
            clipped_bg_domain_gdf.reset_index().grid_x.to_numpy(),
            clipped_bg_domain_gdf.reset_index().grid_y.to_numpy(),
        ),
        axis=-1
    )

    return bg_barycentric_coords, bg_utm_coords, bg_grid_coords


def define_mesh_prior(
    vertex_utm_coords: npt.NDArray[np.float_],
    vertex_mus: npt.NDArray[np.float_],
    vertex_sigmas: npt.NDArray[np.float_],
    correlation_decay_rate: float,
) -> tuple[npt.NDArray[np.float_], npt.NDArray[np.float_]]:
    """
    Set up the multivariate normal (MVN) prior on the vertices of the mesh.

    An exponential correlation decay model is used, and `correlation_decay_rate` is the
    user-specified rate of correlation decay.

    A `MVN(u, Sigma)` parameterization is used, with a vector of mu's, and a matrix of covariances.
    Each entry in covariance matrix is as follows:

        Sigma(i, j) = sigma_i * sigma_j * exp(-correlation_decay_rate * D(i, j))

        D(i, j) = distance in meters between vertex `i` and vertex `j`

    How to set the correlation decay rate:

        Let's say we want correlation to decay by a factor of two after X meters. This is equivalent
        to saying that we want correlation to decline by 0.5, because at D(i, j) = 0, correlation is
        always equal to 1:

        > 0.5 = exp(-correlation_decay_rate * X)
        > ln(0.5) = -correlation_decay_rate * X
        > correlation_decay_rate = -ln(0.5) / X

        Say X is 400km. Then

        > correlation_decay_rate = -ln(0.5) / 400e3 = 0.000001732

        This is a sensitive parameter that encodes important assumptions about how much the
        background concentration of methane is likely to vary over space.

    Inputs:
        vertex_utm_coords: (v,2) The UTM x/y-coordinates of each vertex
        vertex_mus: (v,) The prior mean at each vertex. Likely 0s.
        vertex_sigmas: (v,) The prior sigma at each vertex.
        correlation_decay_rate: The decay rate in an exponential covariance decay model

    Outputs:
        mus: (v,) vector of mean variation at each vertex
        Sigma: (v, v) covariance matrix for each pair of vertices
    """
    # shape (v, v)
    D = distance_matrix(x_coords=vertex_utm_coords[:, 0], y_coords=vertex_utm_coords[:, 1])

    correlation = np.exp(-correlation_decay_rate * D)  # shape (v, v)
    Sigma = np.outer(vertex_sigmas, vertex_sigmas) * correlation  # shape (v, v)

    return vertex_mus, Sigma
