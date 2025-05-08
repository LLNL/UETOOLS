"""
# Calculations of bolometer sensitivity.

This can be used to quickly perform multiple bolometer calculations,
without re-doing raytracing. It can also be used to analyse the
null space of the measurements, and invert bolometer profiles.
"""

import numpy as np

import logging

logger = logging.getLogger(__name__)


class BolometerSensitivity:
    """
    Describes the sensitivity of bolometer measurements to
    volumetric radiated power.
    """

    def __init__(self, cherab, world, channels, ray_count: int = 10000):
        """

        # Parameters

        ray_count : int
            The number of rays to use for each foil
        """
        # Get the mesh triangulation
        triangulation = cherab.triangulation

        # Get the triangulation as a set of voxels
        voxels = triangulation.to_voxel_grid(parent=world)

        logger.info("Calculating bolometer sensitivities")
        matrix = []
        for foil in channels:
            logger.info(f"    Foil {foil}")
            # Compute a row in the sensitivity matrix
            triangle_sensitivity = foil.calculate_sensitivity(
                voxels, ray_count=ray_count
            )
            ntriangles = triangle_sensitivity.size

            if ntriangles % 2 != 0:
                raise ValueError("Expected an even number of triangles")

            # Sum pairs of triangles: Reshape then sum
            cell_sensitivity = np.sum(
                triangle_sensitivity.reshape((ntriangles // 2, 2)), axis=-1
            )
            matrix.append(cell_sensitivity)

        self.nx = triangulation.nx
        self.ny = triangulation.ny
        self._matrix = np.array(matrix)
        self._measurement_projection = None

    @property
    def matrix(self):
        """The NumPy array sensitivity matrix. This has shape MxN
        where M is the number of bolometers and N is the number of
        grid cells in the domain.

        """
        return self._matrix

    @property
    def measurement_projection(self):
        """An NxN matrix, where N is the number of cells in the
        grid. This operator projects a radiation pattern onto the
        measurement space of the bolometers i.e. Removes any
        linear components that can't be measured.

        A tutorial on these operators is here:
        https://pmc.ncbi.nlm.nih.gov/articles/PMC10560448/
        """
        if self._measurement_projection is not None:
            return self._measurement_projection

        # Perform Singular Value Decomposition
        matrix_svd = np.linalg.svd(self.matrix, full_matrices=False)

        # Find the significant singular values
        svs = np.where(matrix_svd[1] > 1e-6 * matrix_svd[1][0])

        # Get the corresponding vectors that span the measurement space
        u = matrix_svd[2][svs[0], :]

        # Projection operator is NxN
        p = u.T @ u

        # Cache to re-use next time
        self._measurement_projection = p
        return p

    def measure(self, power):
        """
        Calculate the bolometer measurements by multiplying the
        input power by the sensitivity matrix.

        This will be much faster than performing ray-tracing.

        # Returns

        A 1D array of power to each foil
        """
        if power.shape != (self.nx, self.ny):
            raise ValueError(
                f"power has shape {power.shape} but must have shape {(self.nx, self.ny)}"
            )

        u = power.flatten()

        return self.matrix @ u

    def project_nullspace(self, power):
        """
        Project a 2D (nx, ny) array of radiated power onto the
        null space of the bolometer system.

        Note: The returned array of radiated power will
        in general have both positive and negative entries.

        # Returns

        A 2D NumPy array (nx, ny)
        """

        if power.shape != (self.nx, self.ny):
            raise ValueError(
                f"power has shape {power.shape} but must have shape {(self.nx, self.ny)}"
            )

        # The measurement space projection operator
        pmeas = self.measurement_projection

        # Turn into a 1D vector
        u = power.flatten()

        # Subtract the measurable components from the input
        return (u - (pmeas @ u)).reshape(power.shape)
