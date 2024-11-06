# Interface to Cherab

from .triangulate import Triangulate
import numpy as np

import logging

logger = logging.getLogger(__name__)


class Cherab:
    """# Interface to Cherab

    Cherab (https://www.cherab.info/) is a python library for forward
    modelling diagnostics based on spectroscopic plasma emission.
    It is based on the Raysect (http://www.raysect.org/) scientific
    ray-tracing framework.

    # Examples

    Load a UEDGE Case
    >>> from uetools import Case
    >>> c = Case('input.yaml')

    Create a description of the device geometry.
    zshift shifts the UEDGE grid in the Z direction.
    >>> experiment = c.cherab.d3d(zshift = 1.6).add_wall('wall1')

    # Check alignment by plotting the triangulated mesh and wall
    >>> ax = c.cherab.plot_triangles()
    >>> experiment.wall.plotRZ(ax=ax)

    # Plot bolometer lines of sight
    >>> experiment.bolometer.plotLinesOfSight(ax=ax)

    # Set the plasma radiation emission function
    >>> experiment.set_emission(c.get('prad'))

    # Calculate power incident on each bolometer foil
    >>> power = experiment.bolometer.power()

    # Plot bolometer power [Watts]
    >>> import matplotlib.pyplot as plt
    >>> plt.errorbar(power[0], power[1], yerr=power[2], marker='x')

    # Plot a subset of bolometer lines of sight
    experiment.bolometer.plotLinesOfSight(channels=['#7', '#11', '#12'], ax=ax, legend=True)

    """

    def __init__(self, case):
        """
        case - A UETOOLS Case object

        Note: Cherab is not imported until the
        triangulation is required.
        """
        self.case = case
        self._triangulation = None
        self.zshift = 0.0

    @property
    def triangulation(self):
        """
        Return the triangulation of the UEDGE mesh.
        (Re-)Compute if needed.
        """

        rm = self.case.get("rm")
        zm = self.case.get("zm") - self.zshift  # Shifted by UEDGE grid generator

        # Recalculate if the mesh has changed
        if (
            self._triangulation is None
            or (rm.shape != self._triangulation.rm.shape)
            or (zm.shape != self._triangulation.zm.shape)
            or (not np.allclose(rm, self._triangulation.rm))
            or (not np.allclose(zm, self._triangulation.zm))
        ):
            logger.info("Triangulating the UEDGE grid")
            self._triangulation = Triangulate(rm, zm)
        return self._triangulation

    def plot_triangles(self, ax=None):
        """
        Plot the triangulated mesh

        ax - Matplotlib axis.
             If not provided then a new figure is created
        """
        return self.triangulation.plot_triangles(ax=ax)

    def make_emitter(self, prad, parent=None, step=1e-3):
        """
        prad - Array of radiation emission [Wm^-3]
        """
        return self.triangulation.with_data(prad).to_emitter(parent=parent, step=step)

    def d3d(self, zshift: float = 0.0):
        """
        Create an object representing the DIII-D device

        zshift - The distance (in m) that the UEDGE grid is shifted upward.
                 The UEDGE grid generator shifts LSN grids by
                 com.zshift = com.zdim/2 – com.zmid
                 which for DIII-D LSN cases is 1.6m
        """
        from .d3d import D3D

        # Set the Z shift so that UEDGE is consistent with device coordinates.
        # This may cause a retriangulation of the mesh.
        self.zshift = zshift
        return D3D(self)
