# Interface to Cherab

from .triangulate import Triangulate


class Cherab:
    """
    # Interface to Cherab

    
    """

    def __init__(self, case):
        """
        case - A UETOOLS Case object

        Note: Cherab is not imported until the
        triangulation is required.
        """
        self.case = case
        self._triangulation = None

    @property
    def triangulation(self):
        """
        Return the triangulation of the UEDGE mesh.
        (Re-)Compute if needed.
        """

        rm = self.case.get("rm")
        zm = self.case.get("zm")

        # Recalculate if the mesh has changed
        if (
            self._triangulation is None
            or (rm.shape != self._triangulation.rm.shape)
            or (zm.shape != self._triangulation.zm.shape)
            or (np.allclose(rm, self._triangulation.rm))
            or (np.allclose(zm, self._triangulation.zm))
        ):
            print("Triangulating the UEDGE grid")
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

    def d3d(self):
        """
        Create an object representing the DIII-D device
        """
        from .d3d import D3D

        return D3D(self)
