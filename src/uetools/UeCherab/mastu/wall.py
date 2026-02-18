import numpy as np


class AxisymmetricWall:
    """
    Represents an axisymmetric wall in Cherab
    """

    def __init__(self, rlim, zlim):
        from cherab.tools.primitives import axisymmetric_mesh_from_polygon

        self.rlim = rlim
        self.zlim = zlim

        wall_polygon = np.stack((rlim, zlim), axis=1)
        self.mesh = axisymmetric_mesh_from_polygon(wall_polygon)

    @property
    def parent(self):
        return self.mesh.parent

    @parent.setter
    def parent(self, value):
        self.mesh.parent = value

    @property
    def material(self):
        return self.mesh.material

    @material.setter
    def material(self, value):
        self.mesh.material = value

    def plotRZ(self, ax=None, show=True):
        """
        Plot as a line in the R-Z plane

        ax - Matplotlib axis [optional]. If not given then a new figure is created.
        show - If true, calls matplotlib.pyplot.show()
        """
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots()
            ax.set_xlabel("Major radius [m]")
            ax.set_ylabel("Height [m]")

        ax.plot(self.rlim, self.zlim, "k", label="Wall")

        if show:
            plt.show()

        return ax


def axisymmetric_wall1():
    """
    Create an axisymmetric wall with open geometry, flat floor top and bottom.
    """
    rlim = np.array(
        [
            1.45,
            1.45,
            1.45,
            1.3214,
            1.1904,
            0.89296,
            0.86938,
            0.83981,
            0.82229,
            0.81974,
            0.81974,
            0.82734,
            0.8548,
            0.89017,
            0.91974,
            0.94066,
            1.555,
            1.7301,
            1.35,
            1.09,
            1.09,
            0.90576,
            0.53594,
            0.5074,
            0.4788,
            0.333,
            0.334,
            0.261,
            0.261,
            0.261,
            0.261,
            0.261,
            0.261,
            0.261,
            0.261,
            0.261,
            0.261,
            0.261,
            0.261,
            0.334,
            0.333,
            0.333,
            0.4788,
            0.5074,
            0.53594,
            0.90576,
            1.09,
            1.09,
            1.35,
            1.7301,
            1.555,
            0.94066,
            0.91974,
            0.89017,
            0.8548,
            0.82734,
            0.81974,
            0.81974,
            0.82229,
            0.83981,
            0.86938,
            0.89296,
            1.1904,
            1.3214,
            1.45,
            1.45,
            1.45,
        ]
    )
    zlim = np.array(
        [
            0.0,
            0.5,
            0.82,
            0.82,
            1.007,
            1.304,
            1.3312,
            1.3826,
            1.4451,
            1.4812,
            1.4936,
            1.5318,
            1.5696,
            1.5891,
            1.5936,
            1.5936,
            1.567,
            1.68,
            2.06,
            2.06,
            2.06,
            1.8786,
            1.5017,
            1.4738,
            1.4458,
            1.303,
            1.1,
            0.502,
            0.348,
            0.348,
            0.146,
            0.146,
            0.0,
            -0.0,
            -0.146,
            -0.146,
            -0.348,
            -0.348,
            -0.502,
            -1.1,
            -1.1,
            -1.303,
            -1.4458,
            -1.4738,
            -1.5017,
            -1.8786,
            -2.06,
            -2.06,
            -2.06,
            -1.68,
            -1.567,
            -1.5936,
            -1.5936,
            -1.5891,
            -1.5696,
            -1.5318,
            -1.4936,
            -1.4812,
            -1.4451,
            -1.3826,
            -1.3312,
            -1.304,
            -1.007,
            -0.82,
            -0.82,
            -0.5,
            0.0,
        ]
    )

    return AxisymmetricWall(rlim, zlim)
