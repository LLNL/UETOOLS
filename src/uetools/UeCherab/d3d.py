# Definitions specific to the DIII-D device

import numpy as np


class AxisymmetricWall:
    """
    Represents an axisymmetric wall in Cherab
    """

    def __init__(self, rlim, zlim):
        from cherab.tools.primitives import axisymmetric_mesh_from_polygon

        wall_polygon = np.stack((rlim, zlim), axis=1)
        self.wall_mesh = axisymmetric_mesh_from_polygon(wall_polygon)


def axisymmetric_wall1():
    """
    Create an axisymmetric wall with open geometry, flat floor top and bottom.
    """
    rlim = np.array(
        [
            1.01680005,
            1.01680005,
            1.01680005,
            1.01680005,
            1.15090001,
            1.70959997,
            1.90550005,
            2.14059997,
            2.20350003,
            2.32929993,
            2.39059997,
            2.3756001,
            2.3756001,
            2.36169982,
            2.33959985,
            2.34979987,
            2.35249996,
            2.3527,
            2.36029983,
            2.3664999,
            2.36689997,
            2.36879992,
            2.37079978,
            2.36879992,
            2.36689997,
            2.3664999,
            2.36029983,
            2.3527,
            2.35249996,
            2.34979987,
            2.33959985,
            2.36169982,
            2.3756001,
            2.3756001,
            2.39059997,
            2.3289001,
            2.20239997,
            2.13910007,
            1.83669996,
            1.82210004,
            1.70000005,
            1.56299996,
            1.39999998,
            1.27499998,
            1.15390003,
            1.01680005,
            1.01680005,
            1.01680005,
            1.01680005,
            1.01680005,
        ]
    )
    zlim = np.array(
        [
            0.0,
            0.39976001,
            0.41475999,
            1.21730006,
            1.35140002,
            1.35140002,
            1.19299996,
            1.00269997,
            0.85009998,
            0.5449,
            0.39219999,
            0.39219999,
            0.3655,
            0.35839999,
            0.3249,
            0.2189,
            0.20100001,
            0.2,
            0.15000001,
            0.1097,
            0.1,
            0.05,
            0.0,
            -0.05,
            -0.1,
            -0.1097,
            -0.15000001,
            -0.2,
            -0.20100001,
            -0.2189,
            -0.3249,
            -0.35839999,
            -0.3655,
            -0.39219999,
            -0.39219999,
            -0.54579997,
            -0.85290003,
            -1.00639999,
            -1.36590004,
            -1.36590004,
            -1.36590004,
            -1.36590004,
            -1.36590004,
            -1.36590004,
            -1.36590004,
            -1.22880006,
            -0.80000001,
            -0.41475999,
            -0.39976001,
            0.0,
        ]
    )

    return AxisymmetricWall(rlim, zlim)


def axisymmetric_wall2():
    """
    This wall was extracted from a GEQDSK file for shot 163518
    """
    rlim = np.array(
        [
            1.01600003,
            1.01600003,
            1.01600003,
            1.01600003,
            1.01600003,
            1.01600003,
            1.01600003,
            1.01600003,
            1.01600003,
            1.01600003,
            1.01600003,
            1.01199996,
            1.00100005,
            1.02900004,
            1.04200006,
            1.046,
            1.05599999,
            1.097,
            1.10800004,
            1.11600006,
            1.13399994,
            1.148,
            1.16199994,
            1.18099999,
            1.18200004,
            1.18499994,
            1.19000006,
            1.19500005,
            1.20099998,
            1.20899999,
            1.21500003,
            1.222,
            1.22800004,
            1.23399997,
            1.23899996,
            1.24199998,
            1.24800003,
            1.25800002,
            1.26300001,
            1.27999997,
            1.27999997,
            1.27999997,
            1.30999994,
            1.32799995,
            1.36099994,
            1.38,
            1.41900003,
            1.41900003,
            1.37199998,
            1.37199998,
            1.60800004,
            1.64699996,
            1.78499997,
            2.06999993,
            2.12800002,
            2.24499989,
            2.32299995,
            2.37700009,
            2.36249995,
            2.36443496,
            2.36496902,
            2.36532402,
            2.36512089,
            2.36484194,
            2.36433411,
            2.36249995,
            2.37700009,
            2.13400006,
            1.78600001,
            1.76800001,
            1.76800001,
            1.68200004,
            1.37199998,
            1.37199998,
            1.41999996,
            1.41999996,
            1.273,
            1.153,
            1.01600003,
            1.01600003,
            1.01600003,
            1.01600003,
            1.01600003,
            1.01600003,
            1.01600003,
            1.01600003,
        ]
    )

    zlim = np.array(
        [
            0.00000000e00,
            9.63999987e-01,
            9.67999995e-01,
            1.00100005e00,
            1.01900005e00,
            1.07700002e00,
            1.07000005e00,
            1.09599996e00,
            1.11300004e00,
            1.13800001e00,
            1.14699996e00,
            1.16499996e00,
            1.21700001e00,
            1.21700001e00,
            1.16240001e00,
            1.16237998e00,
            1.16260004e00,
            1.16450000e00,
            1.16594005e00,
            1.16591001e00,
            1.16895998e00,
            1.17174995e00,
            1.17556000e00,
            1.18299997e00,
            1.18350005e00,
            1.18499994e00,
            1.18799996e00,
            1.19099998e00,
            1.19599998e00,
            1.20200002e00,
            1.20799994e00,
            1.21399999e00,
            1.22099996e00,
            1.23099995e00,
            1.23800004e00,
            1.24399996e00,
            1.25399995e00,
            1.27800000e00,
            1.28999996e00,
            1.33099997e00,
            1.34700000e00,
            1.34800005e00,
            1.34800005e00,
            1.34800005e00,
            1.34800005e00,
            1.34800005e00,
            1.34800005e00,
            1.30999994e00,
            1.30999994e00,
            1.29200006e00,
            1.09500003e00,
            1.07700002e00,
            1.07700002e00,
            1.03999996e00,
            9.92999971e-01,
            7.08999991e-01,
            5.18999994e-01,
            3.88999999e-01,
            4.00000006e-01,
            2.22000003e-01,
            1.33000001e-01,
            4.39999998e-02,
            -4.39999998e-02,
            -1.33000001e-01,
            -2.22000003e-01,
            -4.00000006e-01,
            -3.88999999e-01,
            -9.72999990e-01,
            -1.17400002e00,
            -1.21099997e00,
            -1.25000000e00,
            -1.25000000e00,
            -1.25000000e00,
            -1.32900000e00,
            -1.32900000e00,
            -1.36300004e00,
            -1.36300004e00,
            -1.36300004e00,
            -1.22300005e00,
            -1.22300005e00,
            -8.29999983e-01,
            -8.00000012e-01,
            -4.14999992e-01,
            -4.00000006e-01,
            -1.00000005e-03,
            0.00000000e00,
        ]
    )

    return AxisymmetricWall(rlim, zlim)


class D3D:
    """

    Creates a Raysect World, adds a DIII-D wall
    """

    def __init__(self, cherab):
        self.cherab = cherab

        from raysect.optical import World

        self.world = World()

        self.wall = None
        self.emitter = None
        self._camera = None

    def add_wall(self):
        """
        Add a wall to the Cherab world
        """
        self.wall = axisymmetric_wall2()
        self.wall.parent = self.world
        return self

    def set_emission(self, prad):
        """
        Set the emission function

        prad - 2D [nx,ny] array of [W/m^3] emission
        """
        assert len(prad.shape) == 2
        self.emitter = self.cherab.make_emitter(prad, parent=self.world)
        return self

    @property
    def camera(self):
        """
        Return pinhole camera setup.

        To use this camera:
            camera = d3d.camera
            plt.ion()
            camera.observe()
            plt.ioff()
            plt.show()
        """
        if self._camera is None:
            from raysect.core import translate, Vector3D, rotate_basis
            from raysect.optical.observer import PinholeCamera, PowerPipeline2D

            self._camera = PinholeCamera(
                (256, 256), pipelines=[PowerPipeline2D()], parent=self.world
            )
            self._camera.transform = translate(-1.5, -1.5, 0.67) * rotate_basis(
                Vector3D(1, 0, 0), Vector3D(0, 0, 1)
            )
            self._camera.pixel_samples = 1
        return self._camera
