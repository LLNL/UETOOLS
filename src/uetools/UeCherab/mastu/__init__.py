# Definitions specific to the MAST-U device

from . import wall
from . import bolometer


class MASTU:
    """

    Creates a Raysect World, adds a MAST-U wall
    """

    def __init__(self, cherab):
        self.cherab = cherab

        from raysect.optical import World

        self.world = World()

        self.wall = None
        self.emitter = None
        self._camera = None
        self._bolometer = None

    # Dictionary of available walls
    walls = {
        "wall1": wall.axisymmetric_wall1,
    }

    def add_wall(self, name="wall1"):
        """
        Add a wall to the Cherab world
        """
        from raysect.optical.material import AbsorbingSurface

        if not (name in self.walls):
            raise KeyError(
                f"Wall '{name}' not found. Available walls: {self.walls.keys()}"
            )

        self.wall = self.walls[name]()
        self.wall.parent = self.world
        self.wall.material = AbsorbingSurface()
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
            camera = mastu.camera
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
            self._camera.transform = translate(-3.5, -0.0, -0.3) * rotate_basis(
                Vector3D(1, 0, 0), Vector3D(0, 0, 1)
            )
            self._camera.pixel_samples = 1
        return self._camera

    @property
    def bolometer(self):
        """
        Bolometer arrays.
        """
        if self._bolometer is None:
            self._bolometer = bolometer.BolometerArrays(self.cherab, self.world)

        return self._bolometer
