# Routines to turn a UEDGE grid into triangles

import numpy as np
import matplotlib.pyplot as plt


class TriangularData:
    """
    Represents a set of triangles with data constant on them.
    Creates a Cherab Discrete2DMesh, and can then convert
    that to a 3D (axisymmetric) emitting material.
    """

    def __init__(self, vertices, triangles, data):
        self.vertices = vertices
        self.triangles = triangles
        self.data = data

        from raysect.core.math.function.float import Discrete2DMesh

        self.mesh = Discrete2DMesh(
            self.vertices, self.triangles, self.data, limit=False, default_value=0.0
        )

    def to_emitter(
        self,
        parent=None,
        cylinder_zmin=None,
        cylinder_zmax=None,
        cylinder_rmin=None,
        cylinder_rmax=None,
        step: float = 0.01,
    ):
        """
        Returns a 3D Cherab emitter, by rotating the 2D mesh about the Z axis

        step: Volume integration step length [m]
        
        """
        from raysect.core import translate
        from raysect.primitive import Cylinder, Subtract
        from raysect.optical.material import VolumeTransform
        from cherab.core.math import AxisymmetricMapper
        from cherab.tools.emitters import RadiationFunction

        if cylinder_zmin is None:
            cylinder_zmin = np.amin(self.vertices[:, 1])
        if cylinder_zmax is None:
            cylinder_zmax = np.amax(self.vertices[:, 1])
        if cylinder_rmin is None:
            cylinder_rmin = np.amin(self.vertices[:, 0])
        if cylinder_rmax is None:
            cylinder_rmax = np.amax(self.vertices[:, 0])

        rad_function_3d = AxisymmetricMapper(self.mesh)

        shift = translate(0, 0, cylinder_zmin)
        emitting_material = VolumeTransform(
            RadiationFunction(rad_function_3d, step=step), shift.inverse()
        )

        # Create an annulus by removing the middle from the cylinder.
        return Subtract(
            Cylinder(cylinder_rmax, cylinder_zmax - cylinder_zmin),
            Cylinder(cylinder_rmin, cylinder_zmax - cylinder_zmin),
            transform=shift,
            parent=parent,
            material=emitting_material,
        )

    def plot_2d(self, ax=None, nr: int = 150, nz: int = 150):
        """
        Make a 2D plot of the data

        nr, nz - Number of samples in R and Z
        """
        if ax is None:
            fig, ax = plt.subplots()

        Rmin, Zmin = np.amin(self.vertices, axis=0)
        Rmax, Zmax = np.amax(self.vertices, axis=0)

        from cherab.core.math import sample2d

        r, z, emiss_sampled = sample2d(self.mesh, (Rmin, Rmax, nr), (Zmin, Zmax, nz))

        image = ax.imshow(
            emiss_sampled.T, origin="lower", extent=(r.min(), r.max(), z.min(), z.max())
        )
        fig.colorbar(image)
        ax.set_xlabel("r")
        ax.set_ylabel("z")

        return ax


class Triangulate:
    """
    Computes a triangulation based on UEDGE R and Z coordinates
    """

    def __init__(self, rm, zm):
        """
        rm : [nx, ny, 5]
        zm : [nx, ny, 5]
        """

        # Check that the arrays have the right shape
        assert zm.shape == rm.shape
        assert len(rm.shape) == 3
        nx, ny, n = rm.shape
        assert n == 5

        # Build a list of vertices, and a list of triangles
        vertices = []
        triangles = []

        def vertex_index(R, Z):
            """
            Return the vertex index at given (R,Z) location.
            Note: This is inefficient linear search
            """
            for i, v in enumerate(vertices):
                vr, vz = v
                d2 = (vr - R) ** 2 + (vz - Z) ** 2
                if d2 < 1e-10:
                    return i
            # Not found
            vertices.append((R, Z))
            return len(vertices) - 1

        for ix in range(nx):
            for jy in range(ny):
                # Adding cell (ix, jy)
                # Get the vertex indices of the 4 corners
                vertex_inds = [
                    vertex_index(rm[ix, jy, n], zm[ix, jy, n]) for n in range(1, 5)
                ]
                # Choose corners so triangles have the same sign
                triangles.append(vertex_inds[0:3])  # Corners 1,2,3
                triangles.append(vertex_inds[:0:-1])  # Corners 4,3,2

        self.nx = nx
        self.ny = ny
        self.rm = np.copy(rm)
        self.zm = np.copy(zm)
        self.vertices = np.array(vertices)
        self.triangles = np.array(triangles)

    def plot_triangles(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()

        rs = self.vertices[self.triangles, 0]
        zs = self.vertices[self.triangles, 1]

        # Close the triangles
        rs = np.concatenate((rs, rs[:, 0:1]), axis=1)
        zs = np.concatenate((zs, zs[:, 0:1]), axis=1)

        ax.plot(rs.T, zs.T, "k")

        return ax

    def with_data(self, data2d):
        """
        Returns a new object containing vertices, triangles, and data
        """
        assert data2d.shape == (self.nx, self.ny)
        # Note: Two triangles per quad, so repeat data twice
        return TriangularData(
            self.vertices, self.triangles, np.repeat(data2d.flatten(), 2)
        )

    def to_voxel_grid(self, parent=None):
        """
        Returns a Cherab ToroidalVoxelGrid
        https://www.cherab.info/tools/tomography.html#cherab.tools.inversions.voxels.ToroidalVoxelGrid

        This can be used to calculate diagnostic geometry matrices.

        # Examples

        Triangulate the UEDGE mesh and turn it into voxels
        >>> voxels = c.cherab.triangulation.to_voxel_grid()

        Plot the voxels in the R-Z plane
        >>> voxels.plot()

        >>> experiment = c.cherab.d3d(zshift = 1.6).add_wall('wall1_bol')

        Put the voxels in the same world as the bolometers
        >>> voxels.parent = experiment.world

        for foil in experiment.bolometer.channels:
            print(f"{foil} : {np.amax(foil.calculate_sensitivity(voxels))}")

        """

        from cherab.tools.inversions.voxels import ToroidalVoxelGrid

        # Get arrays rs[ntriangle, 3] and zs[ntriangle, 3]
        rs = self.vertices[self.triangles, 0]
        zs = self.vertices[self.triangles, 1]

        # [ntriangles, 3, 2] array of coordinates
        voxel_data = np.stack((rs, zs), axis=-1)
        return ToroidalVoxelGrid(voxel_data, parent=parent)
