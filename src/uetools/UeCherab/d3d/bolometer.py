"""
48 channel bolometer system, in 4 arrays (A,B,C,D)

"""

import numpy as np

# Definitions of the bolometer arrays, their apertures and channels
#
# This was extracted from ExistingBoloGeometry.xlsx
#
# Note: X, Y, Z are (Radial, Height above midplane, Toroidal)
#       which is not the same as Cherab's
#
# Detector area 0.1208 cm^2, roughly square
D3D_BOLOMETER_ARRAYS = {
    "A": {
        # Box size [Width (cm), Height (cm), Depth (cm)]
        # Note: A and B arrays are stacked next to each other on the Y direction
        "Box": [15, 2, 15],
        # Aperture size [Tor (cm), Pol (cm)]
        "Aper": [0.6985, 0.70104],
        # Radius of curvature of foils [cm]
        "foil_curvature_radius": 0.0,
        # Size of the foil [cm]
        "foil_width": 0.3476,
        "foil_length": 0.3476,
        # Locations are [X (cm), Y (cm), Z (cm), angle (degrees)]
        "Aper Loc": [234.9529342, 72.99039831, 2.584120859, 236.2962527],
        "channels": {
            "#1": [235.0679877, 83.90805959, 2.12571092, 269.3249741],
            "#2": [235.7820876, 83.87845364, 2.149504633, 265.5673286],
            "#3": [236.4936999, 83.80071492, 2.176761665, 261.8109311],
            "#4": [237.1977394, 83.67737094, 2.207142263, 258.0558907],
            "#5": [237.891665, 83.50841888, 2.240569707, 254.302293],
            "#6": [239.3909403, 85.80259035, 2.117292732, 250.8002907],
            "#7": [240.0544622, 85.55254349, 2.155761723, 247.7996396],
            "#8": [240.7053121, 85.27208468, 2.196082936, 244.7999791],
            "#9": [241.3409573, 84.95361162, 2.238738604, 241.8013051],
            "#10": [241.9613894, 84.60472378, 2.283169775, 238.8036025],
            "#11": [242.5589848, 84.22541271, 2.329146293, 235.8068457],
            "#12": [243.1388316, 83.81061773, 2.377194231, 232.8109987],
            "#13": [239.0367475, 77.53903349, 2.714090652, 227.9395759],
            "#14": [239.5430048, 77.03043002, 2.766807512, 221.2047288],
            "#15": [239.9832541, 76.46602385, 2.821628682, 214.4729266],
        },
    },
    "B": {
        # Box size [Width (cm), Height (cm), Depth (cm)]
        "Box": [15, 2, 15],
        # Aperture size [Tor (cm), Pol (cm)]
        "Aper": [0.99314, 0.99314],
        # Radius of curvature of foils [cm]
        "foil_curvature_radius": 0.0,
        # Size of the foil [cm]
        "foil_width": 0.3476,
        "foil_length": 0.3476,
        # Locations are [X (cm), Y (cm), Z (cm), angle (degrees)]
        "Aper Loc": [231.1358145, 82.26133448, -2.498925026, 180.0636009],
        "channels": {
            "#16": [238.1197237, 86.00674368, -2.90243385, 208.2117133],
            "#17": [238.5323865, 85.11130651, -2.824235269, 201.1040492],
            "#18": [238.8078774, 84.1701374, -2.746825767, 193.9961377],
            "#19": [241.2759586, 83.61815856, -2.63171482, 187.6944178],
            "#20": [241.3609019, 82.64138187, -2.557450757, 182.2086101],
            "#21": [241.3518279, 81.65948789, -2.485652642, 176.7206564],
            "#22": [241.2487253, 80.68262533, -2.417065747, 171.2301139],
            "#23": [241.0515829, 79.72094292, -2.352435344, 165.7366337],
            "#24": [240.7654717, 78.78459172, -2.292353266, 160.2399739],
        },
    },
    "D": {
        # Box size [Width (cm), Height (cm), Depth (cm)]
        "Box": [15, 2, 15],
        # Aperture size [Tor (cm), Pol (cm)]
        "Aper": [0.762, 0.762],
        # Radius of curvature of foils [cm]
        "foil_curvature_radius": 0.0,
        # Size of the foil [cm]
        "foil_width": 0.3476,
        "foil_length": 0.3476,
        # Locations are [X (cm), Y (cm), Z (cm), angle (degrees)]
        "Aper Loc": [231.89438, -77.2541, -2.1438, 192.29],
        "channels": {
            "#25": [242.29314, -70.32244, -2.4841, 213.69],
            "#26": [242.697, -70.97014, -2.4841, 210.19],
            "#27": [243.06022, -71.6407, -2.4841, 206.69],
            "#28": [243.3828, -72.33412, -2.4841, 203.18],
            "#29": [243.6622, -73.04532, -2.4841, 199.69],
            "#30": [237.97006, -75.68946, -2.4841, 194.44],
            "#31": [238.11484, -76.44384, -2.4841, 187.44],
            "#32": [238.16818, -77.20838, -2.4841, 180.44],
            "#33": [238.12754, -77.97292, -2.4841, 173.42],
            "#34": [237.99292, -78.7273, -2.4841, 166.42],
            "#35": [237.7694, -79.45882, -2.4841, 159.43],
        },
    },
    "C": {
        # Box size [Width (cm), Height (cm), Depth (cm)]
        "Box": [20, 2, 20],
        # Aperture size [Tor (cm), Pol (cm)]
        "Aper": [0.9144, 0.70866],
        # Radius of curvature of foils [cm]
        "foil_curvature_radius": 0.0,
        # Size of the foil [cm]
        "foil_width": 0.3476,
        "foil_length": 0.3476,
        # Locations are [X (cm), Y (cm), Z (cm), angle (degrees)]
        "Aper Loc": [234.93222, -66.88074, 2.1438, 123.09],
        "channels": {
            "#36": [240.4999, -69.37248, 2.4841, 155.96],
            "#37": [240.1697, -70.00748, 2.4841, 149.23],
            "#38": [239.76838, -70.5993, 2.4841, 142.42],
            "#39": [239.29848, -71.14032, 2.4841, 135.78],
            "#40": [243.55044, -77.32268, 2.4841, 129.6],
            "#41": [242.99164, -77.75956, 2.4841, 126.6],
            "#42": [242.41252, -78.16596, 2.4841, 123.6],
            "#43": [241.81054, -78.54188, 2.4841, 120.6],
            "#44": [241.19078, -78.88732, 2.4841, 117.6],
            "#45": [240.55578, -79.1972, 2.4841, 114.6],
            "#46": [239.903, -79.47406, 2.4841, 111.6],
            "#47": [239.23752, -79.7179, 2.4841, 108.6],
            "#48": [236.18698, -72.86752, 2.4841, 101.91],
        },
    },
}


def cm2m(arr):
    """
    Convert a list from cm to m
    """
    return [e * 1e-2 for e in arr]


def create_array(array, world):
    """
    Follows example:
    https://www.cherab.info/demonstrations/bolometry/observing_radiation_function.html#bolometer-observing-radiation

    Returns a BolometerCamera object, containing an array of foils
    defined in the input `array` structure.
    
    Notes:
    - The poloidal plane in Cherab is X-Z : The torus is rotated around the Z axis.
    - 
    """

    from raysect.core import Point3D, Vector3D, translate, rotate_y
    from raysect.optical.material import AbsorbingSurface
    from raysect.primitive import Box, Subtract

    from cherab.tools.observers.bolometry import (
        BolometerCamera,
        BolometerSlit,
        BolometerFoil,
    )

    # Convenient constants
    XAXIS = Vector3D(1, 0, 0)
    YAXIS = Vector3D(0, 1, 0)
    ZAXIS = Vector3D(0, 0, 1)
    ORIGIN = Point3D(0, 0, 0)

    # Extract box shape. Convert to meters.
    BOX_WIDTH, BOX_HEIGHT, BOX_DEPTH = cm2m(array["Box"])

    # Aperture / Slit size. Convert to meters.
    SLIT_WIDTH, SLIT_HEIGHT = cm2m(array["Aper"])

    # In its local coordinate system, the camera's slit is located at the
    # origin and the foils below the z=0 plane, looking up towards the slit.
    camera_box = Box(
        lower=Point3D(-BOX_WIDTH / 2, -BOX_HEIGHT / 2, -BOX_DEPTH),
        upper=Point3D(BOX_WIDTH / 2, BOX_HEIGHT / 2, 0),
    )
    # Hollow out the box
    inside_box = Box(
        lower=camera_box.lower + Vector3D(1e-5, 1e-5, 1e-5),
        upper=camera_box.upper - Vector3D(1e-5, 1e-5, 1e-5),
    )
    camera_box = Subtract(camera_box, inside_box)
    # The slit is a hole in the box
    aperture = Box(
        lower=Point3D(-SLIT_WIDTH / 2, -SLIT_HEIGHT / 2, -1e-4),
        upper=Point3D(SLIT_WIDTH / 2, SLIT_HEIGHT / 2, 1e-4),
    )
    camera_box = Subtract(camera_box, aperture)
    camera_box.material = AbsorbingSurface()

    # Instance of the bolometer camera
    bolometer_camera = BolometerCamera(camera_geometry=camera_box)
    # The bolometer slit in this instance just contains targeting information
    # for the ray tracing, since we have already given our camera a geometry
    # The slit is defined in the local coordinate system of the camera

    slit = BolometerSlit(
        slit_id="Aperture",
        centre_point=ORIGIN,
        basis_x=XAXIS,
        dx=SLIT_WIDTH,
        basis_y=YAXIS,
        dy=SLIT_HEIGHT,
        parent=bolometer_camera,
        csg_aperture=False,
    )

    # The bolometer positions and orientations are given in the local coordinate
    # system of the camera, just like the slit / aperture

    aper_pos = cm2m(array["Aper Loc"][:3])  # (x,y,z) in m
    aper_degrees = array["Aper Loc"][-1]  # Angle in degrees

    # Transformation from local (camera) coordinates to global coordinates.
    # Operations are performed right to left: First a rotation then translation
    # The aper_degrees angle is anticlockwise, starting horizontal radially outward.
    # The camera starts orientated upwards i.e. at 90 degrees.
    camera_to_global = translate(aper_pos[0], aper_pos[2], aper_pos[1]) * rotate_y(
        90 - aper_degrees
    )

    bolometer_camera.transform = camera_to_global
    bolometer_camera.parent = world

    # The foil locations and angles are in the global coordinates,
    # so we need to map into camera coordinates
    # (alternatively, we could set the foil parent to `world`)
    global_to_camera = camera_to_global.inverse()

    for name, location in array["channels"].items():
        channel_pos = cm2m(location[:3])  # (x,y,z) in m
        channel_degrees = location[-1]

        # Transform from channel to global coordinates
        channel_to_global = translate(
            channel_pos[0], channel_pos[2], channel_pos[1]
        ) * rotate_y(90 - channel_degrees)

        # Get the channel transform in camera coordinates
        channel_to_camera = global_to_camera * channel_to_global

        foil = BolometerFoil(
            detector_id=name,
            centre_point=ORIGIN.transform(channel_to_camera),
            basis_x=XAXIS.transform(channel_to_camera),
            dx=array["foil_width"] * 1e-2,
            basis_y=YAXIS.transform(channel_to_camera),
            dy=array["foil_length"] * 1e-2,
            slit=slit,
            parent=bolometer_camera,
            units="Power",
            accumulate=False,
            curvature_radius=array["foil_curvature_radius"],
        )
        bolometer_camera.add_foil_detector(foil)

    return bolometer_camera


class BolometerArrays:
    def __init__(self, cherab, world):
        """
        # Parameters

        cherab : uetools.UeCherab.Cherab
            Interface betweem Cherab and a Case object

        world : raysect.optical.World
            The root node of the scene graph

        """
        self.cherab = cherab
        self.world = world

        self.arrays = {
            name: create_array(array, world)
            for name, array in D3D_BOLOMETER_ARRAYS.items()
        }
        channels = []
        for name, array in self.arrays.items():
            for channel in array:
                channels.append(channel)
        self._channels = channels

    def plotLinesOfSight(self, channels=None, ax=None, legend=False, show=True):
        """
        Plot lines of sight for all channels

        channels - A list of channels to plot e.g ["#1", "#2"].
                   If None then all channels are plotted

        ax - Matplotlib axis [optional]. If not given then a new figure is created.
        show - If true, calls matplotlib.pyplot.show()
        
        Returns the matplotlib axis
        """
        import math
        from raysect.core import Point2D
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots()
            ax.set_xlabel("Major radius [m]")
            ax.set_ylabel("Height [m]")

        def _point3d_to_rz(point):
            return Point2D(math.hypot(point.x, point.y), point.z)

        # Add a floor
        from raysect.core import Point3D
        from raysect.optical.material import AbsorbingSurface
        from raysect.primitive import Box

        floor = Box(
            lower=Point3D(-10, -10, -1.5),
            upper=Point3D(10, 10, -1.5),
            parent=self.world,
            material=AbsorbingSurface(),
            name="Z=-1.5 plane",
        )
        ceiling = Box(
            lower=Point3D(-10, -10, 1.5),
            upper=Point3D(10, 10, 1.5),
            parent=self.world,
            material=AbsorbingSurface(),
            name="Z=+1.5 plane",
        )
        left = Box(
            lower=Point3D(0.5, -10, -10),
            upper=Point3D(0.5, 10, 10),
            parent=self.world,
            material=AbsorbingSurface(),
            name="X=+0.5 plane",
        )

        for foil in self.channels:
            if (channels is not None) and (foil.name not in channels):
                continue  # Skip this channel
            slit_centre = foil.slit.centre_point
            slit_centre_rz = _point3d_to_rz(slit_centre)
            ax.plot(slit_centre_rz[0], slit_centre_rz[1], "ko")
            origin, hit, _ = foil.trace_sightline()
            centre_rz = _point3d_to_rz(foil.centre_point)
            ax.plot(centre_rz[0], centre_rz[1], "kx")
            origin_rz = _point3d_to_rz(origin)
            hit_rz = _point3d_to_rz(hit)
            ax.plot([origin[0], hit[0]], [origin[2], hit[2]], label=foil.name)

            normal = foil.normal_vector
            ax.arrow(origin_rz[0], origin_rz[1], normal[0] * 0.01, normal[2] * 0.01)

            x_width = foil.x_width
            basis_x = foil.basis_x
            ax.plot(
                [
                    origin_rz[0] - 0.5 * x_width * basis_x.x,
                    origin_rz[0] + 0.5 * x_width * basis_x.x,
                ],
                [
                    origin_rz[1] - 0.5 * x_width * basis_x.z,
                    origin_rz[1] + 0.5 * x_width * basis_x.z,
                ],
                "k",
            )

        floor.parent = None  # Remove from the world
        ceiling.parent = None
        left.parent = None

        if legend:
            ax.legend()

        if show:
            plt.show()

        return ax

    @property
    def channels(self):
        """
        List of bolometer channels (BolometerFoil objects) from all arrays
        """
        return self._channels

    def power(self):
        """
        Returns three lists:
        - Name of each channel (e.g. "#1")
        - Power incident on detector [Watts]
        - Error in power [Watts]
        """
        names = []
        power = []
        error = []
        for foil in self.channels:
            foil.pixel_samples = 100000
            foil.units = "Power"
            foil.observe()
            names.append(foil.name)
            power.append(foil.pipelines[0].value.mean)
            error.append(foil.pipelines[0].value.error())
        return names, power, error

    def sensitivity(self, ray_count: int = 10000):
        """Calculate the sensitivity / geometry matrix of the bolometer
        system.  This can be slow to calculate. The result is cached
        and reused if the mesh triangulation is unchanged.

        Note: Implementation makes assumption that triangulation
              uses two triangles per UEDGE cell. Sums sensitivity
              of the triangles.

        # Parameters

        ray_count : int
            The number of rays to use for each foil

        # Returns

        A BolometerSensitivity object

        """
        from ..bolometer_sensitivity import BolometerSensitivity

        return BolometerSensitivity(
            self.cherab, self.world, self.channels, ray_count=ray_count
        )
