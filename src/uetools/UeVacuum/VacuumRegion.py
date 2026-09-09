"""
Vacuum Neutral Model (VNM) for UEDGE.

This module implements a vacuum neutral model for UEDGE simulations using the
Transport Matrix Method (TMM) to calculate particle transport between surfaces
in vacuum regions. The model handles neutral particle recycling, pumping, and
puffing in both main scrape-off layer (SOL) and private flux region (PFR) geometries.

Classes
-------
VacuumTests
    Test and visualization functions for vacuum region geometry.
VNM_interface
    Interface between UEDGE and vacuum region models.
VacuumRegion
    Core vacuum region implementation using TMM.
Surface
    Individual surface element in a vacuum region.

VNM Input Block (YAML Configuration)
-------------------------------------
The VNM setup is configured through the 'vnm' block in UEDGE YAML input files.

Structure:
    vnm:
        plot: bool, optional
            Generate geometry plots during setup (default: False).
        maxlength: float, optional
            Maximum surface segment length in meters (default: 0.01).
        regions: list of dict
            List of vacuum regions to configure. Each region has:

Region Configuration
~~~~~~~~~~~~~~~~~~~~
Each region in the 'regions' list must contain:

Required fields:
    name: str
        Unique identifier for the region.
    location: str
        Region location, either 'inner' (PFR) or 'outer' (SOL).
    isvacuummodel: int or dict
        Species indices to apply VNM. Can be:
        - int: Single species index
        - dict: {species_index: 1/0} to enable/disable per species

Mode-specific fields:
    mode: str, optional
        Operating mode (default: 'generate'):
        - 'generate': Create new vacuum region from UEDGE geometry or custom nodes
        - 'restore_surfaces': Load pre-computed surfaces from pickle, then compute
          transport matrices (allows adding pump/puff to existing geometry)
        - 'restore_boundary': Load complete pre-computed boundary conditions from
          HDF5 (no new calculations, fastest option)

    For mode='generate':
        nodes: str or array, optional
            Custom boundary nodes. Can be:
            - str: Path to .npy file or text file with (R,Z) coordinates
            - array: NumPy array of (R,Z) node pairs
            If not provided, automatically extracted from UEDGE geometry.
        material_recycling: float, optional
            Default recycling coefficient for material surfaces (default: 1.0).
        r_offset_plasma: float, optional
            Distribution offset for plasma surfaces (0=uniform, 1=cosine, default: 1).
        r_offset_material: float, optional
            Distribution offset for material surfaces (default: 1).
        reflections: int, optional
            Number of particle reflections to consider (default: 1e6).
        multiprocess: bool, optional
            Use multiprocessing for surface calculations (default: True).
        ncores: int, optional
            Number of CPU cores for multiprocessing (default: all available).
        write: bool, optional
            Save surfaces to pickle file after creation (default: False).
        savename: str, optional
            Filename for pickle output if write=True (required when write=True).

    For mode='restore_surfaces':
        surface_file: str, required
            Path to pickle (.pkl) or HDF5 (.hdf5) file containing pre-computed
            VacuumRegion surfaces and geometry.

        All parameters from 'generate' mode can also be used:
        material_recycling, r_offset_plasma, r_offset_material, reflections,
        multiprocess, ncores, write, savename

        Note: Restoring surfaces skips the expensive surface coupling calculation
        but still computes transport matrices. Use this to add pump/puff sources
        to existing geometry without recalculating view factors.

    For mode='restore_boundary':
        boundary_file: str, optional
            Path to HDF5 file containing saved boundary conditions (telematrices,
            puff arrays, VNM flags). If not provided, uses current UEDGE save file.

        Note: This mode performs no calculations - only loads pre-computed data.
        pump and puff configurations are ignored in this mode.

Pumping Configuration
~~~~~~~~~~~~~~~~~~~~~
    pump: list of dict, optional
        Pumping regions that reduce recycling. Each pump entry contains:
        - type: str, required
            Pump type. Currently supported: 'region'
        - nodes: str or array, required
            Polygon nodes defining pump region:
            - str: Path to .npy file or text file with (R,Z) coordinates
            - array: NumPy array of (R,Z) node pairs (minimum 3 nodes)
        - recycling: float, required
            Recycling coefficient for surfaces intersecting pump region
            (0 = perfect pump, 1 = no pumping).

        Available for: mode='generate', mode='restore_surfaces'
        Not available for: mode='restore_boundary'

Puffing Configuration
~~~~~~~~~~~~~~~~~~~~~
    puff: list of dict, optional
        Gas puffing sources. Each puff entry contains:
        - type: str, required
            Puff type. Currently supported: 'point'
        - location: tuple, required
            (R, Z) coordinates of puff location [meters].
            Puff automatically assigned to nearest material surface.
        - current: float, required
            Puff rate in particles/second.
        - igsp: int, required
            Gas species index (0-based) to puff.

        Available for: mode='generate', mode='restore_surfaces'
        Not available for: mode='restore_boundary'

Example Configurations
~~~~~~~~~~~~~~~~~~~~~~
1. Generate new SOL region with puffing:
    vnm:
        maxlength: 0.01
        regions:
            - name: 'sol'
              mode: 'generate'
              location: 'outer'
              isvacuummodel: {0: 1}
              material_recycling: 1.0
              puff:
                - type: "point"
                  location: [2.0, 2.7]
                  current: 1.e20
                  igsp: 0

2. Restore surfaces and add pumping (efficient for geometry reuse):
    vnm:
        regions:
            - name: 'pfr'
              mode: 'restore_surfaces'
              surface_file: 'pfr_vacuum.pkl'  # Pre-computed surfaces
              location: 'inner'
              isvacuummodel: {0: 1}
              pump:
                - type: 'region'
                  nodes: 'pump_geometry.npy'
                  recycling: 0.0  # New pump added to restored geometry

3. Restore from previous simulation:
    vnm:
        regions:
            - name: 'sol'
              mode: 'restore_boundary'
              boundary_file: 'previous_run.hdf5'
              location: 'outer'
              isvacuummodel: {0: 1}

4. Multi-species with custom nodes:
    vnm:
        maxlength: 0.005
        plot: True
        regions:
            - name: 'pfr'
              mode: 'generate'
              location: 'inner'
              nodes: 'custom_boundary.npy'
              isvacuummodel: {0: 1, 1: 1}  # Apply to species 0 and 1
              material_recycling: 0.95
              write: True
              savename: 'pfr_surfaces.pkl'  # Save for future reuse
              pump:
                - type: 'region'
                  nodes: [[1.4, 0.2], [1.5, 0.2], [1.5, 0.3], [1.4, 0.3]]
                  recycling: 0.1

5. Restore surfaces with modified pump/puff (saves computation time):
    vnm:
        regions:
            - name: 'sol'
              mode: 'restore_surfaces'
              surface_file: 'sol_surfaces.pkl'  # Reuse geometry
              location: 'outer'
              isvacuummodel: {0: 1}
              material_recycling: 0.98  # Changed from original
              reflections: 1e7  # Higher accuracy
              puff:
                - type: 'point'
                  location: [2.1, 2.8]  # Different puff location
                  current: 5.e19
                  igsp: 0
              pump:
                - type: 'region'
                  nodes: 'new_pump.npy'  # Add new pump
                  recycling: 0.05

Notes
-----
**Geometry and Implementation:**
- The VNM model is only implemented for single-null ('snull') and double-null
  ('dnull') geometries.
- The 'inner' location corresponds to the private flux region (PFR) and 'outer'
  to the scrape-off layer (SOL).
- Vacuum regions are automatically saved to HDF5 files when UEDGE cases are saved.

**Performance and Reusability:**
- mode='generate': Slowest, full computation (~seconds to minutes for complex geometry)
- mode='restore_surfaces': Medium, skips view factor calculation but computes
  transport matrices (~fraction of generate time). Allows adding/modifying pump/puff.
- mode='restore_boundary': Fastest, no computation (<1 second). Cannot modify
  pump/puff but useful for exact restart.
- Surface files (.pkl) can be reused across runs with the same geometry.
- Use write=True with mode='generate' to create reusable surface files.

**Pump and Puff Configuration:**
- Multiple pump and puff sources can be specified per region.
- When multiple pumps overlap, later pumps in the list take precedence.
- pump/puff available for mode='generate' and mode='restore_surfaces'.
- pump/puff are ignored for mode='restore_boundary' (uses saved values).

See Also
--------
VacuumRegion : Core vacuum region implementation
VNM_interface : Interface for UEDGE integration
"""


class VacuumTests:

    def twoSurfacePlot(self):
        """Plots a source and a receivng surface, including
        the flux triangle, normal vector, and distribution circle."""
        S1 = Surface((4, 2), (1, 6), 0)
        S2 = Surface((5, 9), (6, 8), 1)
        S1.showTwoSurfacePlot(S2, r_offset=1)

    def outerCirclePlot(self):
        """Plots the source surface and the outer circle
        used for obtaining the plot for comparison to
        analytic distributions."""
        S1 = Surface(
            (2, 5),
            (
                4,
                1,
            ),
            0,
        )
        S1.showOuterCirclePlot(r_offset=1)

    def analyticPlot(self, ax):
        """Generates the flux v. angle plot for analytic
        comparison of the distributions."""
        S1 = Surface((2, 5), (4, 1), 0)
        S1.showAnalyticPlot(ax, comparison=True, r_offset=1, showBothDist=True)

    def UanalyticPlot(self, ax):
        """Generates the flux v. angle plot for a uniform
        distribution."""
        S1 = Surface((2, 5), (4, 1), 0)
        S1.analyticUniform(ax)

    def trianglePlot(self):
        from shapely import Point, LineString
        import math
        import numpy as np

        """Plots a triangle geometry."""
        S1 = Surface((1, 3), (2, 6), 0)
        height = math.sqrt(3) * S1.surfaceLength / 2  # height of the et
        vertex = S1.normal.interpolate(height)
        test = VacuumRegion(
            [(S1.start.x, S1.start.y), (S1.end.x, S1.end.y), (vertex.x, vertex.y)],
            multiprocess=False,
        )
        f = test.plotGeometry(labels=True, showCircle=True)

    def squarePlot(self):
        from shapely import Point, LineString

        """Plots a square geometry."""
        S1 = Surface((1, 2), (3, 6), 0)
        # # # Make the sides of the square perpendicular to self # # #
        side1Start = (S1.end.x, S1.end.y)
        side1End = S1.normalHelper(
            abs(S1.dx), abs(S1.dy), side1Start[0], side1Start[1], False
        )

        side3End = (S1.start.x, S1.start.y)
        side3Start = S1.normalHelper(
            abs(S1.dx), abs(S1.dy), side3End[0], side3End[1], False
        )

        test = VacuumRegion(
            [side1Start, side1End, side3Start, side3End], multiprocess=False
        )
        f = test.plotGeometry(labels=True, showCircle=True)

    def shadedSquarePlot(self):
        from shapely import Point, LineString

        """Plots a shaded square geometry. 
            To be used alongside the source surface S1-- Surface((2, 1), (1, 1))-- 
            which is defined in the lineOfSightPlot function in the test functions."""

        geometryVertices = [
            (3, 1),
            (1, 1),
            (1, 5),
            (2, 6),
            (1, 7),
            (7, 7),
            (7, 1),
            (5, 1),
            (5, 3),
            (3, 3),
        ]

        test = VacuumRegion(geometryVertices, P=1, multiprocess=False)
        test.saveVacuumRegion("SavedVacuumRegion")
        # for i in test.errors:
        #     f = test.plotGeometry(labels=True, testsurf=i)
        #     test.surfaces[i].plotSelf(ax=f)
        #     test.surfaces[i].printReport()
        f = test.plotGeometry(labels=True, testsurf=6, showCircle=True)
        return test


class VNM_interface:
    """Interface between UEDGE Case and Vacuum Neutral Model.

    This class manages vacuum region setup, restoration, and integration with UEDGE.
    It handles multiple vacuum regions (SOL and PFR), supports various input formats,
    and coordinates between UEDGE variables and vacuum model outputs.

    Parameters
    ----------
    case : uetools.Case
        UEDGE case object with active simulation.
    vnm_setup : dict
        Configuration dictionary for vacuum regions. See module docstring for
        detailed structure. Must contain at least:
        - 'regions': list of region configurations

    Attributes
    ----------
    coupling : uetools.UeCase.UeCoupling
        Grid coupling utilities from UEDGE case.
    regions : dict
        Dictionary of VacuumRegion objects, keyed by region name.
    nodes : dict
        Boundary node arrays for 'inner' (PFR) and 'outer' (SOL) regions.
    output : dict
        Telematrices and puff vectors for each region:
        - 'telematrix': (nx+2, nx+2, 6) transport matrix array
        - 'puff': (nx+2, ngsp) puffing input array
    uevars : dict
        UEDGE variable names for 'inner' and 'outer' regions:
        - 'isvacuummodelpf'/'isvacuummodelw': VNM enable flags
        - 'cftelematrixpf'/'cftelematrixw': Transport matrices
        - 'cfteleoutpf'/'cfteleoutw': Output scaling factors
        - 'fngyi_use'/'fngyo_use': Puffing arrays

    Methods
    -------
    generate(regions, restore_surfaces=False, maxlength=0.01, plot=False)
        Generate or restore vacuum regions from geometry or files.
    restore(save_file, restore_vars=None, location='')
        Restore telematrices and boundary conditions from HDF5 file.
    plot_grid(sol_plot=True, pfr_plot=True, **kwargs)
        Plot vacuum region geometries.
    save_hdf5(file, **kwargs)
        Save vacuum setup to HDF5 file.

    Raises
    ------
    Exception
        If geometry is not 'snull' or 'dnull'.
    TypeError
        If vnm_setup is not a dictionary.
    KeyError
        If required configuration keys are missing.

    Examples
    --------
    Basic usage::

        from uetools import Case
        vnm_config = {
            'regions': [{
                'name': 'sol',
                'location': 'outer',
                'mode': 'generate',
                'isvacuummodel': {0: 1}
            }]
        }
        case = Case('input.yaml')
        vnm = VNM_interface(case, vnm_config)

    Notes
    -----
    - Automatically populates UEDGE variables after setup
    - Supports three modes: 'generate', 'restore_surfaces', 'restore_boundary'
    - Can handle multiple regions simultaneously
    - Telematrices are automatically expanded/padded for PFR geometry
    """

    def __init__(self, case, vnm_setup):
        from numpy import load, array
        from collections import defaultdict

        if case.get("geometry")[0].strip().decode("UTF-8") not in ["snull", "dnull"]:
            raise Exception(
                "VNM model only implemented for single-null and double-null geometries!"
            )
        self.coupling = case.coupling
        self.set = case.setue
        self.populate = case.populate
        self.info = case.info
        self.tools = case.tools
        self.getue = case.getue
        self.hdf5search = case.tools.hdf5search
        self.uevars = {
            "inner": [
                "isvacuummodelpf",
                "cfteleoutpf",
                "cftelematrixpf",
                "fngyi_use",
            ],
            "outer": [
                "isvacuummodelw",
                "cfteleoutw",
                "cftelematrixw",
                "fngyo_use",
            ],
        }

        # TODO: Test HDF5 restore with UETOOLS HDF5 saves

        def read_txt(file):
            nodes = []
            with open(file) as f:
                for line in f:
                    nodes.append(
                        tuple(
                            [
                                float(x)
                                for x in line.replace(",", " ").strip().split(" ")
                            ]
                        )
                    )
            return nodes

        # Parse VNM setup - determine whether to pass to restore or generate
        if vnm_setup is not None:
            if not isinstance(vnm_setup, dict):
                raise TypeError("vnm_setup must be dict")
            if "regions" not in vnm_setup:
                raise Exception("At least one region must be defined")
            # Check that regions are present and satifsy requirements
            for region in vnm_setup["regions"]:
                if "location" not in region:
                    raise Exception(
                        "Region must have 'location' set to 'inner'/'outer'"
                    )
                else:
                    if region["location"].lower() not in ["inner", "outer"]:
                        raise Exception("Region 'location' must be 'inner'/'outer'")
                if "isvacuummodel" not in region:
                    raise Exception("Specify VNM model for species using isvacuummodel")
            regions = vnm_setup.pop("regions")
            # Detect whether to generate, restore surfaces, or restore BCs
            groups = defaultdict(list)
            for region in regions:
                mode = region.get("mode", "generate")
                groups[mode].append(region)
            generate = groups["generate"]
            restore_surfaces = groups["restore_surfaces"]
            restore_boundary = groups["restore_boundary"]
            for region in restore_surfaces:
                if "surface_file" not in region:
                    raise KeyError(
                        f"'surface_file' to restore not specified for region {region['name']}"
                    )
            for region in restore_boundary:
                if "boundary_file" not in region:

                    print(
                        f"'boundary_file' to restore not specified for region {region['name']}, using save file '{self.info['savefile']}'"
                    )
                    region["boundary_file"] = self.info["savefile"]
                # Restore regions from save file as requested
                self.restore(
                    region["boundary_file"],
                    self.uevars[region["location"]],
                    region["name"],
                )
            # Check generated regions satisfy conditions
            for region in generate + restore_surfaces:
                if "nodes" in region:
                    if "savefile" in region:
                        raise Exception(
                            "Either specify save file or node list"
                            + f" for {region['name']}, not both!"
                        )
                    if isinstance(region["nodes"], str):
                        try:
                            region["nodes"] = load(region["nodes"])
                        except:
                            region["nodes"] = array(read_txt(region["nodes"]))
                if "pump" in region:
                    for pump in region["pump"]:
                        if isinstance("nodes", str):
                            try:
                                pump["nodes"] = load(pump["nodes"])
                            except:
                                pump["nodes"] = array(read_txt(pump["nodes"]))

        regions = (
            [x.copy() for x in generate]
            + [x.copy() for x in restore_boundary]
            + [x.copy() for x in restore_surfaces]
        )
        if len(generate) > 0:
            self.generate(generate, **vnm_setup)
        if len(restore_surfaces) > 0:
            self.generate(restore_surfaces, restore_surfaces=True, **vnm_setup)
        vnm_setup["regions"] = regions
        self.populate()

    def restore(self, save_file, restore_vars=None, location=""):
        """Restore telematrices and boundary conditions from HDF5 save file.

        Loads pre-computed vacuum region data (telematrices, puff arrays, VNM flags)
        from an HDF5 file and populates UEDGE variables. Useful for restarting
        simulations with existing VNM setup.

        Parameters
        ----------
        save_file : str or None
            Path to HDF5 file containing VNM data. If None, uses current case
            save file (self.info['savefile']).
        restore_vars : list of str, optional
            UEDGE variable names to restore. If None, restores all standard
            VNM variables for both inner and outer regions:
            - isvacuummodelpf, isvacuummodelw
            - cftelematrixpf, cftelematrixw
            - cfteleoutpf, cfteleoutw
            - fngyi_use, fngyo_use
        location : str, optional
            Region identifier string for print message (default: '').

        Returns
        -------
        None

        Notes
        -----
        - Variables not found in save file are skipped (assumes defaults)
        - Automatically calls self.populate() after restoration
        - Does not restore Surface objects or pump/puff configurations
        - Only restores boundary condition data needed for UEDGE execution

        Examples
        --------
        Restore from a specific file::

            vnm.restore('previous_run.hdf5', location='SOL')
            # Output: Successfully restored SOL VNM from previous_run.hdf5

        Use current save file::

            vnm.restore(None)
        """
        import h5py
        import warnings
        import numpy
        from uedge import com, bbb

        if restore_vars is None:
            restore_vars = self.uevars["inner"] + self.uevars["outer"]
        if save_file is None:
            save_file = self.info["savefile"]

        for var in restore_vars:
            val = self.hdf5search(save_file, var)
            if val is None:
                # Assume defaults used and not written to save
                pass
            else:
                self.set(var, val)
        self.populate()
        print(f"Successfully restored {location} VNM from {save_file}")
        return

    def generate(self, regions, restore_surfaces=False, maxlength=0.01, plot=False):
        """Generate vacuum regions and populate UEDGE transport matrices.

        Creates VacuumRegion objects from UEDGE geometry or user-provided nodes,
        computes transport matrices via TMM, and populates UEDGE boundary condition
        arrays. Handles automatic grid extraction, matrix padding for PFR geometry,
        and configuration of multiple regions simultaneously.

        Parameters
        ----------
        regions : list of dict
            Region configurations. Each dict must contain:
            - 'name': Region identifier
            - 'location': 'inner' (PFR) or 'outer' (SOL)
            - 'isvacuummodel': Species control (int or dict)
            - 'nodes': (optional) Custom boundary nodes
            - 'pump': (optional) List of pump configurations
            - 'puff': (optional) List of puff configurations
            - 'material_recycling': (optional) Default recycling coefficient
            - 'surface_file': (required if restore_surfaces=True) Path to pickle
            Additional VacuumRegion parameters can be included.
        restore_surfaces : bool, optional
            If True, restore VacuumRegion from 'surface_file' instead of
            generating from scratch (default: False).
        maxlength : float, optional
            Maximum surface segment length [meters] for automatic grid
            discretization (default: 0.01). Smaller values increase accuracy
            but raise computational cost.
        plot : bool, optional
            Generate geometry plots after setup (default: False).

        Returns
        -------
        None

        Side Effects
        ------------
        - Creates self.regions dict with VacuumRegion objects
        - Populates self.output dict with telematrices and puff vectors
        - Sets UEDGE variables:
          * cftelematrix{pf,w}: Transport matrices (nx+2, nx+2, 6)
          * isvacuummodel{pf,w}: VNM enable flags
          * cfteleout{pf,w}: Output scaling (set to 1.0)
          * fngy{i,o}_use: Puffing arrays (nx+2, ngsp)
        - Calls self.populate() to update UEDGE state

        Notes
        -----
        **Automatic Grid Extraction:**
        - For 'outer' (SOL): Extracts outer wall and plate surfaces
        - For 'inner' (PFR): Extracts inner wall and plate surfaces
        - Uses self.coupling.get_snull_vacuum_regions(maxlength)
        - Can be overridden with custom 'nodes' in region config

        **PFR Matrix Padding:**
        Inner regions require special handling because PFR surfaces are
        non-contiguous in UEDGE's poloidal indexing (broken by core).
        Telematrix is expanded and zero-padded to match full (nx+2, nx+2)
        grid structure.

        **Species Control:**
        isvacuummodel can be:
        - int: Apply to all species
        - dict: {species_index: 1/0} for per-species control
        - list: Applied to sequential indices

        **Multiple Regions:**
        All regions are generated in a single call. Regions are independent
        except for shared UEDGE grid parameters (nx, ixpt1, ixpt2, ngsp).

        Examples
        --------
        Generate two regions with custom settings::

            regions = [
                {
                    'name': 'sol',
                    'location': 'outer',
                    'isvacuummodel': {0: 1},
                    'material_recycling': 1.0,
                    'puff': [{'type': 'point', 'location': (2, 2.7),
                              'current': 1e20, 'igsp': 0}]
                },
                {
                    'name': 'pfr',
                    'location': 'inner',
                    'isvacuummodel': {0: 1},
                    'pump': [{'type': 'region', 'nodes': pump_nodes,
                              'recycling': 0.1}]
                }
            ]
            vnm.generate(regions, maxlength=0.005, plot=True)

        Restore from pre-computed surfaces::

            regions = [{
                'name': 'sol',
                'location': 'outer',
                'surface_file': 'sol_surfaces.pkl',
                'isvacuummodel': {0: 1}
            }]
            vnm.generate(regions, restore_surfaces=True)

        Raises
        ------
        KeyError
            If required configuration keys are missing from region dict.

        See Also
        --------
        VacuumRegion : Core vacuum region implementation
        restore : Restore from complete HDF5 boundary conditions
        """

        import h5py
        from numpy import zeros, hstack, vstack, pad
        from uedge import com, bbb

        self.populate(verbose=False)
        self.regions = []

        self.nx = self.getue("nx")
        self.ixpt1 = self.getue("ixpt1")[0]
        self.ixpt2 = self.getue("ixpt2")[0]
        self.ngsp = self.getue("ngsp")
        self.nx = self.getue("nx")
        (sol_nodes, pfr_nodes) = self.coupling.get_snull_vacuum_regions(
            maxlength=maxlength
        )
        self.nodes = {"inner": pfr_nodes, "outer": sol_nodes}
        self.regions = {}
        self.output = {}
        dimension = self.nx + 2
        for region in regions:
            if "name" in region:
                name = region.pop("name")
            else:
                name = len(regions) + 1
            isvacuummodel = region.pop("isvacuummodel")
            # Perform inner/outer setup
            location = region.pop("location").lower()
            (vnm, P) = self.nodes[location]
            if location == "inner":
                P = self.ixpt1 + (self.nx - self.ixpt2)
                savekey = ("pf", "i")
            else:
                P = self.nx
                savekey = ("w", "o")
            self.output[name] = {
                "telematrix": zeros((dimension, dimension, 6)),
                "puff": zeros((dimension, self.ngsp)),
            }
            if restore_surfaces:
                print(f"Restoring vacuum region '{name}' from pickle/HDF5")
                vnm = region.pop("surface_file")
            elif "nodes" in region:
                print(f"Generating vacuum region '{name}' from user-defined nodes")
                vnm = region.pop("nodes")
            else:
                print(f"Generating new vacuum region '{name}'")
            self.regions[name] = VacuumRegion(vnm, P=P, **region)

            if location == "inner":
                # Expand PF matrix  along core cut: get dimensions and mismatch
                dim_expand = dimension - P - 2
                tele = pad(self.regions[name].telematrix.transpose(), pad_width=1)
                # Expand PF matrix along vertical axis
                tele = vstack(
                    [
                        tele[: self.ixpt1 + 1],
                        zeros((dim_expand, P + 2)),
                        tele[self.ixpt1 + 1 :],
                    ]
                )
                # Expand PF matrix along horizontal axis
                tele = hstack(
                    [
                        tele[:, : self.ixpt1 + 1],
                        zeros((dimension, dim_expand)),
                        tele[:, self.ixpt1 + 1 :],
                    ]
                )
                for j in range(6):
                    self.output[name]["telematrix"][:, :, j] = tele
                # Expand and pad PF puffing array
                self.output[name]["puff"] = vstack(
                    [
                        zeros((1, 6)),
                        self.regions[name].puff_vector[: self.ixpt1],
                        zeros((dim_expand, 6)),
                        self.regions[name].puff_vector[self.ixpt1 :],
                        zeros((1, 6)),
                    ]
                )
            else:
                # Populate local pump
                self.output[name]["puff"] = vstack(
                    [
                        zeros((1, 6)),
                        self.regions[name].puff_vector,
                        zeros((1, 6)),
                    ]
                )
                for j in range(6):
                    self.output[name]["telematrix"][:, :, j] = pad(
                        self.regions[name].telematrix.transpose(), pad_width=1
                    )
            # Populate UEDGE cftelematrixw array
            self.set(f"cftelematrix{savekey[0]}", self.output[name]["telematrix"])
            # Turn on the VNM model in UEDGE
            if isinstance(isvacuummodel, dict):
                for key, value in isvacuummodel.items():
                    if isinstance(value, int):
                        self.getue(f"isvacuummodel{savekey[0]}", cp=False)[key] = value
                    elif isinstance(value, list):
                        listlen = len(value)
                        self.getue(f"isvacuummodel{savekey[0]}", cp=False)[
                            key : key + listlen
                        ] = value
            else:
                self.set(f"isvacuummodel{savekey[0]}", isvacuummodel)
            self.set(f"cfteleout{savekey[0]}", 1.0)
            # Populate the puffing array
            self.set(f"fngy{savekey[1]}_use", self.output[name]["puff"][:, : self.ngsp])

        # Plot setup if requested
        if plot:
            self.plot_grid(pf_test_surf=[], label=False)

    def plot_grid(
        self,
        sol_plot=True,
        pfr_plot=True,
        sol_test_surf=[],
        pf_test_surf=[],
        label=False,
        **kwargs,
    ):
        from matplotlib.pyplot import subplots

        f, ax = subplots()
        for regionname, region in self.regions.items():
            region.plotGeometry(labels=label, ax=ax, **kwargs)

    def save_hdf5(self, file, **kwargs):
        """Saves VNM data to HDF5"""
        from h5py import File

        if not isinstance(file, File):
            raise TypeError("file must be an open HDF5 File object")
        vnm = file.require_group("vnm")
        if hasattr(self, "sol"):
            sol = vnm.require_group("sol")
            self.sol.writeVacuumSetup(sol, **kwargs)
        if hasattr(self, "pfr"):
            pfr = vnm.require_group("pfr")
            self.pfr.writeVacuumSetup(pfr, **kwargs)
        for var in self.uevars:
            try:
                val = self.getue(var)
            except:
                val = None
            if val is not None:
                if var in vnm:
                    del vnm[var]
                vnm.create_dataset(var, data=val)
        return


class VacuumRegion:
    """Vacuum region implementation using Transport Matrix Method (TMM).

    Represents a closed vacuum region composed of plasma-facing and material surfaces.
    Calculates particle transport between surfaces using view factors, reflection
    coefficients, and iterative transport matrices to model neutral particle behavior.

    The Transport Matrix Method computes the probability that a particle leaving
    surface i will eventually reach surface j after multiple reflections. This is
    encoded in the telematrix T[i,j], which gives the fraction of particles emitted
    from surface i that are absorbed by surface j.

    Parameters
    ----------
    nodeList : str, array-like, or None
        Boundary geometry specification:
        - str ending in '.hdf5': Restore from HDF5 file using hdf5location
        - str ending in '.pkl': Restore from pickle file
        - str (other): Path to text file with (R,Z) coordinates
        - array-like: List/array of (R,Z) node coordinate tuples
        - None: Not allowed, raises TypeError
    P : int, optional
        Number of plasma-facing surfaces (default: 0). These surfaces appear
        first in the node list and use plasma distribution settings.
    r_offset_plasma : float, optional
        Distribution offset for plasma surfaces (default: 1):
        - 0: Uniform angular distribution
        - 1: Cosine (Knudsen) distribution
        - (0,1): Intermediate distribution
    r_offset_material : float, optional
        Distribution offset for material surfaces (default: 1).
    multiprocess : bool, optional
        Use multiprocessing for surface coupling calculations (default: True).
        Significantly speeds up setup for large geometries.
    ncores : int, optional
        Number of CPU cores for multiprocessing (default: all available).
    verbose : bool, optional
        Print progress messages during setup (default: True).
    material_recycling : float, optional
        Default recycling coefficient for material surfaces (default: 1.0).
        Fraction of incident particles that are re-emitted.
    pump : list of dict, optional
        Pumping region configurations (default: None). Each entry must have:
        - 'type': 'region' (only supported type)
        - 'nodes': Array of (R,Z) coordinates defining pump polygon (≥3 nodes)
        - 'recycling': Recycling coefficient for surfaces in pump region
        Later pumps override earlier ones for overlapping surfaces.
    puff : list of dict, optional
        Gas puffing source configurations (default: None). Each entry must have:
        - 'type': 'point' (only supported type)
        - 'location': (R, Z) tuple for puff location
        - 'current': Particle injection rate [particles/s]
        - 'igsp': Gas species index (0-based)
        Puff assigned to nearest material surface automatically.
    reflections : int, optional
        Number of particle reflections to compute (default: 1e6).
        Higher values more accurately capture transport in high-recycling regions
        but increase computation time as (AB)^reflections.
    hdf5location : str, optional
        Path within HDF5 file for restore (default: "vnm").
    savename : str, optional
        Filename for pickle output if write=True (default: None).
    isvacuummodel : int or dict, optional
        UEDGE species control (deprecated, handled by VNM_interface).
    write : bool, optional
        Save region to pickle file after creation (default: False).
        Requires savename to be specified.
    **kwargs
        Additional keyword arguments (currently unused).

    Attributes
    ----------
    surfaces : dict
        Dictionary of Surface objects, keyed by surface index.
    geometry : shapely.Polygon
        Polygon representation of vacuum region boundary.
    telematrix : ndarray, shape (P, P)
        Transport matrix T[i,j] giving fraction of particles from plasma
        surface i that reach plasma surface j after all reflections.
    puff_vector : ndarray, shape (P, 6)
        Puffing input array for each plasma surface and gas species.
    R_array : ndarray, shape (numSurfaces,)
        Recycling coefficient for each surface.
    C_array : ndarray, shape (numSurfaces, numSurfaces)
        View factor matrix C[i,j]: fraction of particles from i seeing j
        on first flight (no reflections).
    R_matrix : scipy.sparse.csr_array
        Sparse diagonal matrix of recycling coefficients.
    C_matrix : scipy.sparse.csr_array
        Sparse view factor matrix.
    A_matrix, B_matrix, AB_matrix : scipy.sparse.csr_array
        TMM intermediate matrices for transport calculation.
    AB_power_A : scipy.sparse.csr_array
        (AB)^reflections @ A, the final transport operator.
    pumping_regions : list of dict
        Created pump region information.
    puffs : list of dict
        Created puff source information with computed puff vectors.
    time : float
        Wall-clock time [seconds] for surface coupling calculation.
    numSurfaces : int
        Total number of surfaces (plasma + material).

    Methods
    -------
    matrices()
        Construct TMM matrices (R, C, A, B, AB).
    createTeleMatrix()
        Compute final transport matrix from (AB)^M @ A.
    create_pump_region(pump_setup)
        Create pumping region from configuration dict.
    create_puff(puff_setup)
        Create puff source from configuration dict.
    saveVacuumRegion(savename)
        Save region to pickle file.
    writeVacuumSetup(obj, write_matrices=True)
        Write region to HDF5 file or group.
    checkContinuity(verbose=True)
        Verify flux conservation (sum of view factors ≈ 1).
    plotGeometry(ax=None, labels=False, testsurf=[], **kwargs)
        Plot vacuum region geometry with surfaces and sources.
    heatmapPlot()
        Generate heatmap visualizations of transport matrices.
    matrixPower(matrix, power)
        Compute sparse matrix to integer power.

    Raises
    ------
    TypeError
        If nodeList is None or has invalid type.
    ValueError
        If write=True but savename not specified.
    KeyError
        If required pump/puff configuration keys missing.
    AttributeError
        If pump region has fewer than 3 nodes.

    Examples
    --------
    Create simple vacuum region::

        nodes = [(1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0)]
        region = VacuumRegion(nodes, P=2, material_recycling=0.95)

    Create region with pumping::

        pump_cfg = [{
            'type': 'region',
            'nodes': [(1.8, 0.1), (1.9, 0.1), (1.9, 0.2), (1.8, 0.2)],
            'recycling': 0.0
        }]
        region = VacuumRegion(nodes, P=2, pump=pump_cfg)

    Create region with gas puff::

        puff_cfg = [{
            'type': 'point',
            'location': (1.5, 0.5),
            'current': 1e20,
            'igsp': 0
        }]
        region = VacuumRegion(nodes, P=2, puff=puff_cfg)

    Restore from file::

        region = VacuumRegion('saved_region.pkl')

    Notes
    -----
    - Plasma surfaces must appear first in node list (indices 0 to P-1)
    - Material surfaces follow (indices P to numSurfaces-1)
    - Transport matrix is only computed for plasma surfaces
    - View factors include geometric visibility and distribution shape
    - Multiprocessing dramatically speeds up large geometries (>50 surfaces)
    - Flux conservation should be verified with checkContinuity()
    - Higher reflections increase accuracy but scale as O(reflections)

    See Also
    --------
    Surface : Individual surface element implementation
    VNM_interface : UEDGE integration interface
    """

    def __init__(
        self,
        nodeList,
        P=0,
        r_offset_plasma=1,
        r_offset_material=1,
        multiprocess=True,
        ncores=None,
        verbose=True,
        material_recycling=1,
        pump=None,
        puff=None,
        reflections=1e6,
        hdf5location="vnm",
        savename=None,
        isvacuummodel=None,
        write=False,
        **kwargs,
    ):
        from shapely import Point, Polygon
        from tqdm import tqdm
        from pickle import load
        from multiprocessing import Process, Pipe, Manager
        from os import cpu_count, environ
        from itertools import islice
        from numpy import array, array_split, zeros
        from time import time
        from copy import deepcopy
        from h5py import is_hdf5, File

        self.surfaces = {}
        self.material_recycling = material_recycling
        self.reflections = int(reflections)
        self.nodeList = nodeList
        self.r_offset_plasma = r_offset_plasma
        self.r_offset_material = r_offset_material

        if write:
            if savename is None:
                raise ValueError("'savename' must be specified when write is True")
        if isinstance(nodeList, type(None)):
            raise TypeError("nodeList cannot be None!")

        if isinstance(nodeList, str):
            if is_hdf5(nodeList):  # Restoring from HDF5 file
                with File(nodeList, "r") as f:
                    vnm = f[hdf5location]
                    for var in [
                        "nodeList",
                        "P",
                        "r_offset_material",
                        "r_offset_plasma",
                        "reflections",
                    ]:
                        self.__setattr__(var, vnm[var][()])
                    if "puff" in vnm:
                        setup = {}
                        for puffname in vnm["puffs"].keys():
                            puff[puffname] = {}
                            for var in vnm["puffs"][puffname].keys():
                                puff[puffname][var] = vnm["puffs"][puffname][var][()]
                    if "pump" in vnm:
                        pump = {}
                        for pumpname in vnm["pumps"].keys():
                            pump[pumpname] = {}
                            for var in vnm["pumps"][pumpname].keys():
                                pump[pumpname][var] = vnm["pumps"][pumpname][var][()]

        starttime = time()

        if isinstance(self.nodeList, str):
            # TODO: Check whether requested file exists and is pickle
            with open(self.nodeList, "rb") as f:
                save = load(f)
                self.surfaces = save["surfaces"]
                self.P = save["P"]
                self.nodeList = save["nodeList"]
            P = self.P

            # Create Polygon of Vacuum region for intersect checks
            self.geometry = Polygon(self.nodeList)
        else:
            # Set up surfaces of geometry and the polygon object
            self.P = P
            for i in range(len(self.nodeList)):
                startNode = Point(self.nodeList[i])
                if i == len(self.nodeList) - 1:
                    endNode = Point(self.nodeList[0])
                else:
                    endNode = Point(self.nodeList[i + 1])

                # Have plasma surfaces use a uniform dist., while wall surfaces use a cosine dist.
                if i < self.P:  # Plasma surfaces
                    self.surfaces[i] = Surface(
                        (startNode.x, startNode.y),
                        (endNode.x, endNode.y),
                        i,
                        r_offset=r_offset_plasma,
                    )
                else:  # non-plasma surfaces
                    self.surfaces[i] = Surface(
                        (startNode.x, startNode.y),
                        (endNode.x, endNode.y),
                        i,
                        r_offset=r_offset_material,
                    )

            # Create Polygon of Vacuum region for intersect checks
            self.geometry = Polygon(self.nodeList)

            if multiprocess:
                if ncores is None:
                    ncores = cpu_count()
                else:
                    ncores = min(cpu_count(), ncores)
                environ["UETOOLS_SILENT"] = "1"
                parent_conn, child_conn = Pipe()
                manager = Manager()
                surface_chunks = manager.dict()
                # Create list of surfaces to be calculated by each subprocess
                sublist = array_split(array(range(len(self.surfaces))), ncores)
                # Spawn subprocesses
                subprocesses = []
                print(
                    f"Calculating {len(self.surfaces)} surface couplings on {ncores} threads..."
                )
                for subprocess in sublist:
                    subprocesses.append(
                        Process(
                            target=self.subprocess_execute,
                            args=(
                                surface_chunks,
                                child_conn,
                                list(subprocess),
                                self.surfaces,
                                self.geometry,
                            ),
                            kwargs=({"verbose": verbose}),
                        )
                    )
                    subprocesses[-1].start()
                for subprocess in subprocesses:
                    subprocess.join()
                self.surfaces = deepcopy(surface_chunks)
                environ["UETOOLS_SILENT"] = "0"

            else:
                # Iterate surfaces to identify surface neigbors
                for _, surface in tqdm(self.surfaces.items()):
                    surface.getNeighbors(self.surfaces, self.geometry)

        self.time = time() - starttime
        self.numSurfaces = len(self.surfaces)

        # Dictionary of surface reflection coefficients
        """ Set up recycling coefficients """
        self.R_array = zeros(self.numSurfaces)
        self.R_array[self.P :] = self.material_recycling

        """ Loop through any pumping regions """
        self.pumping_regions = []
        if pump is not None:
            self.pump = pump
            for _pump in pump:
                self.pumping_regions.append(self.create_pump_region(_pump))

        """ Compile Transport Matrix Method matrices """
        self.matrices()
        self.createTeleMatrix()

        """ Loop through any puffs """
        self.puffs = []
        self.puff_vector = zeros((self.P, 6))
        if puff is not None:
            self.puff = puff
            for _puff in puff:
                self.puffs.append(self.create_puff(_puff))
                self.puff_vector[:, _puff["igsp"]] += self.puffs[-1][
                    "puff_vector"
                ].flatten()

        """ Save to pickle if requested """
        if write:
            self.saveVacuumRegion(savename)

    @staticmethod
    def subprocess_execute(output, conn, surflist, surfaces, geometry, verbose=True):
        from os import getpid

        count = 0
        msg = [0.25, 0.50, 0.75, 1]
        # Iterate Process surfaces
        for surfid in surflist:
            if verbose:
                if count / len(surflist) > msg[0]:
                    print(f"Process {getpid()} {msg[0]*100}% completed.")
                    msg.pop(0)
            surfaces[surfid].getNeighbors(surfaces, geometry)
            output[surfid] = surfaces[surfid]
            count += 1
        if verbose:
            print(
                f"Process {getpid()} completed surfaces {surflist[0]}-{surflist[-1]}."
            )

    def matrixPower(self, matrix, power):
        """Raises a given matrix to the specified power."""
        from numpy import zeros, identity
        from scipy.sparse import csr_array, block_array, linalg

        resultMatrix = linalg.matrix_power(matrix, power)

        return resultMatrix

    def matrices(self):
        """Construct Transport Matrix Method (TMM) matrices.

        Builds sparse matrices for iterative particle transport calculation:
        - R: Diagonal recycling coefficient matrix
        - C: View factor (first-flight coupling) matrix
        - A, B: TMM propagation matrices
        - AB: Combined transport operator
        - (AB)^M @ A: Final transport operator with M reflections

        The TMM formulation tracks particle populations on surfaces:
        γ_out[n+1] = (AB)^n @ A @ γ_in[0]

        where γ_in/γ_out are incident/emitted particle flux vectors.

        Side Effects
        ------------
        Creates and stores sparse matrices:
        - self.R_matrix : Diagonal recycling coefficients
        - self.C_matrix : View factors (transposed from C_array)
        - self.A_matrix : [C, 0; 0, I] block matrix
        - self.B_matrix : [R, 0; I-R, I] block matrix
        - self.AB_matrix : A @ B
        - self.AB_power_A : (AB)^reflections @ A

        Notes
        -----
        - All matrices are scipy.sparse.csr_array for memory efficiency
        - Size is 2*numSurfaces (separate incident/emitted populations)
        - Matrix power uses sparse linear algebra for efficiency
        - self.C_array is transposed before creating C_matrix

        See Also
        --------
        createTeleMatrix : Extract plasma-to-plasma transport from full TMM result
        matrixPower : Compute sparse matrix power
        """
        import numpy
        from numpy import zeros, identity, percentile, log, diag
        from scipy.sparse import csr_array, block_array
        import seaborn as sns
        import matplotlib.pyplot as plt

        # Array representations of R and C
        self.C_array = zeros((self.numSurfaces, self.numSurfaces))
        # Populate C array and take transpose
        for surfaceID, surface in self.surfaces.items():  # self.surfaces.items()
            for outputID in surface.neighbors.keys():
                self.C_array[surfaceID][outputID] = surface.neighbors[outputID]["flux"]
        self.C_array = self.C_array.transpose()
        # R and C into sparse matrices
        self.R_matrix = csr_array(diag(self.R_array))
        self.C_matrix = csr_array(self.C_array)
        # Zero and identity sparse matrices
        Zero_matrix = csr_array(zeros((self.numSurfaces, self.numSurfaces)))
        Identity_matrix = csr_array(identity(self.numSurfaces))
        # Create A, B, and AB sparse matrices
        self.A_matrix = block_array(
            [[self.C_matrix, Zero_matrix], [Zero_matrix, Identity_matrix]]
        )
        self.B_matrix = block_array(
            [
                [self.R_matrix, Zero_matrix],
                [Identity_matrix - self.R_matrix, Identity_matrix],
            ]
        )
        # A * B
        self.AB_matrix = self.A_matrix @ self.B_matrix
        # (AB)^M * A
        self.AB_power_A = (
            self.matrixPower(self.AB_matrix, self.reflections) @ self.A_matrix
        )

    def createTeleMatrix(self):
        """Extract plasma-to-plasma transport matrix from full TMM calculation.

        Computes the telematrix T[i,j] giving the fraction of particles emitted
        from plasma surface i that ultimately reach plasma surface j after all
        reflections off material surfaces.

        This is the core output used by UEDGE's vacuum model to couple plasma
        surfaces through neutral transport.

        Side Effects
        ------------
        Creates self.telematrix : ndarray, shape (P, P)
            Transport probability matrix where T[i,j] is the fraction of
            particles leaving plasma surface i that reach plasma surface j.

        Algorithm
        ---------
        For each plasma surface i:
        1. Set γ_in[i] = 1, all other γ_in = 0
        2. Compute γ_out = (AB)^M @ A @ γ_in
        3. Extract T[:, i] = γ_out[P:2P] (emitted from plasma surfaces)

        Notes
        -----
        - Only plasma surfaces (0 to P-1) are included in telematrix
        - Material surfaces (P to numSurfaces-1) are not in output
        - Rows sum to ≤1 (particle loss to material surfaces)
        - Used directly as UEDGE's cftelematrix{pf,w} array

        See Also
        --------
        matrices : Construct TMM matrices including AB_power_A
        """

        from numpy import zeros, transpose

        # Final transport matrix
        self.telematrix = zeros((self.P, self.P))

        gamma_array = zeros((self.numSurfaces * 2, 1))
        for i in range(0, self.P):
            gamma_array[i, 0] = 1
            rowCalculation = self.AB_power_A @ gamma_array
            gamma_array[i, 0] = 0
            gammaOut = rowCalculation[self.numSurfaces :]
            gammaFinal = gammaOut[0 : self.P].flatten()
            for j in range(self.P):
                self.telematrix[j, i] = gammaFinal[j]

    def create_puff(self, puff_setup):
        """Create gas puffing source from configuration dictionary.

        Processes puff configuration, locates nearest material surface, and
        computes puff vector (contribution to each plasma surface after transport).

        Parameters
        ----------
        puff_setup : dict
            Puff configuration with required keys:
            - 'type': 'point' (only supported type)
            - 'location': (R, Z) tuple of puff coordinates [meters]
            - 'current': Injection rate [particles/s]
            - 'igsp': Gas species index (0-based)

        Returns
        -------
        dict
            Puff information with computed transport:
            - 'type': Puff type ('point')
            - 'point': shapely.Point of puff location
            - 'material_surface_index': Index of injection surface
            - 'current': Injection rate [particles/s]
            - 'igsp': Species index
            - 'location': Original (R, Z) tuple
            - 'puff_vector': (P, 1) array of contributions to plasma surfaces

        Raises
        ------
        KeyError
            If required configuration keys missing or invalid puff type.

        Notes
        -----
        - Puff automatically assigned to nearest material surface
        - Uses (AB)^M @ A to propagate injection to plasma surfaces
        - Multiple puffs are additive (summed into puff_vector)
        - Only material surfaces (P ≤ index < numSurfaces) can receive puffs
        - 'point' is currently the only supported puff type

        See Also
        --------
        create_pump_region : Create pumping region
        """

        from shapely import Point
        from numpy import argsort, zeros

        if "type" not in puff_setup:
            raise KeyError(f"Define a puff type")
        ret = {"type": puff_setup["type"]}

        if puff_setup["type"] == "point":
            for key in ["location", "current", "igsp"]:
                if key not in puff_setup:
                    raise KeyError(
                        "Required 'point' puff setup entry" + f" '{key}' not found."
                    )
            try:
                ret["point"] = Point(puff_setup["location"])
            except Exception:
                raise
            nodes = list(self.geometry.exterior.coords)[self.P : -1]
            dists = [self.P + ret["point"].distance(Point(n)) for n in nodes]
            ret["material_surface_index"] = argsort(dists)[:2].max() + self.P
            ret["current"] = puff_setup["current"]
            ret["igsp"] = puff_setup["igsp"]
            ret["location"] = puff_setup["location"]
            drive = zeros((self.numSurfaces * 2, 1))
            drive[ret["material_surface_index"], 0] = ret["current"]
            ret["puff_vector"] = (self.AB_power_A @ drive)[self.numSurfaces :][: self.P]

        else:
            raise KeyError(
                f"Puff type '{pump_setup['type']}' not recognized!"
                + "\nAvailable options are: 'point'"
            )
        return ret

    def create_pump_region(self, pump_setup):
        """Create pumping region from configuration dictionary.

        Processes pump configuration, identifies intersecting surfaces, and
        modifies their recycling coefficients. Later pumps override earlier
        ones for overlapping surfaces.

        Parameters
        ----------
        pump_setup : dict
            Pump configuration with required keys:
            - 'type': 'region' (only supported type)
            - 'nodes': Array-like of (R, Z) node tuples defining pump polygon
            - 'recycling': Recycling coefficient for surfaces in pump region

        Returns
        -------
        dict
            Pump information:
            - 'type': Pump type ('region')
            - 'polygon': shapely.Polygon of pump region
            - 'nodes': Node array as provided
            - 'pumped_segments': List of surface indices with modified recycling
            - 'recycling': Applied recycling coefficient

        Raises
        ------
        KeyError
            If required configuration keys missing or invalid pump type.
        AttributeError
            If fewer than 3 nodes provided (cannot form polygon).

        Notes
        -----
        - Only material surfaces (P ≤ index < numSurfaces) can be pumped
        - Recycling modified in self.R_array for intersecting surfaces
        - Multiple pumps: later in list take precedence for overlapping regions
        - Recycling coefficient ∈ [0, 1]:
          * 0: Perfect pump (no particle return)
          * 1: No pumping (full recycling)
        - 'region' is currently the only supported pump type
        - Geometric intersection tested between surface segments and pump polygon

        Examples
        --------
        Create a pumping region::

            pump_cfg = {
                'type': 'region',
                'nodes': [(1.8, 0.1), (1.9, 0.1), (1.9, 0.2), (1.8, 0.2)],
                'recycling': 0.0
            }
            pump_info = region.create_pump_region(pump_cfg)
            print(f"Pumping {len(pump_info['pumped_segments'])} surfaces")
            # Output: Pumping 2 surfaces

        See Also
        --------
        create_puff : Create gas puffing source
        """

        from shapely import Polygon, intersects, Point

        if "type" not in pump_setup:
            raise KeyError(f"Define a pump type")
        ret = {"type": pump_setup["type"]}

        if pump_setup["type"] == "region":
            for key in ["nodes", "recycling"]:
                if key not in pump_setup:
                    raise KeyError(
                        "Required 'region' pump setup entry" + f" '{key}' not found."
                    )
            if len(pump_setup["nodes"]) < 3:
                raise AttributeError(
                    "Too few nodes provided: provide "
                    + "at least three nodes for Polygon!"
                )
            ret["polygon"] = Polygon(pump_setup["nodes"])
            ret["nodes"] = pump_setup["nodes"]
            ret["pumped_segments"] = []
            ret["recycling"] = pump_setup["recycling"]
            for i in range(self.P, self.numSurfaces):
                if intersects(self.surfaces[i].segment, ret["polygon"]):
                    self.R_array[i] = ret["recycling"]
                    ret["pumped_segments"].append(i)
        else:
            raise KeyError(
                f"Pump type '{pump_setup['type']}' not recognized!"
                + "\nAvailable options are: 'region'"
            )
        return ret

    def saveVacuumRegion(self, savename):
        from pickle import dump

        """Use to save a Vacuum Region to avoid having to generate a new one every time. 
            Be sure to set pf/main (tokamakPlot), r_offset (Surface constructor), and 
            variation (Vacuum Region constructor)."""

        save = {"surfaces": self.surfaces, "P": self.P, "nodeList": self.nodeList}
        with open(savename, "wb") as f:
            dump(save, f)

    def writeVacuumSetup(self, obj, write_matrices=True):
        from h5py import File, Group

        def save_csr(group, mat):
            for var in ["data", "indices", "indptr"]:
                if var in group:
                    del group[var]
                group.create_dataset(var, data=mat.__getattribute__(var))
            group.attrs["shape"] = mat.shape

        variables = [
            "nodeList",
            "P",
            "r_offset_plasma",
            "r_offset_material",
            "material_recycling",
            "reflections",
            "telematrix",
        ]
        csr = [
            "R_matrix",
            "C_matrix",
            "A_matrix",
            "B_matrix",
            "AB_matrix",
            "AB_power_A",
        ]

        if isinstance(obj, str):
            f = File(obj, "a")
            vnm = f.require_group("vnm")
        elif isinstance(obj, Group):
            vnm = obj
        else:
            raise TypeError(
                "Object must be Group or file name as string where to write."
            )

        # Store regular data
        for var in variables:
            if var in vnm:
                del vnm[var]
            vnm.create_dataset(var, data=self.__getattribute__(var))
        if write_matrices:
            # Store CSR matrices
            for mat in csr:
                matgroup = vnm.require_group(mat)
                save_csr(matgroup, self.__getattribute__(mat))
        # Store pump setups
        if len(self.pumping_regions) > 0:
            pumps = vnm.require_group("pumps")
        pumpid = 1
        for setup in self.pumping_regions:
            pumpname = f"pump_{pumpid}"
            pump = pumps.require_group(pumpname)
            for var, data in setup.items():
                if var not in ["polygon"]:
                    if var in pump:
                        del pump[var]
                    pump.create_dataset(var, data=data)
            pumpid += 1
        # Store puff setups
        if len(self.puffs) > 0:
            puffs = vnm.require_group("puffs")
        puffid = 1
        for setup in self.puffs:
            puffname = f"puff_{puffid}"
            puff = puffs.require_group(puffname)
            for var, data in setup.items():
                if var not in ["point"]:
                    if var in puff:
                        del puff[var]
                    puff.create_dataset(var, data=data)

        if isinstance(obj, str):
            f.close()

    def checkContinuity(self, verbose=True):
        """For identifying errors in flux unity for surfaces in the geometry."""
        self.errors = []
        for surfid, surface in self.surfaces.items():
            if abs(surface.totflux - 1) > 1e-6:
                if verbose:
                    surface.printReport()
                self.errors.append(surfid)
        return len(self.errors) == 0

    def plotGeometry(
        self,
        ax=None,
        labels=False,
        testsurf=[],
        showCircle=False,
        markers=None,
        connectionLineWidth=0.5,
        **kwargs,
    ):
        from matplotlib.pyplot import subplots, Figure, Axes, ioff
        import matplotlib.pyplot as plt

        """Plot the full geometrical representation of the VacuumRegion."""

        if isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        elif ax is None:
            f, ax = subplots(figsize=(5, 10))
        elif not isinstance(ax, Axes):
            raise Exception("ax not a valid Figure or Axes object")

        if isinstance(testsurf, int):
            testsurf = [testsurf]

        ioff()

        for surfid, surface in self.surfaces.items():
            color = "k"
            if surface.ID < self.P:
                color = "red"
            else:
                if self.R_array[surfid] != self.material_recycling:
                    color = "grey"
            surface.plotSelf(color=color, ax=ax, label=labels, showCircle=showCircle)

        for itest in testsurf:
            self.surfaces[itest].plotConnections(
                ax=ax, linewidth=connectionLineWidth, **kwargs
            )
        c = 0
        N_regions = len(self.pumping_regions) + len(self.puffs)
        colors = [plt.get_cmap("rainbow")(i / N_regions) for i in range(N_regions)]
        for region in self.pumping_regions:
            for i in region["pumped_segments"]:
                self.surfaces[i].plotSelf(color=colors[c], ax=ax, showCircle=False)
            c += 1
        for puff in self.puffs:
            ax.plot(
                puff["point"].xy[0][0], puff["point"].xy[1][0], "o", color=colors[c]
            )
            self.surfaces[puff["material_surface_index"]].plotSelf(
                color=colors[c], ax=ax, showCircle=False
            )
            c += 1

        for line in ax.lines:
            line.set_marker("")

        ax.set_aspect("equal")
        ax.grid(False)
        plt.xlabel("R [m]")
        plt.ylabel("Z [m]")

        # plt.savefig('fullGeometry62.svg', dpi=300)

        plt.show(block=False)

        return ax.get_figure()

    def heatmapPlot(self):
        import numpy
        from numpy import zeros, identity, percentile, log
        from scipy.sparse import csr_array, block_array
        import seaborn as sns
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm

        """Creates a heatmap of C, R, and Transport matrices."""

        """Print statements to check for unity of transport matrix."""

        # Plotting heatmaps of C, R, and Output (Transport)
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        fig.suptitle("Cosine Distribution", fontsize=16)

        matricesToPlot = [self.C_array, self.R_array, self.telematrix]
        matricesToPlotNames = ["C", "R", "Output"]
        for i, ax in enumerate(axes):
            sns.heatmap(
                matricesToPlot[i],
                cmap="jet",
                annot=False,
                ax=ax,
                norm=LogNorm(vmin=1e-5, vmax=1),
            )
            ax.set_title(matricesToPlotNames[i])
            ax.set_aspect("equal")
            if matricesToPlotNames[i] == "Output":
                ax.set_xlabel("Source Surfaces")
                ax.set_ylabel("Receiving Surfaces")
            else:
                ax.set_xlabel("Receiving Surfaces")
                ax.set_ylabel("Source Surfaces")

        plt.tight_layout()
        plt.show(block=False)

        return self.telematrix


class Surface:
    """Individual surface element for vacuum region modeling.

    Represents a line segment surface in 2D (R,Z) space with associated geometric
    properties for particle emission and view factor calculations. Implements
    angular distribution modeling (uniform or cosine) via a distribution circle
    positioned relative to the surface.

    The distribution circle determines the angular probability distribution for
    particles leaving the surface:
    - r_offset=1 (cosine): Circle tangent to surface → cosine distribution (physical)
    - r_offset=0 (uniform): Circle center on surface → uniform distribution
    - 0<r_offset<1: Intermediate distribution

    View factors to neighboring surfaces are computed via geometric overlap between
    the distribution circle and flux triangles formed by line-of-sight connections.

    Parameters
    ----------
    start : tuple
        (R, Z) coordinates of surface start point [meters].
    end : tuple
        (R, Z) coordinates of surface end point [meters].
    ID : int
        Unique surface identifier within vacuum region.
    material : int, optional
        Material type (currently unused, default: 1).
    emitting : int, optional
        Emission flag (currently unused, default: 0).
    absorbing : int, optional
        Absorption flag (currently unused, default: 0).
    r_offset : float, optional
        Distribution circle offset (default: 1):
        - 0: Uniform distribution
        - 1: Cosine (Knudsen) distribution
        - (0,1): Intermediate distribution

    Attributes
    ----------
    start, end : shapely.Point
        Surface endpoint coordinates.
    segment : shapely.LineString
        Line segment representation of surface.
    ID : int
        Surface identifier.
    surfaceLength : float
        Euclidean length of surface [meters].
    midpoint : shapely.Point
        Surface center point.
    normalStart, normalEnd : shapely.Point
        Start and end points of outward normal vector.
    normal : shapely.LineString
        Outward normal vector from surface midpoint.
    dx, dy : float
        Surface direction components (end - start).
    circle : shapely.Polygon
        Distribution circle for view factor calculations.
    dCircleCenter : shapely.Point
        Center of distribution circle.
    r_offset : float
        Applied distribution offset.
    distType : str
        Distribution description ("Uniform Distribution", "Cosine Distribution", etc.).
    neighbors : dict
        View factors to neighboring surfaces:
        {neighbor_ID: {'flux': float, 'los': shapely.Polygon}}
        where 'flux' is the view factor and 'los' is the line-of-sight triangle.
    totflux : float
        Sum of all view factors (should be ≈1 for flux conservation).
    epsilon : float
        Numerical tolerance for geometric operations (default: 1e-5).

    Methods
    -------
    distributionCircle(r_offset)
        Create distribution circle for given offset.
    intersectionArea(s2)
        Compute view factor to another surface s2.
    getNeighbors(surfaces, geometry)
        Find all visible neighbors and compute view factors.
    getSmallestIntersectAngle(neighbor, geometry, polygon)
        Adjust flux triangle for line-of-sight obstructions.
    drawOuterCircle()
        Create outer comparison circle for analytical validation.
    printReport()
        Print view factors to all neighbors.
    plotSelf(ax=None, color='k', label=False, showCircle=True)
        Plot surface and distribution circle.
    plotConnections(ax=None, linewidth=2, **kwargs)
        Plot view factor triangles to all neighbors.
    showTwoSurfacePlot(s2, r_offset=0)
        Interactive plot showing two-surface geometry and flux triangle.
    showOuterCirclePlot(r_offset=1)
        Plot surface with outer circle for analytical comparison.
    showAnalyticPlot(ax, comparison=True, r_offset=1, showBothDist=False)
        Plot view factor vs angle and compare to analytical distribution.
    analyticUniform(ax)
        Plot analytical uniform distribution for comparison.
    normalHelper(dx, dy, endX, endY, init)
        Compute outward normal direction from surface orientation.
    vectorHelper(start, end)
        Create numpy vector from two points.
    dotProductAngle(v1, v2)
        Compute signed angle between two vectors using dot product.

    Examples
    --------
    Create surface and compute properties::

        s1 = Surface((1.0, 0.0), (2.0, 0.0), ID=0, r_offset=1)
        print(f"Length: {s1.surfaceLength:.3f} m")
        # Output: Length: 1.000 m
        print(f"Midpoint: ({s1.midpoint.x:.3f}, {s1.midpoint.y:.3f})")
        # Output: Midpoint: (1.500, 0.000)

    Compute view factor to another surface::

        s2 = Surface((1.5, 0.5), (2.5, 0.5), ID=1, r_offset=1)
        flux, triangle = s1.intersectionArea(s2)
        print(f"View factor: {flux:.4f}")
        # Output: View factor: 0.1234

    Notes
    -----
    - Distribution circle radius is surfaceLength/84 (empirically chosen)
    - View factors include both geometric visibility and angular distribution
    - Line-of-sight obstructions are handled by geometric intersection tests
    - For physical thermal emission, use r_offset=1 (cosine distribution)
    - Total view factors should sum to 1.0 for flux conservation
    - The epsilon parameter prevents numerical issues in geometric calculations

    See Also
    --------
    VacuumRegion : Container for multiple surfaces
    """

    def __init__(self, start, end, ID, material=1, emitting=0, absorbing=0, r_offset=1):
        from shapely import Point, LineString, plotting
        from matplotlib.pyplot import subplots
        import math

        # Start and end points of the surface and a segment representation of the surface
        self.start = Point(start[0], start[1])
        self.end = Point(end[0], end[1])

        # LineString representing the surface
        self.segment = LineString([start, end])
        self.ID = ID

        self.surfaceLength = math.sqrt(
            (self.end.x - self.start.x) ** 2 + (self.end.y - self.start.y) ** 2
        )

        # Midpoint of surface
        self.midpoint = self.segment.centroid

        # End point of the normal vector
        self.normalEndX = self.midpoint.x
        self.normalEndY = self.midpoint.y

        # Getting the correct slope for the normal vector
        self.dx = self.end.x - self.start.x
        self.dy = self.end.y - self.start.y

        # Creating the normal vector
        self.normalHelper(
            abs(self.dx), abs(self.dy), self.normalEndX, self.normalEndY, True
        )

        # Start and end points, and the normal line itself (use these for reference)
        self.normalStart = Point(self.midpoint.x, self.midpoint.y)
        self.normalEnd = Point(self.normalEndX, self.normalEndY)
        self.normal = LineString([self.normalStart, self.normalEnd])

        # Additional aspects of the surface
        self.circle = None

        # Coupling to surfaces with LOS
        self.neighbors = {}
        self.totflux = 0

        # Creating the distribution circle given an offset (r_offset = 1 for cosine, 0 for uniform)
        self.distributionCircle(r_offset)

        self.epsilon = 1e-5  # use as a reference for buffering works when epsilon = 0.00001 (1e-5) --> Use for adjusting line of sight

        return

    def distributionCircle(self, r_offset):
        """Create angular distribution circle for view factor calculations.

        Constructs a circle whose geometric overlap with flux triangles determines
        the angular probability distribution for particles emitted from this surface.

        Parameters
        ----------
        r_offset : float
            Distribution circle offset relative to radius:
            - 0: Circle center on surface → uniform angular distribution
            - 1: Circle tangent to surface → cosine (Knudsen) distribution
            - (0, 1): Intermediate distribution

        Side Effects
        ------------
        Sets attributes:
        - self.circle : shapely.Polygon of distribution circle (500 points)
        - self.dCircleCenter : shapely.Point at circle center
        - self.r_offset : Stored offset value
        - self.distType : Human-readable distribution name string

        Notes
        -----
        - Circle radius is surfaceLength / 84 (empirically chosen)
        - Circle center is offset along outward normal by r_offset * radius
        - 500 points used to approximate circle (high precision)
        - For cosine distribution (r_offset=1), midpoint added explicitly
        - Circle winding order matches surface direction (clockwise/counterclockwise)
        - Physical thermal emission uses r_offset=1 (cosine law)

        Algorithm
        ---------
        1. Compute circle center: midpoint + r_offset * radius * normal_direction
        2. Generate 500 points around circle at angles 0 to 2π
        3. Special case: If cosine (r_offset=1), insert surface midpoint
        4. Create shapely.Polygon from circle points

        See Also
        --------
        intersectionArea : Use distribution circle to compute view factor
        """

        from shapely import Point, LineString, plotting, LinearRing, Polygon
        from matplotlib.pyplot import subplots
        import math

        # Finding radius and center of circle
        self.r_offset = r_offset
        radius = self.surfaceLength / 84
        self.dCircleCenter = self.normal.interpolate(self.r_offset * radius)

        # Visual plot labeling
        if r_offset == 0:
            self.distType = "Uniform Distribution"
        elif r_offset == 1:
            self.distType = "Cosine Distribution"
        else:
            self.distType = "Distribution"

        # Getting the Points of the circle
        numPoints = 500
        circlePoints = []

        if (
            self.end.x < self.start.x
        ):  # builds circle clockwise rather than counter clockwise
            for i in range(numPoints, -1, -1):
                # Calculates every angle from 0-2pi, placing 500 points to create the circle
                angle = (2 * math.pi) * (i / numPoints)
                x = self.dCircleCenter.x + radius * math.cos(angle)
                y = self.dCircleCenter.y + radius * math.sin(angle)
                if angle == math.pi and self.distType == "Cosine Distribution":
                    circlePoints.append((self.midpoint.x, self.midpoint.y))
                    circlePoints = sorted(
                        circlePoints, key=lambda x: x[0], reverse=True
                    )
                    i += 1
                circlePoints.append((x, y))
        else:
            for i in range(numPoints + 1):
                angle = (2 * math.pi) * (i / numPoints)
                x = self.dCircleCenter.x + radius * math.cos(angle)
                y = self.dCircleCenter.y + radius * math.sin(angle)
                if angle == math.pi and self.distType == "Cosine Distribution":
                    circlePoints.append((self.midpoint.x, self.midpoint.y))
                    circlePoints = sorted(
                        circlePoints, key=lambda x: x[0], reverse=True
                    )
                    i += 1
                circlePoints.append((x, y))

        # Object representation of the distribution circle
        self.circle = Polygon(circlePoints)

        return

    def printReport(self):
        print("Surface {}: {}".format(self.ID, self.totflux))
        for neighid, neighbor in self.neighbors.items():
            print("    -> {}: {}".format(f"{neighid}".rjust(4), neighbor["flux"]))

    def intersectionArea(self, s2):
        """Compute view factor from this surface to another surface.

        Calculates the fraction of particles emitted from this surface (self)
        that directly reach another surface (s2) without intermediate reflections.
        View factor is determined by geometric overlap between the distribution
        circle and the flux triangle connecting the two surfaces.

        Parameters
        ----------
        s2 : Surface
            Target surface to compute view factor toward.

        Returns
        -------
        fractionalArea : float
            View factor (0 to 1): fraction of emitted particles reaching s2.
        triangle : shapely.Polygon
            Flux triangle connecting surfaces (for visualization/debugging).

        Side Effects
        ------------
        Sets temporary attributes (overwritten in subsequent calls):
        - self.triangle : Flux triangle from self.midpoint to s2 endpoints
        - self.leg1, self.leg2 : LineStrings forming triangle sides
        - self.overlapShape : Intersection of triangle and distribution circle
        - self.vLeg1, self.vLeg2 : Vector representations of triangle legs
        - self.totalAreaCircle : Relevant hemisphere of distribution circle

        Raises
        ------
        None
            Prints warning if distribution circle not created.

        Algorithm
        ---------
        1. Form flux triangle: self.midpoint → s2.start → s2.end → self.midpoint
        2. Compute intersection: overlapShape = triangle ∩ distribution_circle
        3. For offset distributions (0 < r_offset < 1):
           - Split circle at self.segment
           - Use only hemisphere facing s2
        4. fractionalArea = overlapShape.area / relevant_circle_area

        Notes
        -----
        - Requires self.distributionCircle() called first
        - Does not account for line-of-sight obstructions (see getNeighbors)
        - View factor depends on:
          * Distance between surfaces
          * Relative orientation (normal directions)
          * Angular distribution (via r_offset)
        - For conservation: Σ fractionalArea over all s2 should equal 1.0
        - Epsilon buffering handles numerical precision in circle splitting

        Examples
        --------
        Compute view factor between two surfaces::

            s1 = Surface((1, 0), (2, 0), ID=0, r_offset=1)
            s2 = Surface((1.5, 0.5), (2.5, 0.5), ID=1)
            flux, triangle = s1.intersectionArea(s2)
            print(f"View factor from s1 to s2: {flux:.4f}")
            # Output: View factor from s1 to s2: 0.1234

        See Also
        --------
        distributionCircle : Create distribution circle for view factors
        getNeighbors : Account for line-of-sight obstructions
        """

        from shapely import Point, LineString, plotting, Polygon, is_closed
        import math

        if self.circle == None:
            print(
                "Call distributionCircle on Surface before finding intersection area!"
            )
            return

        # Triangle of flux
        self.triangle = Polygon(
            [s2.start, s2.end, self.midpoint, s2.start]
        )  # from midpoint of self to the endpoints of s2

        # Legs of the triangle
        self.leg1 = LineString([self.midpoint, s2.start])
        self.leg2 = LineString([self.midpoint, s2.end])

        # Overlap of triangle and distribution circle
        self.overlapShape = self.triangle.intersection(self.circle)

        # Vector representations of triangle legs
        self.vLeg1 = self.vectorHelper(
            (self.midpoint.x, self.midpoint.y), (s2.start.x, s2.start.y)
        )
        self.vLeg2 = self.vectorHelper(
            (self.midpoint.x, self.midpoint.y), (s2.end.x, s2.end.y)
        )

        # Getting the correct area of the distribution circle on one side of the normal line
        overlapArea = self.overlapShape.area
        if (0 <= self.r_offset and self.r_offset < 1) and self.circle.intersects(
            self.segment
        ):
            # Circle shifted between uniform and cosine
            buffLine = self.segment.buffer(self.epsilon * 1e-7)
            # Split the distribution circle and the surface line
            splitdCircle = self.circle.difference(buffLine)
            if splitdCircle.geoms[0].intersects(self.normal):
                self.totalAreaCircle = splitdCircle.geoms[0]
            else:
                self.totalAreaCircle = splitdCircle.geoms[1]
            circleArea = self.totalAreaCircle.area
        else:
            circleArea = self.circle.area

        # Calculate the flux (fractional area)
        fractionalArea = overlapArea / circleArea

        return fractionalArea, self.triangle

    def getNeighbors(self, surfaces, geometry):
        """Find visible neighbor surfaces and compute view factors with line-of-sight.

        Iterates through all surfaces to identify which are visible from this surface,
        accounting for geometric obstructions. Computes view factors and stores them
        in self.neighbors dictionary. This is the core geometry calculation for TMM.

        Parameters
        ----------
        surfaces : dict
            Dictionary of all Surface objects in vacuum region, keyed by surface ID.
        geometry : shapely.Polygon
            Closed polygon defining vacuum region boundary for obstruction tests.

        Side Effects
        ------------
        Updates attributes:
        - self.neighbors : dict
            {neighbor_ID: {'flux': float, 'los': shapely.Polygon}}
            View factors and line-of-sight triangles for visible neighbors.
        - self.totflux : float
            Sum of all view factors (should be ≈1.0 for conservation).

        Algorithm
        ---------
        For each potential neighbor surface:
        1. Compute initial view factor via intersectionArea()
        2. Check flux triangle intersects neighbor's outward normal (front side)
        3. Check flux triangle intersects this surface's outward normal
        4. Find geometric obstructions (triangle intersections with exterior)
        5. Adjust triangle legs to avoid obstructions (getSmallestIntersectAngle)
        6. Recompute view factor with adjusted geometry
        7. Store in self.neighbors if flux > 0

        Notes
        -----
        - Self-coupling (i == j) is skipped
        - Back-side coupling is rejected (triangle doesn't intersect normal)
        - Line-of-sight obstructions reduce view factors
        - Flux conservation requires Σ flux ≈ 1.0 (check with totflux)
        - Handles complex geometries with MultiPolygon obstructions
        - Uses epsilon buffering for numerical robustness

        See Also
        --------
        intersectionArea : Compute raw geometric view factor
        getSmallestIntersectAngle : Adjust for line-of-sight obstructions
        checkContinuity : Verify flux conservation
        """

        from shapely import (
            intersects,
            difference,
            crosses,
            buffer,
            contains,
            intersection,
        )

        #  Fractional Area Calculations
        for neighid, neighbor in surfaces.items():

            if self.ID != neighid:
                flux, triangle = self.intersectionArea(neighbor)

                # Flux would go to the wrong (back) side of the surface
                try:
                    if (
                        not intersects(triangle, neighbor.normal)
                        or triangle.intersection(neighbor.normal) == neighbor.midpoint
                    ):
                        continue
                except:
                    continue

                if not intersects(triangle, self.normal):
                    continue

                # Polygon/multipolygon of the intersections of the outside area of the geometry w/ the triangle legs
                AllIntersections = difference(triangle, geometry)
                self.AllIntersections = difference(triangle, geometry)

                self.leg1SmallestPoint = (neighbor.start.x, neighbor.start.y)
                self.leg2SmallestPoint = (neighbor.end.x, neighbor.end.y)

                self.vNewLeg1 = self.vLeg1
                self.vNewLeg2 = self.vLeg2

                # Angle between the legs used for line of sight calculations
                self.smallestAngle = self.dotProductAngle(self.vNewLeg1, self.vNewLeg2)

                #  Only one area of intersection with the triangle
                if AllIntersections.geom_type == "Polygon":
                    if not self.getSmallestIntersectAngle(
                        neighbor, geometry, AllIntersections
                    ):
                        continue

                # Multiple intersections with the outside region and the triangle/triangle legs
                elif (
                    AllIntersections.geom_type == "MultiPolygon"
                    or AllIntersections.geom_type == "GeometryCollection"
                ):
                    # Check all possible polygons
                    for polygon in AllIntersections.geoms:
                        if not self.getSmallestIntersectAngle(
                            neighbor, geometry, polygon
                        ):
                            continue

                # Set new S2 Surface after line of sight adjustments and calculate the new flux. Update the total flux out of Surface 1.
                newS2 = Surface(
                    (self.leg1SmallestPoint[0], self.leg1SmallestPoint[1]),
                    (self.leg2SmallestPoint[0], self.leg2SmallestPoint[1]),
                    neighbor.ID,
                )
                flux, triangle = self.intersectionArea(newS2)
                self.totflux += flux
                if flux > 0:
                    if neighid not in self.neighbors:
                        self.neighbors[neighid] = {}
                    self.neighbors[neighid]["flux"] = flux
                    self.neighbors[neighid][
                        "los"
                    ] = triangle  # Get the new triangle shape and add it here

    def getSmallestIntersectAngle(self, neighbor, geometry, polygon):
        from shapely import intersects, crosses, contains, buffer, intersection, Point

        """Adjust the legs of the flux triangle in accordance with flux calculations."""

        if contains(
            buffer(geometry, self.epsilon), polygon
        ):  # When polygon is on the border and not at all outside the geometry
            return True

        if crosses(polygon, buffer(self.leg1, self.epsilon * 1e-1)) and crosses(
            polygon, buffer(self.leg2, self.epsilon * 1e-1)
        ):
            return False

        # If Leg 1 (from S1 midpoint to S2 start) intersects the outside region
        if intersects(polygon, buffer(self.leg1, self.epsilon * 1e-1)) and (
            contains(geometry, self.leg1) == False
        ):

            if polygon.geom_type == "LineString":
                coordinates1 = list(polygon.coords)
            else:
                coordinates1 = list(polygon.exterior.coords)

            leg1Start = buffer(self.midpoint, self.epsilon)
            leg1End = buffer(neighbor.start, self.epsilon)

            if (
                contains(
                    leg1Start,
                    intersection(polygon, buffer(self.leg1, self.epsilon * 1e-2)),
                )
                == False
            ) and (
                contains(
                    leg1End,
                    intersection(polygon, buffer(self.leg1, self.epsilon * 1e-2)),
                )
                == False
            ):
                for pair in coordinates1:
                    if Point(pair) == self.midpoint or contains(
                        buffer(self.midpoint, self.epsilon * 1e2), Point(pair)
                    ):
                        continue
                    self.vNewLeg1 = self.vectorHelper(
                        (self.midpoint.x, self.midpoint.y), pair
                    )
                    self.vNewLeg2 = self.vectorHelper(
                        (self.midpoint.x, self.midpoint.y),
                        (self.leg2SmallestPoint[0], self.leg2SmallestPoint[1]),
                    )
                    newAngle = self.dotProductAngle(self.vNewLeg1, self.vNewLeg2)

                    # If the legs cross over during adjustment
                    if (abs(newAngle * self.smallestAngle) > 0) and (
                        newAngle * self.smallestAngle < 0
                    ):
                        self.leg1SmallestPoint = self.leg2SmallestPoint
                        self.smallestAngle = 0
                        return False

                    if abs(newAngle) < abs(self.smallestAngle):
                        self.smallestAngle = newAngle
                        self.leg1SmallestPoint = pair

        # If Leg 2 (S1 midpoint to S2 end) intersects the outside region
        if intersects(polygon, buffer(self.leg2, self.epsilon * 1e-1)) and (
            contains(geometry, self.leg2) == False
        ):
            if polygon.geom_type == "LineString":
                coordinates2 = list(polygon.coords)
            else:
                coordinates2 = list(polygon.exterior.coords)

            leg2Start = buffer(self.midpoint, self.epsilon)
            leg2End = buffer(neighbor.end, self.epsilon)

            if (
                contains(
                    leg2Start,
                    intersection(polygon, buffer(self.leg2, self.epsilon * 1e-2)),
                )
                == False
            ) and (
                contains(
                    leg2End,
                    intersection(polygon, buffer(self.leg2, self.epsilon * 1e-2)),
                )
                == False
            ):
                for pair in coordinates2:
                    if Point(pair) == self.midpoint or contains(
                        buffer(self.midpoint, self.epsilon * 1e2), Point(pair)
                    ):
                        continue
                    self.vNewLeg1 = self.vectorHelper(
                        (self.midpoint.x, self.midpoint.y),
                        (self.leg1SmallestPoint[0], self.leg1SmallestPoint[1]),
                    )
                    self.vNewLeg2 = self.vectorHelper(
                        (self.midpoint.x, self.midpoint.y), pair
                    )
                    newAngle = self.dotProductAngle(self.vNewLeg1, self.vNewLeg2)

                    # If legs cross over each other during adjustment
                    if (abs(newAngle * self.smallestAngle) > 0) and (
                        newAngle * self.smallestAngle < 0
                    ):
                        self.leg2SmallestPoint = self.leg1SmallestPoint
                        self.smallestAngle = 0
                        return False

                    if abs(newAngle) < abs(self.smallestAngle):
                        self.smallestAngle = newAngle
                        self.leg2SmallestPoint = pair

        return True

    def drawOuterCircle(self):
        from shapely import (
            Point,
            plotting,
            Polygon,
            MultiPoint,
            is_closed,
            get_coordinates,
            LineString,
        )
        import math

        """Creates the distribution circle that will be compared to the analytic cosine (does not show plot)."""
        """Reference Variables/Important:
            self.outerCircle (Polygon)
        """

        # Creating the outer circle of S2 surfaces
        outerRadius = (
            self.surfaceLength / 2
        )  # outer circle has diameter equal to the length of self surface (S1)
        outerCenter = self.midpoint
        outerCirclePoints = []
        numPoints = 500

        for i in range(numPoints):
            angle = (2 * math.pi) * (i / numPoints)
            x = outerCenter.x + outerRadius * math.cos(
                angle
            )  # uses same center as the distribution circle
            y = outerCenter.y + outerRadius * math.sin(angle)
            outerCirclePoints.append((x, y))

        # Polygon object of the full outer circle, before splitting to correct side of normal
        self.fullCircle = Polygon(outerCirclePoints)

        buffLine = self.segment.buffer(self.epsilon * 1e-7)
        splitCircles = self.fullCircle.difference(buffLine)
        if splitCircles.geoms[0].intersects(self.normal):
            self.outerCircle = splitCircles.geoms[0]  # desired circle
        else:
            self.outerCircle = splitCircles.geoms[1]  # desired circle

        return

    # # # # # # # # # # #
    # PLOTTING FUNCTIONS #
    # # # # # # # # # # #

    def plotSelf(self, ax=None, color="k", label=False, showCircle=True):
        from matplotlib.pyplot import subplots, ioff, Figure, Axes
        from shapely import plotting, buffer, Point

        """Plots a surface."""

        ioff()
        if ax is None:
            fig, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        elif not isinstance(ax, Axes):
            raise Exception("ax not a valid Figure or Axes object")
        ax.set_aspect("equal")

        plotting.plot_line(self.segment, ax, color=color, linewidth=2)
        if label:
            ax.text(
                self.midpoint.x,
                self.midpoint.y,
                f"Surface {self.ID}",
                fontsize=8,
                color="black",
            )

        if showCircle:
            plotting.plot_polygon(
                self.circle, ax, add_points=False, color="green", linewidth=1
            )

        return ax

    def plotConnections(self, ax=None, linewidth=2, colorseed=1, **kwargs):
        from matplotlib.pyplot import subplots, ioff, Figure, Axes, get_cmap
        from shapely import plotting
        import random
        from numpy import linspace

        ioff()
        if ax is None:
            fig, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        elif not isinstance(ax, Axes):
            raise Exception("ax not a valid Figure or Axes object")
        ax.set_aspect("equal")

        cmap = get_cmap("jet")
        cols = linspace(0, 1, len(self.neighbors))
        random.Random(colorseed).shuffle(cols)
        colors = iter(cmap(cols))
        for neighid, neighbor in self.neighbors.items():
            plotting.plot_polygon(
                neighbor["los"],
                add_points=False,
                color=next(colors),
                linewidth=linewidth,
                **kwargs,
            )

    def showAnalyticPlot(self, ax, comparison=True, r_offset=1, showBothDist=False):
        from shapely import (
            Point,
            plotting,
            Polygon,
            MultiPoint,
            is_closed,
            get_coordinates,
            LineString,
        )
        from shapely.plotting import plot_points
        from matplotlib.pyplot import subplots, ioff
        import matplotlib.pyplot as plt
        import math
        import numpy as np

        """Plots the curve from the outer circle and compares it to the analytic equation plot."""

        """The 'comparison' variable determines if we are comparing to the analytic cosine distribution.
            Only set to true if modeling a COSINE distribution (offset = 1) 
            and you want to compare it to the analytic."""

        """Set showBothDist to True to display both uniform and cosine (geometric and analytic) on the same plot.
            to do this, set r_offset to 1. Also be sure to call analyticUniform on self."""

        ioff()
        fig = subplots()
        ax.set_aspect("auto")
        plt.xlabel("Angle (Radians)")
        plt.ylabel("Fractional Area")

        self.distributionCircle(r_offset)
        self.drawOuterCircle()

        # CREATING THE PLOT OF ANGLE VS AREA
        s2Points = get_coordinates(
            self.outerCircle
        )  # points defining the outer circle of S2 surfaces
        s2Start = s2Points[0]  # The starting point of the first S2 surface
        plotPoints = []  # Points to plot for the outer circle plot
        pdfArea = 0  # "C" for the outer circle plot

        for i in range(1, len(s2Points), 1):

            # S2 Surface
            s2End = s2Points[i]
            s2Surface = Surface((s2Start[0], s2Start[1]), (s2End[0], s2End[1]), i)

            # Removes some outlier points
            if (s2Surface.segment.length >= (self.surfaceLength - 1)) and (
                s2Surface.segment.length <= (self.surfaceLength + 1)
            ):
                s2Start = s2End
                continue

            # Edge case outliers (mainly for uniform distribution)
            startTuple = (s2Start[0], s2Start[1])
            endTuple = (s2End[0], s2End[1])
            if startTuple not in list(
                self.fullCircle.exterior.coords
            ) or endTuple not in list(self.fullCircle.exterior.coords):
                s2Start = s2End
                continue

            # Vector representations of the normal and the line from the midpoint of S2 to the midpoint of S1
            vNormal = self.vectorHelper(
                (self.normalStart.x, self.normalStart.y),
                (self.normalEnd.x, self.normalEnd.y),
            )
            vS2 = self.vectorHelper(
                (self.midpoint.x, self.midpoint.y),
                (s2Surface.midpoint.x, s2Surface.midpoint.y),
            )

            # Call helper function that uses dot product to calculate angle between vectors
            angle = self.dotProductAngle(vNormal, vS2)
            areaValue, _ = self.intersectionArea(s2Surface)

            dTheta = self.dotProductAngle(self.vLeg1, self.vLeg2)
            pdfArea += areaValue * dTheta

            # Add the point and continue to the next iteration of the loop
            plotPoints.append(Point(angle, areaValue))
            s2Start = s2End

        # Just plotting the outer circle plot points-- no need to normalize
        if comparison == False:
            for point in plotPoints:
                plot_points(point, ax, color="black")

        # Plotting the analytic cosine distribution from the equation y = (1/2pi) * (1 + cos(x)), scaled to be from -pi/2 to pi/2
        # Normalizing the outer circle area value points as well using pdfArea
        else:  # if comparison == True
            cosPoints = []
            adjustedPlotPoints = []
            for plotPoint in plotPoints:
                # Red points (analytic)
                cosXval = plotPoint.x
                cosYval = (1 / math.pi) * (1 + math.cos(cosXval * 2))

                cosPoints.append(Point(cosXval, cosYval))

                # Normalized outer circle points
                adjustedArea = plotPoint.y / abs(pdfArea)
                adjustedPlotPoints.append(Point(plotPoint.x, adjustedArea))

            for adjPoint in adjustedPlotPoints:
                plot_points(adjPoint, ax, color="black")

            for cosPoint in cosPoints:
                plot_points(cosPoint, ax, color="red", marker="1")

        # Generate the plot
        if self.r_offset == 0:
            yUpperLim = plotPoints[-1].y + plotPoints[-1].y * 0.1
            ax.set_ylim(bottom=0, top=yUpperLim)
        plt.show(block=False)

        return

    def analyticUniform(self, ax):
        from shapely import (
            Point,
            plotting,
            Polygon,
            MultiPoint,
            is_closed,
            get_coordinates,
            LineString,
        )
        from shapely.plotting import plot_points
        from matplotlib.pyplot import subplots, ioff
        import matplotlib.pyplot as plt
        import math
        import numpy as np

        ioff()
        fig = subplots()

        self.distributionCircle(0)  # uniform distribution
        self.drawOuterCircle()

        # Creating the plot of angle vs. flux
        s2Points = get_coordinates(
            self.outerCircle
        )  # points defining the outer circle of S2 surfaces
        s2Start = s2Points[0]  # The starting point of the first S2 surface
        plotPoints = []  # Points to plot for the outer circle plot
        pdfArea = 0  # "C" for the outer circle plot

        for i in range(1, len(s2Points), 1):
            # S2 Surface
            s2End = s2Points[i]
            s2Surface = Surface((s2Start[0], s2Start[1]), (s2End[0], s2End[1]), i)

            # Removes some outlier points
            if (s2Surface.segment.length >= (self.surfaceLength - 1)) and (
                s2Surface.segment.length <= (self.surfaceLength + 1)
            ):
                s2Start = s2End
                continue

            # Edge case outliers
            startTuple = (s2Start[0], s2Start[1])
            endTuple = (s2End[0], s2End[1])
            if startTuple not in list(
                self.fullCircle.exterior.coords
            ) or endTuple not in list(self.fullCircle.exterior.coords):
                s2Start = s2End
                continue

            # Vector representations of the normal and the line from the midpoint of S2 to the midpoint of S1
            vNormal = self.vectorHelper(
                (self.normalStart.x, self.normalStart.y),
                (self.normalEnd.x, self.normalEnd.y),
            )
            vS2 = self.vectorHelper(
                (self.midpoint.x, self.midpoint.y),
                (s2Surface.midpoint.x, s2Surface.midpoint.y),
            )

            # Call helper function that uses dot product to calculate angle between vectors
            angle = self.dotProductAngle(vNormal, vS2)  # Plot on x-axis
            areaValue, _ = self.intersectionArea(s2Surface)  # Plot on y-axis

            dTheta = self.dotProductAngle(self.vLeg1, self.vLeg2)
            pdfArea += areaValue * dTheta

            # Add the point and continue to the next iteration of the loop
            plotPoints.append(Point(angle, areaValue))
            s2Start = s2End

        adjustedPlotPoints = []
        for plotPoint in plotPoints:
            # Normalized outer circle points
            adjustedArea = plotPoint.y / abs(pdfArea)
            adjustedPlotPoints.append(Point(plotPoint.x, adjustedArea))

        for adjPoint in adjustedPlotPoints:
            plot_points(adjPoint, ax, color="green")

        # Generate the plot
        plt.show(block=False)

        return

    def showTwoSurfacePlot(self, s2, r_offset=0):
        from shapely import Point, LineString, plotting
        from matplotlib.pyplot import subplots, ioff
        import matplotlib.pyplot as plt

        ioff()
        fig, ax = subplots()
        ax.set_aspect("equal")

        self.distributionCircle(r_offset)
        self.intersectionArea(s2)

        plotting.plot_line(self.segment, ax, color="red", linewidth=2)
        ax.text(self.start.x, self.start.y, "Surface 1", color="red")

        plotting.plot_line(self.normal, ax, color="blue", linewidth=2)
        ax.text(self.normalEnd.x, self.normalEnd.y, "Normal (S1)", color="blue")

        plotting.plot_line(self.dCircleCenter, ax, color="green")
        plotting.plot_polygon(
            self.circle, ax, add_points=False, color="green", linewidth=1
        )
        ax.text(
            self.dCircleCenter.x, self.dCircleCenter.y, self.distType, color="green"
        )

        plotting.plot_polygon(self.triangle, ax, color="orange", linewidth=2)

        plotting.plot_line(s2.segment, ax, color="black", linewidth=2)
        ax.text(s2.midpoint.x, s2.midpoint.y, "Surface 2", color="orange")

        plotting.plot_polygon(
            self.overlapShape, ax, add_points=False, color="black", linewidth=2
        )
        ax.text(
            self.overlapShape.centroid.x,
            self.overlapShape.centroid.y,
            "Overlap Area",
            color="black",
        )

        plt.show(block=False)

        return

    def showOuterCirclePlot(self, r_offset=1):
        from shapely import Point, LineString, plotting
        from matplotlib.pyplot import subplots, ioff
        import matplotlib.pyplot as plt

        ioff()
        fig, ax = subplots()
        ax.set_aspect("equal")

        self.distributionCircle(r_offset)
        self.drawOuterCircle()

        plotting.plot_line(self.segment, ax, color="red", linewidth=2)
        ax.text(self.start.x, self.start.y, "Surface 1", color="red")

        plotting.plot_line(self.normal, ax, color="blue", linewidth=2)
        ax.text(self.normalEnd.x, self.normalEnd.y, "Normal (S1)", color="blue")

        plotting.plot_line(self.dCircleCenter, ax, color="green")
        plotting.plot_polygon(
            self.circle, ax, add_points=False, color="green", linewidth=1
        )
        ax.text(
            self.dCircleCenter.x, self.dCircleCenter.y, self.distType, color="green"
        )

        plotting.plot_polygon(self.outerCircle, ax, add_points=True, color="gray")

        plt.show(block=False)

        return

    # # # # # # # # # # # #
    # # HELPER FUNCTIONS # #
    # # # # # # # # # # # #
    def normalHelper(self, dx, dy, endX, endY, init):
        """Helper function to find the normal direction."""

        if self.end.y > self.start.y:  # +x
            endX += dy
            if self.end.x > self.start.x:
                endY -= dx  # (+x, -y)
            else:
                endY += dx  # (+x, +y)
        else:  # -x
            endX -= dy
            if self.end.x < self.start.x:
                endY += dx  # (-x, +y)
            else:
                endY -= dx  # (-x, -y)
        if init:
            self.normalEndX = endX
            self.normalEndY = endY
        else:
            return (endX, endY)

    def vectorHelper(self, start, end):
        from shapely import Point, LineString
        import numpy as np

        """Finds the vector representation of a segment (surface object). 
            Takes in tuples (x, y) that represent the start and end points of a surface."""

        iSurface = end[0] - start[0]
        jSurface = end[1] - start[1]
        vSurface = np.array([iSurface, jSurface])

        return vSurface

    def dotProductAngle(self, v1, v2):
        import numpy as np

        """Calculates the angle between two vectors using the dot product. 
        v1 and v2 must be np.array objects that represent the surfaces: [iSurface, jSurface]. 
        Use vectorHelper on the start and end points of a surface before passing anything into dotProductAngle."""

        dotProduct = np.dot(v1, v2)
        magnitude1 = np.linalg.norm(v1)
        magnitude2 = np.linalg.norm(v2)

        cosineAngle = dotProduct / (magnitude1 * magnitude2 + 1e-10)

        angle = np.arccos(cosineAngle)

        crossProduct = v1[0] * v2[1] - v1[1] * v2[0]
        if crossProduct < 0:
            angle = -angle

        return angle
