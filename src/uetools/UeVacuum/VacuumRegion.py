class VacuumTests:
 
    def twoSurfacePlot(self):
        '''Plots a source and a receivng surface, including 
            the flux triangle, normal vector, and distribution circle.'''
        S1 = Surface((4, 2), (1, 6), 0)
        S2 = Surface((5, 9), (6, 8), 1)
        S1.showTwoSurfacePlot(S2, r_offset=1)
    
    def outerCirclePlot(self):
        '''Plots the source surface and the outer circle 
            used for obtaining the plot for comparison to 
            analytic distributions.'''
        S1 = Surface((2, 5), (4, 1,), 0)
        S1.showOuterCirclePlot(r_offset=1)

    def analyticPlot(self, ax):
        '''Generates the flux v. angle plot for analytic 
            comparison of the distributions.'''
        S1 = Surface((2, 5), (4, 1), 0)
        S1.showAnalyticPlot(ax, comparison=True, r_offset=1, showBothDist=True)
    
    def UanalyticPlot(self, ax):
        '''Generates the flux v. angle plot for a uniform 
            distribution.'''
        S1 = Surface((2, 5), (4, 1), 0)
        S1.analyticUniform(ax)    

    def trianglePlot(self):
        from shapely import Point, LineString
        import math
        import numpy as np

        '''Plots a triangle geometry.'''
        S1 = Surface((1, 3), (2, 6), 0)
        height = math.sqrt(3) * S1.surfaceLength / 2 # height of the et
        vertex = S1.normal.interpolate(height)
        test = VacuumRegion( [
            (S1.start.x, S1.start.y), 
            (S1.end.x, S1.end.y), 
            (vertex.x, vertex.y)
        ], multiprocess=False)
        f = test.plotGeometry(labels=True, showCircle=True)

    def squarePlot(self):
        from shapely import Point, LineString

        '''Plots a square geometry.'''
        S1 = Surface((1, 2), (3, 6), 0)
        # # # Make the sides of the square perpendicular to self # # #
        side1Start = (S1.end.x, S1.end.y)
        side1End = S1.normalHelper(abs(S1.dx), abs(S1.dy), side1Start[0], side1Start[1], False)

        side3End = (S1.start.x, S1.start.y)
        side3Start = S1.normalHelper(abs(S1.dx), abs(S1.dy), side3End[0], side3End[1], False)

        test = VacuumRegion([side1Start, side1End, side3Start, side3End], multiprocess=False)
        f = test.plotGeometry(labels=True, showCircle=True)

    def shadedSquarePlot(self):
        from shapely import Point, LineString

        """Plots a shaded square geometry. 
            To be used alongside the source surface S1-- Surface((2, 1), (1, 1))-- 
            which is defined in the lineOfSightPlot function in the test functions."""

        geometryVertices = [(3, 1), (1, 1), (1, 5), (2, 6), (1, 7), (7, 7), (7, 1), (5, 1), (5, 3), (3, 3)]

        test = VacuumRegion(geometryVertices, P=1, multiprocess=False)
        test.saveVacuumRegion("SavedVacuumRegion")
        # for i in test.errors:
        #     f = test.plotGeometry(labels=True, testsurf=i)
        #     test.surfaces[i].plotSelf(ax=f)
        #     test.surfaces[i].printReport()
        f = test.plotGeometry(labels=True, testsurf=6, showCircle=True)
        return test

    # def tokamakPlot(self, savefile): # Functionality transferred over to VNM_interface.couple
    #     from uetools import Case
    #     from numpy import zeros
    #     '''Plots the full tokamak geometry.'''
    #     c = Case(savefile, inplace=True)
    #     (main, pf) = c.coupling.get_snull_vacuum_regions(maxlength = 0.0087)
    #     # nobug = zeros((main[0].shape[0]-1, main[0].shape[1]))
    #     # test = VacuumRegion(main[0], P=main[1]) # main geometry
    #     test = VacuumRegion(pf[0], P=pf[1] - 1) # private flux region
    #     '''To plot surfaces that aren't meeting unity:'''
    #     # for i in test.errors:
    #     #     if (i > 50) and (i<90):
    #     #         f = test.plotGeometry(labels=False, testsurf=i, markers='.')
    #     #         f.get_axes()[0].set_title(f"Surface {i}")
    #     f = test.plotGeometry(labels=False, testsurf=27, showCircle=True)
    #     m = test.heatmapPlot() # TO PLOT MATRIX HEATMAPS
    #     # f = test.plotGeometry(labels=False, testsurf=150, showCircle=True)
    #     # f = test.plotGeometry(labels=False, testsurf=4)
    #     return test


class VNM_interface:
    def __init__(self, case, vnm_setup):
        from numpy import load, array
        from collections import defaultdict
        if case.get('geometry')[0].strip().decode('UTF-8') not in ['snull', 'dnull']:
            raise Exception("VNM model only implemented for singe nulls geometries!")
        self.coupling = case.coupling
        self.set = case.setue
        self.populate = case.populate
        self.info = case.info
        self.tools = case.tools
        self.getue = case.getue
        self.hdf5search = case.tools.hdf5search
        self.uevars = {
            'inner': [
                'isvacuummodelpf',
                'cfteleoutpf',
                'cftelematrixpf',
                'fngyi_use',
            ],
            'outer': [
                'isvacuummodelw',
                'cfteleoutw',
                'cftelematrixw',
                'fngyo_use',
            ]
        }
    
        # TODO: Test HDF5 restore with UETOOLS HDF5 saves

        def read_txt(file):
            nodes = []
            with open(file) as f:
                for line in f:
                    nodes.append(
                            tuple(
                                [float(x) for x in line.replace(',',' ').strip().split(' ')]
                    ))
            return nodes

        # Parse VNM setup - determine whether to pass to restore or generate
        if vnm_setup is not None:
            if not isinstance(vnm_setup, dict):
                raise TypeError("vnm_setup must be dict")
            if "regions" not in vnm_setup:
                raise Exception("At least one region must be defined")
            # Check that regions are present and satifsy requirements
            for region in vnm_setup['regions']:
                if "location" not in region:
                    raise Exception("Region must have 'location' set to 'inner'/'outer'")
                else:
                    if region['location'].lower() not in ['inner', 'outer']:
                        raise Exception("Region 'location' must be 'inner'/'outer'")
                if 'isvacuummodel' not in region:
                    raise Exception("Specify VNM model for species using isvacuummodel")
            regions = vnm_setup.pop('regions')
            # Detect whether to generate, restore surfaces, or restore BCs
            groups = defaultdict(list)
            for region in regions:
                mode = region.get("mode", "generate")
                groups[mode].append(region)
            generate = groups["generate"]
            restore_surfaces = groups["restore_surfaces"]
            restore_boundary = groups["restore_boundary"]
            for region in restore_surfaces:
                if 'surface_file' not in region:
                    raise KeyError(f"'surface_file' to restore not specified for region {region['name']}")
            for region in restore_boundary:
                if 'boundary_file' not in region:
                    
                    print(f"'boundary_file' to restore not specified for region {region['name']}, using save file '{self.info['savefile']}'")
                    region['boundary_file'] = self.info['savefile']
                # Restore regions from save file as requested
                self.restore(region['boundary_file'], self.uevars[region['location']], region['name'])
#            # Split into generated/restored regions
#            restore = [region for region in regions if region.get('restore')]
#            generate = [region for region in regions if not region.get('restore')]
            # Check generated regions satisfy conditions
            for region in generate + restore_surfaces:
                if 'nodes' in region:
                    if "savefile" in region:
                        raise Exception("Either specify save file or node list" +
                        f" for {region['name']}, not both!")
                    if isinstance(region['nodes'], str):
                        try:
                            region['nodes'] = load(region['nodes'])
                        except:
                            region['nodes'] = array(read_txt(region['nodes']))
                if 'pump' in region:
                    for pump in region['pump']:
                        if isinstance('nodes', str):
                            try:
                                pump['nodes'] = load(pump['nodes'])
                            except:
                                pump['nodes'] = array(read_txt(pump['nodes']))
                        else:
                            # TODO: Assert pump nodes are OK
                            1 
    
        regions =   [x.copy() for x in generate] + \
                    [x.copy() for x in restore_boundary] + \
                    [x.copy() for x in restore_surfaces]
        if len(generate) > 0: 
            self.generate(generate, **vnm_setup)
        if len(restore_surfaces) > 0: 
            self.generate(restore_surfaces, restore_surfaces=True, **vnm_setup)
        vnm_setup['regions'] = regions
        self.populate()

    def restore(self, save_file, restore_vars=None, location=''):
        """Restores existing telematrices from the provided save file, and calculates puffing input arrays if desired.
        Uses the current save file if none is provided.
        """
        import h5py
        import warnings
        import numpy
        from uedge import com, bbb

        if restore_vars is None:
            restore_vars = self.uevars['inner'] + self.uevars['outer']
        if save_file is None:
            save_file = self.info['savefile']
    
        for var in restore_vars:
            val = self.hdf5search(save_file, var)
            if val is None:
#                raise Exception(f"Variable '{var}' not found in {save_file}!")
                # Assume defaults used and not thus written to save
                pass
            else:
                self.set(var, val)
        self.populate()
        print(f"Successfully restored {location} VNM from {save_file}")
        return



    def generate(self, regions, restore_surfaces=False, maxlength=0.01, plot=False):
        ''' 
            region_setup,
            write = False,
            maxlength = 0.01,
                sol=False, 
                sol_savename=None, 
                sol_nodes=None, 
                sol_hdf5location=None, 
                sol_puff=None, 
                sol_pump=None, 
                vnm_sol=True, 
                pfr=False, 
                pfr_pump=None, 
                pfr_savename=False, 
                plot_setup=False, 
                pfr_nodes=None,
                pfr_hdf5location=None,
                pfr_puff=None,
                vnm_pfr=True, 
                kwargs_sol={}):
        '''
        """Generates telematrices for given VacuumRegions.
        
                  Keyword arguments:
        sol - dict/None/False (default = False)
            Setup for main-SOL, definingt the main-SOL model. If False,
            VNM is not used for SOL. If None, SOL VNM is generated from UEDGE
            data. If dictionary, VNM is created based on the dictionary settings.
            Dictionary keys available:
                savefile - path to pickle/HDF5 file containing SOL VNM save
                hdf5location - (default: vnm/sol) 
                        string pointing to the location of the SOL VNM setup in 
                        the HDF5 if savefile is an HDF5
                picklename - name of file where VNM for SOL is pickled
                pump_setup - nested setup dictionary for SOL pumping surface
                    Pump setting keys:
                    name - defines the pump key name, contains dict with follwoing keys:
                        type - defines the pump setup type. Available options "region"
                        region options:
                        nodes (required) - nodes defining polygon of region: must define 
                                open plygon shape, minimum 3 node (x,y) pairs
                        recycling (required) - recycling coefficient for surfaces intersecting
                                with pump region 
                puff_setup - nested setup dictionary for puff setup type. Available options "point"
                    Puff settings keys:
                    name - defines the puff key name, contains the following keys
                    type - defines the puff type. Available options: 'point'
                        point options:
                        location - (R, Z) coordinate of puff. Puff automatically assigned
                            to the surface closest to location.
                        current - puff strength in part/s
                nodes - list of SOL nodes to replace the automatically generated ones
                    from the UEDGE case
                recycling - recycling coefficient for SOL (default = 1)
        pfr - ditto for the PFR vacuum region
        write -- decides whether or not to write generated matrices into save file (default False).
        """

        import h5py
        from numpy import zeros, hstack, vstack, pad
        from uedge import com, bbb
        self.populate(verbose=False)
        self.regions = []


        self.nx = self.getue('nx')
        self.ixpt1 = self.getue('ixpt1')[0]
        self.ixpt2 = self.getue('ixpt2')[0]
        self.ngsp = self.getue('ngsp')  
        self.nx = self.getue('nx')
        (sol_nodes, pfr_nodes) = self.coupling.get_snull_vacuum_regions(maxlength=maxlength)
        self.nodes = {'inner': pfr_nodes, 'outer': sol_nodes}
        self.regions = {}
        self.output = {}
        dimension = self.nx+2
        for region in regions:
            if 'name' in region:
                name = region.pop('name')
            else:
                name = len(regions)+1
            isvacuummodel = region.pop('isvacuummodel')
            # Perform inner/outer setup
            location = region.pop('location').lower()
            (vnm, P) = self.nodes[location]
            if location == 'inner':
                P = self.ixpt1 + (self.nx - self.ixpt2)
                savekey = ('pf', 'i')
            else:
                P = self.nx
                savekey = ('w', 'o')
            self.output[name] = {
                'telematrix': zeros((dimension, dimension, 6)),
                'puff': zeros((dimension, self.ngsp))
            }
            if restore_surfaces:
                print(f"Restoring vacuum region '{name}' from pickle/HDF5")
                vnm = region.pop('surface_file')
            elif 'nodes' in region:
                print(f"Generating vacuum region '{name}' from user-defined nodes")
                vnm = region.pop('nodes')
            else:
                print(f"Generating new vacuum region '{name}'")
            self.regions[name] = VacuumRegion(vnm, P=P, **region)
                 
            if location == 'inner':
                # Expand PF matrix  along core cut: get dimensions and mismatch
                dim_expand = dimension - P - 2
                tele = pad(self.regions[name].telematrix.transpose(), pad_width=1)
                # Expand PF matrix along vertical axis
                tele = vstack([
                    tele[:self.ixpt1+1],
                    zeros((dim_expand, P + 2)), 
                    tele[self.ixpt1+1:]   
                ])
                # Expand PF matrix along horizontal axis
                tele = hstack([
                    tele[:, :self.ixpt1+1],
                    zeros((dimension, dim_expand)),
                    tele[:, self.ixpt1+1:],
                ])
                for j in range(6):
                    self.output[name]['telematrix'][:,:,j] = tele
                # Expand and pad PF puffing array
                self.output[name]['puff'] = vstack([
                    zeros((1,6)),
                    self.regions[name].puff_vector[:self.ixpt1],
                    zeros((dim_expand, 6)),
                    self.regions[name].puff_vector[self.ixpt1:],
                    zeros((1,6)),
                ])        
            else:
                # Populate local pump
                self.output[name]['puff'] = vstack([
                        zeros((1,6)),
                        self.regions[name].puff_vector,
                        zeros((1,6)),
                    ])  
                for j in range(6):
                    self.output[name]['telematrix'][:,:,j] = pad(
                                    self.regions[name].telematrix.transpose(), 
                                    pad_width=1
                    )
            # Populate UEDGE cftelematrixw array
            self.set(f'cftelematrix{savekey[0]}', self.output[name]['telematrix'])
            # Turn on the VNM model in UEDGE
            if isinstance(isvacuummodel, dict):
                for key, value in isvacuummodel.items():
                    if isinstance(value, int):
                        self.getue(f'isvacuummodel{savekey[0]}', cp=False)[key] = value
                    elif isinstance(value, list):
                        listlen = len(value)
                        self.getue(f'isvacuummodel{savekey[0]}', cp=False)[key:key+listlen] = value
            else:
                self.set(f'isvacuummodel{savekey[0]}', isvacuummodel)
            self.set(f'cfteleout{savekey[0]}', 1.0)
            # Populate the puffing array
            self.set(f'fngy{savekey[1]}_use', self.output[name]['puff'][:,:self.ngsp])

        # Plot setup if requested
        if plot:
            self.plot_grid(pf_test_surf=[], label=False)

    def plot_grid(self, sol_plot=True, pfr_plot=True, sol_test_surf=[], pf_test_surf=[], label=False, **kwargs):
        from matplotlib.pyplot import subplots
        f, ax = subplots()
        for regionname, region in self.regions.items():
            region.plotGeometry(labels=label, ax=ax, **kwargs)


    def save_hdf5(self, file, **kwargs):
        ''' Saves VNM data to HDF5 '''
        from h5py import File
        if not isinstance(file, File):
            raise TypeError("file must be an open HDF5 File object")
        vnm = file.require_group('vnm')
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
    def __init__(self, nodeList, P=0, r_offset_plasma=1, r_offset_material=1, multiprocess=True, ncores=None, verbose=True, material_recycling=1, pump=None, puff=None, reflections=1e6, hdf5location="vnm", savename=None, isvacuummodel=None, write=False, **kwargs):
        """
        nodeList - str, list of nodes, or HDF5 file name
                HDF5 - populates data based on hdf5location pointing to the vnm setup in
                    the HDF5 file
                list of nodes - generates from scratch based on input
                str - reads from pickle and populates based on setup 
        P - number of plasma surfaces in region: leads arrays/matrices
        r_offset_plasma - offset of distribution circle relative to circle radius. 
            1 - cosine
            0 - uniform
        r_offset_material - ditto, for material surfaces
        multiprocess - Multiprocessing when constructing objects
        ncores - multiprocessing cores
        verbose - False supresses output
        material_recycling - recycling coefficient on material surfaces
        pump_setup - dictionary defining pumping setup
            Nested dict of pumping regions. Region name and properties. Sorted by
            order of appearance, later regions overwrite earlier ones.
            The required structure is:
            pump_setup[region_name] = {
                'type': "region",
                'recycling': recycling_coefficient,
                'nodes': nodelist
            }
            The required "data" entry is determined by the available types, 
            listed below:
                "region" - A closed polygon is created based on the supplied
                    nodes. All polygons intersecting with the polygon are
                    assigned a reccyling coefficient as defined by 'recycling'.
                "recycling" - recycling coefficient to be applied to intersecting
                    surfaces
                "nodes" - nested list of (X,Y) nodes that define the pumping surface
        puff_setup - dictionary defining puffing surfaces 
        reflections - the number of reflections to be considered: the exponent of the TMM 
        """
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
            if is_hdf5(nodeList): # Restoring from HDF5 file
                with File(nodeList, 'r') as f:
                    vnm = f[hdf5location]
                    for var in ['nodeList', 'P', 'r_offset_material', 'r_offset_plasma', 'reflections']:
                        self.__setattr__(var, vnm[var][()])
                    if 'puff' in vnm:
                        setup = {}
                        for puffname in vnm['puffs'].keys():
                            puff[puffname] = {}
                            for var in vnm['puffs'][puffname].keys():
                                puff[puffname][var] = vnm['puffs'][puffname][var][()]
                    if 'pump' in vnm:
                        pump = {}
                        for pumpname in vnm['pumps'].keys():
                            pump[pumpname] = {}
                            for var in vnm['pumps'][pumpname].keys():
                                pump[pumpname][var] = vnm['pumps'][pumpname][var][()]

        starttime = time()
                
        if isinstance(self.nodeList, str): 
            # TODO: Check whether requested file exists and is pickle
            with open(self.nodeList, 'rb') as f:
                save = load(f)
                self.surfaces = save['surfaces']
                self.P = save['P']
                self.nodeList = save['nodeList']
            P=self.P

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
                if i < self.P: # Plasma surfaces
                    self.surfaces[i] = Surface((startNode.x, startNode.y), (endNode.x, endNode.y), i, r_offset=r_offset_plasma)
                else: # non-plasma surfaces
                    self.surfaces[i] = Surface((startNode.x, startNode.y), (endNode.x, endNode.y), i, r_offset=r_offset_material)
                
            # Create Polygon of Vacuum region for intersect checks
            self.geometry = Polygon(self.nodeList) 

            if multiprocess:
                if ncores is None:
                    ncores = cpu_count()
                else:
                    ncores = min(cpu_count(), ncores)
                environ['UETOOLS_SILENT'] = "1"
                parent_conn, child_conn = Pipe()
                manager = Manager()
                surface_chunks = manager.dict()
                # Create list of surfaces to be calculated by each subprocess
                sublist = array_split(array(range(len(self.surfaces))), ncores)
                # Spawn subprocesses
                subprocesses = []
                print(f"Calculating {len(self.surfaces)} surface couplings on {ncores} threads...")
                for subprocess in sublist:
                    subprocesses.append(Process(
                        target=self.subprocess_execute,
                        args=(surface_chunks, child_conn, list(subprocess), self.surfaces, self.geometry),
                        kwargs=({"verbose": verbose})
                    ))
                    subprocesses[-1].start()
                for subprocess in subprocesses:
                    subprocess.join()
                self.surfaces = deepcopy(surface_chunks)
                environ['UETOOLS_SILENT'] = "0"

            else:
                # Iterate surfaces to identify surface neigbors
                for _, surface in tqdm(self.surfaces.items()):
                    surface.getNeighbors(self.surfaces, self.geometry)

        self.time = time() - starttime
        self.numSurfaces = len(self.surfaces)

        # Dictionary of surface reflection coefficients
        """ Set up recycling coefficients """
        self.R_array = zeros(self.numSurfaces)
        self.R_array[self.P:] = self.material_recycling

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
        self.puff_vector = zeros((self.P,6))
        if puff is not None:
            self.puff = puff
            for _puff in puff:
                self.puffs.append(self.create_puff(_puff))
                self.puff_vector[:,_puff['igsp']] += self.puffs[-1]['puff_vector'].flatten()
        
        """ Save to pickle if requested """
        if write:
            self.saveVacuumRegion(savename)

        '''Print statements to use if surfaces are not conserving flux via line of sight.'''
        # if not self.checkContinuity(False): # BRING BACK AFTER TESTING
        #     print("Warning! Continuity violated for surfaces:", self.errors)
        #     print(f"Fluxes: {[(s, self.surfaces[s].totflux) for s in self.errors]}")

    @staticmethod
    def subprocess_execute(output, conn, surflist, surfaces, geometry, verbose=True):
        from os import getpid
        count = 0
        msg = [.25, .50, .75, 1]
        # Iterate Process surfaces
        for surfid in surflist:
            if verbose:
                if count/len(surflist) > msg[0]:
                    print(f"Process {getpid()} {msg[0]*100}% completed.")
                    msg.pop(0)
            surfaces[surfid].getNeighbors(surfaces, geometry)
            output[surfid] = surfaces[surfid]
            count += 1
        if verbose:
            print(f"Process {getpid()} completed surfaces {surflist[0]}-{surflist[-1]}.")

    def matrixPower(self, matrix, power):
        '''Raises a given matrix to the specified power.'''
        from numpy import zeros, identity
        from scipy.sparse import csr_array, block_array, linalg

        resultMatrix = linalg.matrix_power(matrix, power)

        return resultMatrix


    def matrices(self):
        '''Creates R (self.R_matrix), C (self.C_matrix), A (self.A_matrix), 
            B (self.B_matrix), and AB (self.AB_matrix) matrices.'''
        import numpy
        from numpy import zeros, identity, percentile, log, diag
        from scipy.sparse import csr_array, block_array
        import seaborn as sns
        import matplotlib.pyplot as plt

        # Array representations of R and C
        self.C_array = zeros((self.numSurfaces, self.numSurfaces))
        # Populate C array and take transpose
        for surfaceID, surface in self.surfaces.items(): # self.surfaces.items()
            for outputID in surface.neighbors.keys():
                self.C_array[surfaceID][outputID] = surface.neighbors[outputID]['flux']
        self.C_array = self.C_array.transpose()
        # R and C into sparse matrices
        self.R_matrix = csr_array(diag(self.R_array))
        self.C_matrix = csr_array(self.C_array)
        # Zero and identity sparse matrices
        Zero_matrix = csr_array(zeros((self.numSurfaces, self.numSurfaces)))
        Identity_matrix = csr_array(identity(self.numSurfaces))
        # Create A, B, and AB sparse matrices
        self.A_matrix = block_array([[self.C_matrix, Zero_matrix], [Zero_matrix, Identity_matrix]])
        self.B_matrix = block_array([[self.R_matrix, Zero_matrix], [Identity_matrix - self.R_matrix, Identity_matrix]])
        # A * B
        self.AB_matrix = self.A_matrix @ self.B_matrix 
        # (AB)^M * A
        self.AB_power_A = self.matrixPower(self.AB_matrix, self.reflections) @ self.A_matrix 

    def createTeleMatrix(self):
        from numpy import zeros, transpose

        # Final transport matrix
        self.telematrix = zeros((self.P, self.P))

        gamma_array = zeros((self.numSurfaces * 2, 1))
        for i in range(0, self.P):
            gamma_array[i, 0] = 1
            rowCalculation = self.AB_power_A @ gamma_array
            gamma_array[i, 0] = 0
            gammaOut = rowCalculation[self.numSurfaces:]
            gammaFinal = gammaOut[0:self.P].flatten()
            for j in range(self.P):
                self.telematrix[j, i] = gammaFinal[j]

    def create_puff(self, puff_setup):
        from shapely import Point
        from numpy import argsort, zeros
        if "type" not in puff_setup:
            raise KeyError(f"Define a puff type")
        ret = {'type': puff_setup["type"]}
        
        if puff_setup['type'] == "point":
            for key in ['location', 'current', 'igsp']:
                if key not in puff_setup:
                    raise KeyError("Required 'point' puff setup entry"+
                        f" '{key}' not found.")
            try:
                ret['point'] = Point(puff_setup['location'])
            except Exception:
                raise  
            nodes = list(self.geometry.exterior.coords)[self.P:-1]
            dists = [self.P+ret['point'].distance(Point(n)) for n in nodes]
            ret['material_surface_index'] = argsort(dists)[:2].max() + self.P
            ret['current'] = puff_setup['current']
            ret['igsp'] = puff_setup['igsp']
            ret['location'] = puff_setup['location']
            drive = zeros((self.numSurfaces*2,1))
            drive[ret['material_surface_index'], 0] = ret['current']
            ret['puff_vector'] = (self.AB_power_A @ drive)[self.numSurfaces:][:self.P]
            
        else:
            raise KeyError(f"Puff type '{pump_setup['type']}' not recognized!" +
                "\nAvailable options are: 'point'")
        return ret

    def create_pump_region(self, pump_setup):
        from shapely import Polygon, intersects, Point
        if "type" not in pump_setup:
            raise KeyError(f"Define a pump type")
        ret = {"type": pump_setup['type']}

        if pump_setup['type'] == "region":
            for key in ['nodes', 'recycling']:
                if key not in pump_setup:
                    raise KeyError("Required 'region' pump setup entry"+
                        f" '{key}' not found.")
            if len(pump_setup['nodes'])<3:
                raise AttributeError("Too few nodes provided: provide "+
                    "at least three nodes for Polygon!")
            ret["polygon"] = Polygon(pump_setup['nodes'])
            ret['nodes'] = pump_setup['nodes']
            ret["pumped_segments"] = []
            ret["recycling"] = pump_setup["recycling"]
            for i in range(self.P, self.numSurfaces):
                if intersects(self.surfaces[i].segment, ret["polygon"]):
                    self.R_array[i] = ret["recycling"]
                    ret['pumped_segments'].append(i) 
        else:
            raise KeyError(f"Pump type '{pump_setup['type']}' not recognized!" +
                "\nAvailable options are: 'region'")
        return ret

    def saveVacuumRegion(self, savename):
        from pickle import dump
        '''Use to save a Vacuum Region to avoid having to generate a new one every time. 
            Be sure to set pf/main (tokamakPlot), r_offset (Surface constructor), and 
            variation (Vacuum Region constructor).'''

        save = {
            'surfaces': self.surfaces,
            'P': self.P,
            'nodeList': self.nodeList
        }
        with open(savename, 'wb') as f:
            dump(save, f)

    def writeVacuumSetup(self, obj, write_matrices=True):
        from h5py import File, Group

        def save_csr(group, mat):
            for var in ['data', 'indices', 'indptr']:
                if var in group:
                    del group[var]
                group.create_dataset(var, data=mat.__getattribute__(var))
            group.attrs['shape'] = mat.shape


        variables = [
            "nodeList",
            "P",
            "r_offset_plasma",
            "r_offset_material",
            "material_recycling",
            "reflections",
            "telematrix"
        ]
        csr = [
            "R_matrix",
            "C_matrix",
            "A_matrix",
            "B_matrix",
            "AB_matrix",
            "AB_power_A",
        ]
        # TODO: Storing CSR matrices?

        if isinstance(obj, str):
            f = File(obj, 'a')
            vnm = f.require_group('vnm')
        elif isinstance(obj, Group):
            vnm = obj
        else:
            raise TypeError("Object must be Group or file name as string where to write.")

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
#        for pumpname, setup in self.pumping_regions.items():
            pump = pumps.require_group(pumpname)
            for var, data in setup.items():
                if var not in ["polygon"]:
                    if var in pump:
                        del pump[var]
                    pump.create_dataset(var, data=data)
            pumpid += 1
        # Stor puff setups
        if len(self.puffs) > 0:
            puffs = vnm.require_group("puffs")
        puffid = 1
        for setup in self.puffs:
            puffname = f"puff_{puffid}"
            puff = puffs.require_group(puffname)
            for var, data in setup.items():
                if var not in ['point']:
                    if var in puff:
                        del puff[var]
                    puff.create_dataset(var, data=data)    
    
        if isinstance(obj, str):
            f.close()

    def checkContinuity(self, verbose=True):
        '''For identifying errors in flux unity for surfaces in the geometry.'''
        self.errors = []
        for surfid, surface in self.surfaces.items():
            if abs(surface.totflux - 1) > 1e-6:
                if verbose:
                    surface.printReport()
                self.errors.append(surfid)
        return len(self.errors)==0

    def plotGeometry(self, ax=None, labels=False, testsurf=[], 
        showCircle=False, markers=None, connectionLineWidth=0.5,**kwargs):
        from matplotlib.pyplot import subplots, Figure, Axes, ioff
        import matplotlib.pyplot as plt


        '''Plot the full geometrical representation of the VacuumRegion.'''

        if isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        elif ax is None:
            f, ax = subplots(figsize=(5,10))
        elif not isinstance(ax, Axes):
            raise Exception('ax not a valid Figure or Axes object')
        
        if isinstance(testsurf, int):
            testsurf = [testsurf]
        
        ioff()
            
        for surfid, surface in self.surfaces.items():
            color = 'k'
            if (surface.ID <  self.P):
                color = 'red'
            else:
                if self.R_array[surfid] != self.material_recycling:
                    color='grey'
            surface.plotSelf(color=color, ax=ax, label=labels, showCircle=showCircle)

        for itest in testsurf:
            self.surfaces[itest].plotConnections(
                            ax=ax, 
                            linewidth=connectionLineWidth,
                            **kwargs
            )
        c = 0
        N_regions = len(self.pumping_regions) + len(self.puffs)
        colors = [plt.get_cmap("rainbow")(i/N_regions) for i in range(N_regions)]
        for region in self.pumping_regions:
            for i in region['pumped_segments']:
                self.surfaces[i].plotSelf(color=colors[c], ax=ax, showCircle=False)
            c += 1
        for puff in self.puffs:
            ax.plot(puff['point'].xy[0][0], puff['point'].xy[1][0], 'o', color=colors[c])
            self.surfaces[puff['material_surface_index']].plotSelf(color=colors[c], ax=ax, showCircle=False)
            c += 1
    
        


        for line in ax.lines:
            line.set_marker("")

        ax.set_aspect('equal')
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
        '''Creates a heatmap of C, R, and Transport matrices.'''

        '''Print statements to check for unity of transport matrix.'''

        # Plotting heatmaps of C, R, and Output (Transport)
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        fig.suptitle("Cosine Distribution", fontsize=16)

        matricesToPlot = [self.C_array, self.R_array, self.telematrix]
        matricesToPlotNames = ["C", "R", "Output"]
        for i, ax in enumerate(axes):
            sns.heatmap(matricesToPlot[i], cmap='jet', annot=False, ax=ax, norm=LogNorm(vmin=1e-5, vmax=1))
            ax.set_title(matricesToPlotNames[i])
            ax.set_aspect('equal')
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

    def __init__(self, start, end, ID, material=1, emitting=0, absorbing=0, r_offset=1):
        from shapely import Point, LineString, plotting
        from matplotlib.pyplot import subplots
        import math

        '''Reference Variables/Important:
            self.start (Point), self.end (Point), self.ID, self.segment (LineString), 
            self.surfaceLength, self.midpoint (Point), self.normalStart (Point), 
            self.normalEnd (Point), self.normal (LineString).'''

        '''Start and end passed into the constructor are tuples (x, y)'''

        # Start and end points of the surface and a segment representation of the surface 
        self.start = Point(start[0], start[1])
        self.end = Point(end[0], end[1])

        # LineString representing the surface
        self.segment = LineString([start, end]) 
        self.ID = ID 

        self.surfaceLength = math.sqrt(
                                (self.end.x - self.start.x)**2 \
                                + (self.end.y - self.start.y)**2
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
                abs(self.dx), 
                abs(self.dy),
                self.normalEndX, 
                self.normalEndY, 
                True
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

        self.epsilon = 1e-5 # use as a reference for buffering works when epsilon = 0.00001 (1e-5) --> Use for adjusting line of sight

        return

    def distributionCircle(self, r_offset): # creates the distribution circle
        from shapely import Point, LineString, plotting, LinearRing, Polygon
        from matplotlib.pyplot import subplots
        import math

        """Reference Variables/Important:
        self.dCircleCenter (Point), self.circle (Polygon), self.r_offset
        """

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
    
        if self.end.x < self.start.x: # builds circle clockwise rather than counter clockwise
            for i in range(numPoints, -1, -1):
                # Calculates every angle from 0-2pi, placing 500 points to create the circle
                angle = (2 * math.pi) * (i / numPoints) 
                x = self.dCircleCenter.x + radius * math.cos(angle)
                y = self.dCircleCenter.y + radius * math.sin(angle)
                if angle == math.pi and self.distType == "Cosine Distribution":
                    circlePoints.append((self.midpoint.x, self.midpoint.y))
                    circlePoints = sorted(circlePoints, key=lambda x: x[0], reverse = True)
                    i += 1
                circlePoints.append((x, y))
        else: 
            for i in range(numPoints + 1):
                angle = (2 * math.pi) * (i / numPoints)
                x = self.dCircleCenter.x + radius * math.cos(angle)
                y = self.dCircleCenter.y + radius * math.sin(angle)
                if angle == math.pi and self.distType == "Cosine Distribution":
                    circlePoints.append((self.midpoint.x, self.midpoint.y))
                    circlePoints = sorted(circlePoints, key=lambda x: x[0], reverse = True)
                    i += 1
                circlePoints.append((x, y))

        # Object representation of the distribution circle
        self.circle = Polygon(circlePoints) 

        return
    
    def printReport(self):
        print("Surface {}: {}".format(self.ID, self.totflux))
        for neighid, neighbor in self.neighbors.items():
            print("    -> {}: {}".format(f"{neighid}".rjust(4), neighbor['flux']))

    def intersectionArea(self, s2):
        from shapely import Point, LineString, plotting, Polygon, is_closed
        import math

        '''Finds the overlapping area (flux) between two surfaces, given that self has a distribution circle generated.
            Creates/draws the relevant shapes for finding the flux (fractional area) and other reference.'''

        '''Reference Variables/Important:
        self.triangle, self.leg1, self.leg2, self.overlapShape, overlapArea, fractionalArea'''

        if self.circle == None:
            print("Call distributionCircle on Surface before finding intersection area!")
            return

        # Triangle of flux
        self.triangle = Polygon([s2.start, s2.end, self.midpoint, s2.start]) # from midpoint of self to the endpoints of s2

        # Legs of the triangle
        self.leg1 = LineString([self.midpoint, s2.start])
        self.leg2 = LineString([self.midpoint, s2.end])

        # Overlap of triangle and distribution circle
        self.overlapShape = self.triangle.intersection(self.circle)

        # Vector representations of triangle legs
        self.vLeg1 = self.vectorHelper((self.midpoint.x, self.midpoint.y), (s2.start.x, s2.start.y))
        self.vLeg2 = self.vectorHelper((self.midpoint.x, self.midpoint.y), (s2.end.x, s2.end.y))
        
        # Getting the correct area of the distribution circle on one side of the normal line 
        overlapArea = self.overlapShape.area
        if (0 <= self.r_offset and self.r_offset < 1) and self.circle.intersects(self.segment): # A circle shifted between uniform and cosine
            buffLine = self.segment.buffer(self.epsilon * 1e-7)
            splitdCircle = self.circle.difference(buffLine) # Split the distribution circle and the surface line
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
        from shapely import intersects, difference, crosses, buffer, contains, intersection

        #  Fractional Area Calculations
        for neighid, neighbor in surfaces.items():

            if self.ID != neighid:
                flux, triangle = self.intersectionArea(neighbor)
                    
                # Flux would go to the wrong (back) side of the surface
                try:
                    if not intersects(triangle, neighbor.normal) or triangle.intersection(neighbor.normal) == neighbor.midpoint:
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
                if AllIntersections.geom_type == 'Polygon':
                    if not self.getSmallestIntersectAngle(neighbor, geometry, AllIntersections):
                        continue

                # Multiple intersections with the outside region and the triangle/triangle legs
                elif AllIntersections.geom_type == 'MultiPolygon' or AllIntersections.geom_type == 'GeometryCollection':
                    # Check all possible polygons
                    for polygon in AllIntersections.geoms:
                        if not self.getSmallestIntersectAngle(neighbor, geometry, polygon):
                            continue
                        
                # Set new S2 Surface after line of sight adjustments and calculate the new flux. Update the total flux out of Surface 1. 
                newS2 = Surface(
                                (self.leg1SmallestPoint[0], self.leg1SmallestPoint[1]), 
                                (self.leg2SmallestPoint[0], self.leg2SmallestPoint[1]), 
                                neighbor.ID
                )
                flux, triangle = self.intersectionArea(newS2)
                self.totflux += flux
                if flux > 0:
                    if neighid not in self.neighbors:
                        self.neighbors[neighid] = {}
                    self.neighbors[neighid]['flux'] = flux
                    self.neighbors[neighid]['los'] = triangle # Get the new triangle shape and add it here
    
    def getSmallestIntersectAngle(self, neighbor, geometry, polygon):
        from shapely import intersects, crosses, contains, buffer, intersection, Point 

        '''Adjust the legs of the flux triangle in accordance with flux calculations.'''

        if contains(buffer(geometry, self.epsilon), polygon): # When polygon is on the border and not at all outside the geometry
            return True

        if crosses(polygon, buffer(self.leg1, self.epsilon * 1e-1)) and crosses(polygon, buffer(self.leg2, self.epsilon * 1e-1)):
            return False

        # If Leg 1 (from S1 midpoint to S2 start) intersects the outside region
        if (    intersects(polygon, buffer(self.leg1, self.epsilon * 1e-1)) \
                and (contains(geometry, self.leg1) == False)
        ):

            if polygon.geom_type == 'LineString':
                coordinates1 = list(polygon.coords)
            else:
                coordinates1 = list(polygon.exterior.coords)

            leg1Start = buffer(self.midpoint, self.epsilon) 
            leg1End = buffer(neighbor.start, self.epsilon)

            if (    (contains(leg1Start, intersection(polygon, buffer(self.leg1, self.epsilon * 1e-2))) == False) \
                    and (contains(leg1End, intersection(polygon, buffer(self.leg1, self.epsilon * 1e-2))) == False)
            ):
                for pair in coordinates1:
                    if Point(pair) == self.midpoint or contains(buffer(self.midpoint, self.epsilon * 1e2), Point(pair)):
                        continue
                    self.vNewLeg1 = self.vectorHelper(
                                (self.midpoint.x, self.midpoint.y),
                                pair
                    )
                    self.vNewLeg2 = self.vectorHelper(
                                (self.midpoint.x, self.midpoint.y), 
                                (self.leg2SmallestPoint[0], self.leg2SmallestPoint[1])
                    )
                    newAngle = self.dotProductAngle(self.vNewLeg1, self.vNewLeg2)

                     # If the legs cross over during adjustment
                    if (abs(newAngle * self.smallestAngle) > 0) and (newAngle * self.smallestAngle < 0):
                        self.leg1SmallestPoint = self.leg2SmallestPoint
                        self.smallestAngle = 0
                        return False

                    if abs(newAngle) < abs(self.smallestAngle):
                        self.smallestAngle = newAngle
                        self.leg1SmallestPoint = pair
            
        # If Leg 2 (S1 midpoint to S2 end) intersects the outside region
        if (    intersects(polygon, buffer(self.leg2, self.epsilon * 1e-1)) \
                and (contains(geometry, self.leg2) == False)
        ):
            if polygon.geom_type == 'LineString':
                coordinates2 = list(polygon.coords)
            else:
                coordinates2 = list(polygon.exterior.coords)

            leg2Start = buffer(self.midpoint, self.epsilon)
            leg2End = buffer(neighbor.end, self.epsilon)

            if (    (contains(leg2Start, intersection(polygon, buffer(self.leg2, self.epsilon * 1e-2))) == False) \
                    and (contains(leg2End, intersection(polygon, buffer(self.leg2, self.epsilon * 1e-2))) == False)
            ):
                for pair in coordinates2:
                    if Point(pair) == self.midpoint or contains(buffer(self.midpoint, self.epsilon * 1e2), Point(pair)):
                        continue
                    self.vNewLeg1 = self.vectorHelper(
                            (self.midpoint.x, self.midpoint.y), 
                            (self.leg1SmallestPoint[0], self.leg1SmallestPoint[1])
                    )
                    self.vNewLeg2 = self.vectorHelper(
                            (self.midpoint.x, self.midpoint.y), 
                            pair
                    )
                    newAngle = self.dotProductAngle(self.vNewLeg1, self.vNewLeg2)

                    # If legs cross over each other during adjustment
                    if (abs(newAngle * self.smallestAngle) > 0) and (newAngle * self.smallestAngle < 0): 
                        self.leg2SmallestPoint = self.leg1SmallestPoint
                        self.smallestAngle = 0
                        return False

                    if abs(newAngle) < abs(self.smallestAngle):
                        self.smallestAngle = newAngle
                        self.leg2SmallestPoint = pair

        return True

    def drawOuterCircle(self):
        from shapely import Point, plotting, Polygon, MultiPoint, is_closed, get_coordinates, LineString
        import math
        
        '''Creates the distribution circle that will be compared to the analytic cosine (does not show plot).'''
        """Reference Variables/Important:
            self.outerCircle (Polygon)
        """

        # Creating the outer circle of S2 surfaces
        outerRadius = self.surfaceLength / 2 # outer circle has diameter equal to the length of self surface (S1)
        outerCenter = self.midpoint 
        outerCirclePoints = []
        numPoints = 500

        for i in range(numPoints):
                angle = (2 * math.pi) * (i / numPoints)
                x = outerCenter.x + outerRadius * math.cos(angle) # uses same center as the distribution circle
                y = outerCenter.y + outerRadius * math.sin(angle)
                outerCirclePoints.append((x, y))
    
        # Polygon object of the full outer circle, before splitting to correct side of normal
        self.fullCircle = Polygon(outerCirclePoints) 

        buffLine = self.segment.buffer(self.epsilon * 1e-7)
        splitCircles = self.fullCircle.difference(buffLine)
        if splitCircles.geoms[0].intersects(self.normal):
            self.outerCircle = splitCircles.geoms[0] # desired circle
        else:
            self.outerCircle = splitCircles.geoms[1] # desired circle

        return


    # # # # # # # # # # # 
    # PLOTTING FUNCTIONS #
    # # # # # # # # # # # 

    def plotSelf(self, ax=None, color='k', label=False, showCircle=True):
        from matplotlib.pyplot import subplots, ioff, Figure, Axes
        from shapely import plotting, buffer, Point

        '''Plots a surface.'''

        ioff()
        if ax is None:
            fig, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        elif not isinstance(ax, Axes):
            raise Exception('ax not a valid Figure or Axes object')
        ax.set_aspect('equal')

        plotting.plot_line(self.segment, ax, color=color, linewidth=2) # plots surface
        if label:
            ax.text(self.midpoint.x, self.midpoint.y, f"Surface {self.ID}", fontsize=8, color='black')

        if showCircle:
            plotting.plot_polygon(self.circle, ax, add_points=False, color='green', linewidth=1) # plots the distribution circle
            
        return ax

    def plotConnections(self,  ax=None, linewidth=2, 
            colorseed=1, **kwargs):
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
            raise Exception('ax not a valid Figure or Axes object')
        ax.set_aspect('equal')

        cmap=get_cmap('jet')
        cols = linspace(0,1,len(self.neighbors))
        random.Random(colorseed).shuffle(cols)
        colors = iter(cmap(cols))
        for neighid, neighbor in self.neighbors.items():
            plotting.plot_polygon(
                    neighbor['los'], 
                    add_points=False,
                    color=next(colors), 
                    linewidth=linewidth,
                    **kwargs
            )

    def showAnalyticPlot(self, ax, comparison=True, r_offset=1, showBothDist=False): 
        from shapely import Point, plotting, Polygon, MultiPoint, is_closed, get_coordinates, LineString
        from shapely.plotting import plot_points
        from matplotlib.pyplot import subplots, ioff
        import matplotlib.pyplot as plt
        import math
        import numpy as np

        '''Plots the curve from the outer circle and compares it to the analytic equation plot.'''

        '''The 'comparison' variable determines if we are comparing to the analytic cosine distribution.
            Only set to true if modeling a COSINE distribution (offset = 1) 
            and you want to compare it to the analytic.'''

        '''Set showBothDist to True to display both uniform and cosine (geometric and analytic) on the same plot.
            to do this, set r_offset to 1. Also be sure to call analyticUniform on self.'''

        ioff()
        fig = subplots()
        ax.set_aspect('auto')
        plt.xlabel("Angle (Radians)")
        plt.ylabel("Fractional Area")

        self.distributionCircle(r_offset)
        self.drawOuterCircle()
        
        # CREATING THE PLOT OF ANGLE VS AREA
        s2Points = get_coordinates(self.outerCircle) # points defining the outer circle of S2 surfaces
        s2Start = s2Points[0] # The starting point of the first S2 surface
        plotPoints = [] # Points to plot for the outer circle plot
        pdfArea = 0 # "C" for the outer circle plot

        for i in range(1, len(s2Points), 1):

            # S2 Surface
            s2End = s2Points[i]
            s2Surface = Surface((s2Start[0], s2Start[1]), (s2End[0], s2End[1]), i)

            # Removes some outlier points 
            if (s2Surface.segment.length >= (self.surfaceLength - 1)) and (s2Surface.segment.length <= (self.surfaceLength + 1)): 
                s2Start = s2End
                continue
            
            # Edge case outliers (mainly for uniform distribution)
            startTuple = (s2Start[0], s2Start[1])
            endTuple = (s2End[0], s2End[1])
            if startTuple not in list(self.fullCircle.exterior.coords) or endTuple not in list(self.fullCircle.exterior.coords):
                s2Start = s2End
                continue

            # Vector representations of the normal and the line from the midpoint of S2 to the midpoint of S1
            vNormal = self.vectorHelper((self.normalStart.x, self.normalStart.y), (self.normalEnd.x, self.normalEnd.y))
            vS2 = self.vectorHelper((self.midpoint.x, self.midpoint.y), (s2Surface.midpoint.x, s2Surface.midpoint.y))

            # Call helper function that uses dot product to calculate angle between vectors
            angle = self.dotProductAngle(vNormal, vS2) # Plot on x-axis

            areaValue, _ = self.intersectionArea(s2Surface) # Plot on y-axis

            dTheta = self.dotProductAngle(self.vLeg1, self.vLeg2)
            pdfArea += areaValue * dTheta

            # Add the point and continue to the next iteration of the loop
            plotPoints.append(Point(angle, areaValue))
            s2Start = s2End

        # Just plotting the outer circle plot points-- no need to normalize
        if comparison == False:
            for point in plotPoints:
                plot_points(point, ax, color='black')
             
        # Plotting the analytic cosine distribution from the equation y = (1/2pi) * (1 + cos(x)), scaled to be from -pi/2 to pi/2
        # Normalizing the outer circle area value points as well using pdfArea
        else: # if comparison == True
            cosPoints = []
            adjustedPlotPoints = []
            for plotPoint in plotPoints:
                # Red points (analytic)
                cosXval = plotPoint.x
                cosYval = (1/math.pi)*(1 + math.cos(cosXval*2)) 

                cosPoints.append(Point(cosXval, cosYval))

                # Normalized outer circle points
                adjustedArea = plotPoint.y / abs(pdfArea)
                adjustedPlotPoints.append(Point(plotPoint.x, adjustedArea))
            
            for adjPoint in adjustedPlotPoints:
                plot_points(adjPoint, ax, color='black')

            for cosPoint in cosPoints:
                plot_points(cosPoint, ax, color='red', marker='1')
        

        '''Print statements for sanity check of total fractional area and the PDF Area.'''
        # total = 0
        # for point in plotPoints:
        #     total += point.y 

        # print("Total Fractional Area: ", total)
        # print("PDF Area: ", pdfArea)

        # Generate the plot
        if self.r_offset == 0:
            yUpperLim = plotPoints[-1].y + plotPoints[-1].y * 0.1
            ax.set_ylim(bottom=0, top=yUpperLim)
        plt.show(block=False)

        return

    def analyticUniform(self, ax):
        from shapely import Point, plotting, Polygon, MultiPoint, is_closed, get_coordinates, LineString
        from shapely.plotting import plot_points
        from matplotlib.pyplot import subplots, ioff
        import matplotlib.pyplot as plt
        import math
        import numpy as np

        ioff()
        fig = subplots()

        self.distributionCircle(0) # uniform distribution
        self.drawOuterCircle()
        
        # Creating the plot of angle vs. flux
        s2Points = get_coordinates(self.outerCircle) # points defining the outer circle of S2 surfaces
        s2Start = s2Points[0] # The starting point of the first S2 surface
        plotPoints = [] # Points to plot for the outer circle plot
        pdfArea = 0 # "C" for the outer circle plot

        for i in range(1, len(s2Points), 1):
            # S2 Surface
            s2End = s2Points[i]
            s2Surface = Surface((s2Start[0], s2Start[1]), (s2End[0], s2End[1]), i)

            # Removes some outlier points
            if (s2Surface.segment.length >= (self.surfaceLength - 1)) and (s2Surface.segment.length <= (self.surfaceLength + 1)): 
                s2Start = s2End
                continue
            
            # Edge case outliers
            startTuple = (s2Start[0], s2Start[1])
            endTuple = (s2End[0], s2End[1])
            if startTuple not in list(self.fullCircle.exterior.coords) or endTuple not in list(self.fullCircle.exterior.coords):
                s2Start = s2End
                continue

            #Vector representations of the normal and the line from the midpoint of S2 to the midpoint of S1
            vNormal = self.vectorHelper((self.normalStart.x, self.normalStart.y), (self.normalEnd.x, self.normalEnd.y))
            vS2 = self.vectorHelper((self.midpoint.x, self.midpoint.y), (s2Surface.midpoint.x, s2Surface.midpoint.y))

            # Call helper function that uses dot product to calculate angle between vectors
            angle = self.dotProductAngle(vNormal, vS2) # Plot on x-axis
            areaValue, _ = self.intersectionArea(s2Surface) # Plot on y-axis

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
            plot_points(adjPoint, ax, color='green')
        

        '''Print statements for sanity check of total fractional area and the PDF Area.''' # # #'''
        # total = 0
        # for point in plotPoints:
        #     total += point.y 

        # print("Total Fractional Area: ", total)
        # print("PDF Area: ", pdfArea)

        # Generate the plot
        plt.show(block=False)

        return

    def showTwoSurfacePlot(self, s2, r_offset=0):
        from shapely import Point, LineString, plotting
        from matplotlib.pyplot import subplots, ioff
        import matplotlib.pyplot as plt

        ioff()
        fig, ax = subplots()
        ax.set_aspect('equal')

        self.distributionCircle(r_offset)
        self.intersectionArea(s2)

        plotting.plot_line(self.segment, ax, color='red', linewidth=2) # plots surface
        ax.text(self.start.x, self.start.y, "Surface 1", color='red')

        plotting.plot_line(self.normal, ax, color='blue', linewidth=2) # plots normal line to surface
        ax.text(self.normalEnd.x, self.normalEnd.y, "Normal (S1)", color='blue')

        plotting.plot_line(self.dCircleCenter, ax, color='green') # plots point at center of distribution circle
        plotting.plot_polygon(self.circle, ax, add_points=False, color='green', linewidth=1) # plots the distribution circle
        ax.text(self.dCircleCenter.x, self.dCircleCenter.y, self.distType, color='green')

        plotting.plot_polygon(self.triangle, ax, color='orange', linewidth=2) # plots the triangle connecting to another surface

        plotting.plot_line(s2.segment, ax, color='black', linewidth=2) # plot surface2
        ax.text(s2.midpoint.x, s2.midpoint.y, "Surface 2", color='orange')

        plotting.plot_polygon(self.overlapShape, ax, add_points=False, color='black', linewidth=2) # displays the overlapping area
        ax.text(self.overlapShape.centroid.x, self.overlapShape.centroid.y, "Overlap Area", color='black')

        plt.show(block=False)

        return

    def showOuterCirclePlot(self, r_offset=1):
        from shapely import Point, LineString, plotting
        from matplotlib.pyplot import subplots, ioff
        import matplotlib.pyplot as plt

        ioff()
        fig, ax = subplots()
        ax.set_aspect('equal')

        self.distributionCircle(r_offset)
        self.drawOuterCircle()

        plotting.plot_line(self.segment, ax, color='red', linewidth=2) # plots surface
        ax.text(self.start.x, self.start.y, "Surface 1", color='red')

        plotting.plot_line(self.normal, ax, color='blue', linewidth=2) # plots normal line to surface
        ax.text(self.normalEnd.x, self.normalEnd.y, "Normal (S1)", color='blue')

        plotting.plot_line(self.dCircleCenter, ax, color='green') # plots point at center of distribution circle
        plotting.plot_polygon(self.circle, ax, add_points=False, color='green', linewidth=1) # plots the distribution circle
        ax.text(self.dCircleCenter.x, self.dCircleCenter.y, self.distType, color='green')

        plotting.plot_polygon(self.outerCircle, ax, add_points=True, color='gray') # displays the "outerCircle"

        plt.show(block=False)

        return



    # # # # # # # # # # # #
    # # HELPER FUNCTIONS # # 
    # # # # # # # # # # # #
    def normalHelper(self, dx, dy, endX, endY, init): 

        '''Helper function to find the normal direction.'''

        if self.end.y > self.start.y: # +x
            endX += dy
            if self.end.x > self.start.x:
                endY -= dx # (+x, -y)
            else:
                endY += dx # (+x, +y)
        else: # -x
            endX -= dy
            if self.end.x < self.start.x:
                endY += dx # (-x, +y)
            else:
                endY -= dx # (-x, -y)
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

        crossProduct = v1[0]*v2[1] - v1[1]*v2[0]
        if crossProduct < 0:
            angle = -angle

        return angle
 


        







