
class DEGAS2runner:
    """ Object for coupling UEDGE to DEGAS2

    """
    def __init__(self, case, inpath, runpath, material='C', recyc_coef=1.,
        overwrite=False, vesselfile='vessel.dat', uefile='uedata.u', **kwargs):
        from os import getenv, makedirs
        from os.path import isdir

        self.degas2path = getenv("DEGAS2_PATH")
        if self.degas2path is None:
            raise OSError("$DEGAS2_PATH not set. Aborting")
        self.inpath = inpath
        if not isdir(self.inpath):
            raise OSError(f"DEGAS2 input file path {self.inpath} does not exist")
        self.runpath = runpath
        makedirs(self.runpath, exist_ok=overwrite)
        self.runfilepath = f"{self.runpath}/degas2infiles"
        makedirs(self.runfilepath, exist_ok=overwrite)
        self.outpath = f"{self.runpath}/degas2output"
        makedirs(self.outpath, exist_ok=overwrite)

        # Store filenames
        self.vesselfile = vesselfile
        self.uefile = uefile

        # TODO: add flexibility to materials and recycling coeffs, etc
        self.material = material
        self.recyc_coef = recyc_coef
        
        self.ionspecies, self.gasspecies = case.about.species_setup(True)

        """ Links class to uetools.Case functions """
        self.get = case.get
        self.plot = case.plot

    def setup_degas2_run(self, vesselfile=None, **kwargs):
        """ Creates DEGAS2 input files and intializes directories """
        from shutil import copytree
        if vesselfile is None:
            vesselfile = self.vesselfile
        else:
            self.vesselfile = vesselfile
        copytree(self.inpath, self.runfilepath, dirs_exist_ok=True)
        self.write_uedge_data()
        self.write_degas2in()
        self.write_rb()
        self.define_boundaries(**kwargs)
        self.setup_dg(plot=False, **kwargs)
        self.write_vessel_file(vesselfile)
        self.write_dgin("dg.in", **kwargs)


    def run_degas2(self):
        """ Exectures DEGAS2 run commands """
        command = {
            "datasetup": "",
            "problemsetup": f"{self.runfilepath}/pr.in",
            "definegeometry2d": f"{self.runfilepath}/dg.in",
            "defineback": f"{self.runfilepath}/db.in",
            "tallysetup": f"{self.runfilepath}/tally.input",
            "flighttest": ""
        }
        for cmd, arg in command.items():
            try:
                subprocess.run(f"{cmd} {arg}", shell=True, check=True,
                                        capture_output=False, text=True)
            except subprocess.CalledProcessError as e:
                print(f"Command '{cmd} {arg}' failed with return code {e.returncode}")
                print(e.stderr)
    
        # TODO: implement degas_runner here
        return

    def write_degas2in(self, problemname='pr', tallyname='tally',
        geometryout='dg', backgroundout='bk', outputfile='degas2_output'):
        print(f"Writing degas2.in to {self.runfilepath}")
        # TODO: harden with lookups? Or standardized locations
        inoutfilepaths = {
            'elements': 'data',
            'species': 'data',
            'materials': 'data',
            'reaction': 'data',
            'pmi': 'data',
        }
        degasin = {
            'problem': problemname,
            'tally': tallyname
        }
        degasout = {
            'geometry': geometryout,
            'background': backgroundout,
            'output': outputfile
        }
        inout = {
            '_infile': ['.input', f"{self.runfilepath}"],
            'file': ['.nc', f"{self.outpath}"]
        }
        with open(f"{self.runfilepath}/degas2.in", 'w') as f:
            for key, path in inoutfilepaths.items():
                for keyapp, fileapp in inout.items():
                    f.write(f"{key}{keyapp} {self.degas2path}/{path}/{key}{fileapp[0]}\n")
                f.write('\n')
            for key, file in degasin.items():
                for keyapp, fileapp in inout.items():
                    f.write(f"{key}{keyapp} {fileapp[1]}/{file}{fileapp[0]}\n")
                f.write('\n')
            for key, file in degasout.items():
                f.write(f"{key}{keyapp} {self.outpath}/{file}.nc\n")
                f.write('\n')
        print(f"    Successfully wrote degas2.in to {self.runfilepath}")
                
    def write_rb(self):
        print(f"Writing db.in to {self.runfilepath}")
        with open(f"{self.runfilepath}/db.in", 'w') as f:
            f.write(f"plasma_file {self.runfilepath}/rb.in\n")
        print(f"    Successfully wrote db.in to {self.runfilepath}")

        ion_species = ''
        for species in self.ionspecies:
            if 'D0' not in species.upper():
                ion_species = f"{ion_species} {species.strip().replace('1','')}"
        print(f"Writing rb.in to {self.runfilepath}")
        with open(f"{self.runfilepath}/rb.in", 'w') as f:
            f.write(f"uedge_file {self.runfilepath}/{self.uefile}\n")
            f.write(f"ion_species {ion_species}\n")
            f.write(f"polygon_file {self.outpath}/polygons.nc")
        print(f"    Successfully wrote rb.in to {self.runfilepath}")

    def write_uedge_data(self, uefile=None):
        from Forthon import packageobject
        if uefile is None:
            uefile = self.uefile
        else:
            self.uefile = uefile
        runid = self.get('runid')
        if runid is None:
            runid = ""
        print(f"Writing {self.uefile} to {self.runfilepath}")
        packageobject('bbb').__getattribute__('writemcnfile')(
            f"{self.runfilepath}/{uefile}", 
            runid
        )
        print(f"    Successfully wrote {self.uefile} to {self.runfilepath}")

    def define_boundaries(self, bounds=None,  **kwargs):
        print("Identifying geometry and setting up zones")
        rm, zm = self.get('rm'), self.get('zm')
        # Define bounding box
        if bounds is None:
            self.bounds = [
                rm.min()*0.95,
                rm.max()*1.05,
                zm.min()*0.95,
                zm.max()*1.05
            ]
        else: 
            self.bounds = bounds


        if self.get('geometry')[0].decode('UTF-8').strip().lower() == 'snull':

            self.boundaries = {
                'SOL': {
                    'target': (
                                    (rm[-1, -1, 1], zm[-1, -1, 1]),
                                    (rm[0, -1, 2], zm[0, -1, 2]) 
                    ),
                    'segment': [None, None],
                    'isegment': [None, None],
                    'intersects': [], 
                    'CW': [False, True], 
                },
                'PFR': { 
                    'target': ( 
                                    (rm[0, 0, 4], zm[0, 0, 4]), 
                                    (rm[-1, 0, 3], zm[-1, 0, 3])
                    ),
                    'segment': [None, None],
                    'isegment': [None, None],
                    'intersects': [], 
                    'CW': [False, True], 
                }
            }

            ####################################### 
            #                                     # 
            #                 ####                # 
            #                #    #               # 
            #               #      #              # 
            #              #   ##   #             # 
            #             #   #  #   #            # 
            #           #    #     #  #           # 
            #          #    #       #  #          #
            #         #   #          #   #        # 
            #        #   #     ##     #   #       # 
            #       #   #     #  #     #   #      # 
            #      #   #     #    #     #   #     # 
            #    #    #     #   1   #    #   #    # 
            #   #    #      # CORE  #      #  #   # 
            #   #    #       #     #      #   #   # 
            #   #     #       #  #       #    #   # 
            #   #   4  #       ##      #      #   #
            #   #  SOL  #             #       #   # 
            #    #       #           #       #    # 
            #      #      #         #       #     #  
            #       #    #     2     #     #      #    
            #   5    #  #   PLASMA    #  #        # 
            # WALLS   ##               ##         #    
            #         ##               ##         # 
            #        #  #      ##     #   #       # 
            #       #    #    #  #   #     #      # 
            #     #       #  # 3  # #       #     # 
            #    #          # PFR #          #    # 
            #   #     6      #   #            #   # 
            #  #   TARGETS    ###              #  # 
            # #                                 # # 
            #######################################
            self.zones = {
                'CORE': {
                    'type': 'exit',
                    'boundaries': {
                        '1-edge': [self.get('ixpt1')[0]+1, self.get('ixpt2')[0], 0, 0]
                    },
                    'triangulation': 'triangulate_to_zones'
                },
#                'PLASMA': {
#                },
                'PFR': {
                    'type': 'plasma',
                    'boundaries': {
                        '1-edge': [0, self.get('ixpt1')[0], 0, 0],
                        '2-edge-xcut': [self.get('ixpt2')[0]+1, self.get('ixpt2')[0]+1, 0, 0],
                        '3-edge': [self.get('ixpt2')[0]+1, self.get('nx'), 0, 0],
                        '4-wall': {
                            'intersects': self.boundaries['PFR'],
                            'connection': True,
                            'start': False,
                            'end': False
                        }
                    },
                    'triangulation': 'triangulate_to_zones'
                },
                'SOL': {
                    'type': 'plasma',
                    'boundaries': {
                        '1-edge-reverse': [self.get('ny'), self.get('ny')],
                        '2-wall': {
                            'intersects': self.boundaries['SOL'],
                            'connection': True,
                            'start': False,
                            'end': False
                        }
                    },
                    'triangulation': 'triangulate_to_zones'
                },
                'WALLS': {
                    'type': 'solid',
                    'boundaries': {
                        '1-outer': [0, 1, 2, 3],
                        '2-wall-reverse': {
                            'intersects': self.boundaries['SOL'],
                            'connection': True,
                            'start': False,
                            'end': False
                        }
                    },
                    'material': self.material,
                    'recyc_coef': self.recyc_coef,
                    'triangulation': 'triangulate_polygon'
                },
                'TARGETS': {
                    'type': 'solid',
                    'boundaries': {
                        '1-outer': [3, 0],
                        '2-wall': {
                            'intersects': self.boundaries['SOL'],
                            'connection': False,
                            'start': True,
                            'end': False
                        },
                        '3-edge-reverse': [0, 0, 0, self.get('ny')],
                        '4-wall-reverse': {
                            'intersects': self.boundaries['PFR'],
                            'connection': True,
                            'start': False,
                            'end': False
                        },
                        '5-edge': [self.get('nx'), self.get('nx'), 0, self.get('ny')],
                        '6-wall': {
                            'intersects': self.boundaries['SOL'],
                            'connection': False,
                            'start': False,
                            'end': True
                        },
                    },
                    'material': self.material,
                    'recyc_coef': self.recyc_coef,
                    'triangulation': 'triangulate_polygon'
                }
            }

        else:
            raise Exception('Only "snull" geometries zoning implemented. Aborting!')

        print("    Zone setup completed")


    def get_limiter(self, limiter=None, maxlength=None, **kwargs):
        from shapely import LinearRing, Point, LineString, Polygon 
        from numpy import array, cross
        print("Getting limiter nodes...")
        # Get Shapely object for vessel
        if limiter is None:
            self.limiter = LinearRing( zip( self.get("xlim"), self.get("ylim")) )
        else: 
            self.limiter = LinearRing( limiter )
        # Identify self-intersecting/folding points
        self.limiter_array = array(self.limiter.coords)
        inter = []
        # Iterate through all points
        for i in range(1,len(self.limiter_array)-1):
            # Check if the following point lies on the line segment
            # made up by all previous point s
            line = LineString(self.limiter_array[:i+1])
            if line.distance(Point(self.limiter_array[i+1])) < 1e-5:
                # If so, the Point is folding/self-intersecting: store index
                inter.append(i+1)
        # Remove folding/self-intersecting points
        if len(inter)>0:
            print(" ...Removing folding nodes")
            # Cast as list, pop points, and cast as array
            self.limiter_array = list(self.limiter_array)
            # Reverse order to avoid index-issues
            for i in inter[::-1]:
                self.limiter_array.pop(i)
            self.limiter_array = array(self.limiter_array)
        # Refine resolution, if requested
        if maxlength is not None:
            self.limiter_array = array(LinearRing(self.limiter_array).segmentize(maxlength).coords)
            print(" ...Refining resolution")
        # Ensure orientation of limiter is CW
        COG = array(Polygon(self.limiter_array).centroid.coords)[0]
        if cross( self.limiter_array[0] - COG, self.limiter_array[1] - COG)>0:
            print(" ...Orienting vessel clock-wise")
            self.limiter_array = self.limiter_array[::-1]
        print("    Vessel node setup complete")
        return


    def setup_dg(self, limiter=None, maxlength=None, plot=False, roll_array=True, **kwargs):
        """ Writes DEGAS2-compatible vessel geometry """
        from matplotlib.pyplot import subplots
        from shapely.ops import nearest_points, orient
        from shapely import Point
        from numpy import array, roll 

        def is_on_segment(p0, p1, p2, epsilon=1e-6):
            return (sum( (p0-p1)**2) + sum((p0-p2)**2) - sum((p1-p2)**2))<epsilon

        # Get limiter array
        self.get_limiter(limiter, maxlength, **kwargs)
        # Get UEDGE grid end-points
        self.define_boundaries(**kwargs)
        # Identify vessel segment containing UE grid corner
        print("Setting up vessel intersects...")
        for zonename, zone in self.boundaries.items():
            for point in range(2):
                # Find closest vessl point to target corner
                pt = nearest_points(self.limiter, Point(zone['target'][point]))[0]
                pt = array((pt.x, pt.y))
                # Loop trough all segments
                for i in range(len(self.limiter_array)-1):
                    p1, p2 = self.limiter_array[i], self.limiter_array[i+1]
                    if is_on_segment(pt, p1, p2):
                        zone['segment'][point] = (p1, p2)
                        zone['isegment'][point] = (i, i+1)
                        break
                if zone['segment'][point] is None:
                    p1, p2 = self.limiter_array[0], self.limiter_array[-1]
                    if is_on_segment(pt, p1, p2):
                        zone['segment'][point] = (p1, p2)
                        zone['isegment'][point] = (0, len(self.limiter_array)-1)
                if zone['segment'][point] is None:
                    raise Exception("Reference point not on vessel")
                # Identify the correct point 
                zone['intersects'].append(zone['isegment'][point][zone['CW'][point]])
        # Roll vessel so first intersect point has index 0
        if roll_array:
            print(" ...Reordering vessel nodes")
            ref = self.boundaries['SOL']['intersects'][1]
            self.limiter_array = roll(self.limiter_array, -ref, axis=0)
            for key, zone in self.boundaries.items():
                for i in range(2):
                    zone['intersects'][i] = zone['intersects'][i] - ref
                    if zone['intersects'][i] < 0:
                        zone['intersects'][i] = len(self.limiter_array) + zone['intersects'][i]
        # Plot geometry
        if plot:
            f, ax = subplots()
            ax.plot(*self.limiter_array[0], "or")
            for key, zone in self.boundaries.items():
                for pt in zone['target']:
                    ax.plot(*pt, 'v'*(key=="PFR")+'d'*(key=="SOL"), color='b')
                for pt in zone['intersects']:
                    ax.plot(*self.limiter_array[pt], 'v'*(key=="PFR")+'d'*(key=="SOL"), color='r')
            ax.plot(*self.limiter_array.T, ".-k")
        print("    Vessel intersect setup completed")
        return

                
    def write_vessel_file(self, outfile, **kwargs):
        print(f"Writing {outfile} to {self.runfilepath}")
        with open(f"{self.runfilepath}/{outfile}", 'w') as f:
            f.write(f'1\n{len(self.limiter_array)}\n')
            for p in self.limiter_array:
                f.write(f'{p[0]:.8f} {p[1]:.8f}\n')
        print(f"    {outfile} written successfully.")
        return 


    def write_dgin(self, outfile, tab=4*' ', **kwargs):
        i = 1
        print(f"Writing {outfile} to {self.runfilepath}")
        with open(f"{self.runfilepath}/{outfile}", 'w') as f:
            # Write initialization block
            f.write("symmetry cylindrical\n")
            f.write(f"uedge_mesh {self.runfilepath}/{self.uefile}\n")
            f.write(f"wallfile {self.runfilepath}/{self.vesselfile}\n")
            f.write(f"bounds {self.bounds[0]:.3f} {self.bounds[1]:.3f} {self.bounds[2]:.3f} {self.bounds[3]:.3f}\n")
            f.write("end_prep\n")
            # Start writing zones
            for zone, data in self.zones.items():
                # Required entries
                f.write(f'\n# {zone.upper()}\n')
                f.write(f'new_zone {data["type"]}\n')
                f.write(f'new_polygon\n{tab}stratum {i}\n')
                # Wall material data
                if data['type'].strip().lower() == 'solid':
                    f.write(f'{tab}material {data["material"]}\n') 
                    f.write(f'{tab}recyc_coef {data["recyc_coef"]:.2f}\n') 
                # Write boundaries of zones to file
                for btype, boundary in data['boundaries'].items():
                    pre = btype.split('-')[1]
                    try:
                        app = btype.split('-')[2]
                    except:
                        app = ''
                    indices = ''
                    # Check if the bounday is defined on UEDGE/outer nodes
                    if not isinstance(boundary, dict):
                        for index in boundary:
                            indices = f'{indices} {index}'
                        f.write(f'{tab}{pre} {indices.strip()} {app}\n')
                    # If not, access vessel indices from structs
                    else:
                        wallrange = ""
                        intersects = boundary['intersects']['intersects']
                        intersects.sort()
                        if boundary['connection']:
                            for i in intersects:
                                wallrange = f"{wallrange} {i}"
                        else:
                            i = min(intersects)*boundary['start'] + max(intersects)*boundary['end']
                            for j in range(2):
                                wallrange = f"{wallrange} {i}"
                        f.write(f'{tab}{pre} 1 {wallrange.strip()} {app}\n')
                        
                f.write(f'{tab}triangulate {data["triangulation"]}\n')

                i += 1
        
            f.write(f'\npolygon_nc_file {self.runfilepath}/polygons.nc\n')
        print(f"{tab}{outfile} written successfully.")


        return
        






















class VacuumTransport_coupling:
    def __init__(self, case):
        """ Links class to uetools.Case functions """
        self.get = case.get
        self.plot = case.plot

    def define_vacuum_region(self, **kwargs):
        """ Initializer for vacuum region plotting 
        
        Calls define_snull_vacuum_region or
        define_dnull_vacuum_region (not implemented) depending
        on the case geometry.
        """
        geo = self.get('geometry')[0].decode('UTF-8').strip() 
        if geo in ('snull', 'uppersn'):
            self.define_snull_vacuum_region(**kwargs)
        else:
            raise Exception("Vacuum region definition not implemented" +
                f"geometry {geo}.")
            
    def get_vacuum(self, boundary, vessel, north):
        """ Function that returns the clockwise vacuum region boundary
        
        Start is defined as the first point of the vessel region
        
        Parameters
        ----------
        boundary - Nx2 array of plasma boundary points
        vessel - Mx2 array of points defining the vessel
        north - boolean defning main vacuum region (True)
                or PFR vacuum region (False)

        Returns
        -------
        vacuum  - Kx2 array of vacuum region boundaries
        vesselind - int defning last point on vessel in vacuum
        intersects  - Tuple of (R,Z) coords of first and second
                        intersects between boundary and vessel
        """
        from numpy import cross, roll, array
        retval = []
        # For northern boundary, calculate polarity from plasma crown
        if north:
            ref = boundary[int(len(boundary)/2)]
        else:
            # For south boundary, reverse order of boundary points
            # to maintain polarity
            boundary = boundary[::-1]
            # Set reference point to SW of the vessel to ensure
            # polarity accuracy regardelss of pf shape
            ref = (vessel[:,0].min()-0.1, vessel[:,1].min()-0.1)
        # Get the two endpoints
        endpoints = (boundary[0], boundary[-1])
        for j in range(2):
            # Calculate the distance from the endpoint
            mindist =   (vessel[:,0] - endpoints[j][0])**2 \
                        + (vessel[:,1] - endpoints[j][1])**2
            # Find th indices of the two closest points
            minind = sorted(range(len(mindist)), key=lambda sub: mindist[sub])[:2]
            # Find the RH polar point (clockwise) closest intersect with start point
            i = 0**((-1)**(not north)*cross(endpoints[j]-ref, vessel[minind[j]]-ref)<0)
            retval.append(vessel[minind[i]])
            if j == 0:
                # Roll the array to align vessel start to start point
                vessel = roll(vessel, -minind[i], axis=0)
        vacuum = array(tuple(vessel[:minind[i]+1]) + tuple(boundary)[::-1])
        return vacuum, minind[i], retval

    def run_triangle(self, fname, regiondict, triaoptions = "-pqDne", maxarea=None, **kwargs):
        """ Creates Triangle input and runs Triangle under triangle_out

        Triangle source code available from:
        http://www.cs.cmu.edu/~quake/triangle.html

        Variables
        ---------
        fname - str of file name for files written
        regiondict - dictionary containing the data for the triangulation
                Consists of subdicts:
                filled - Nx2 array of points bounding region to be filled
                        with triangles
                holes - Mx2 array of regions to be treated as holes
        triaoptions - optional (default: -pqDne)
                Optional arguments passed to Triangle
        maxarea - optional (default: None)
                Optional restriction on max triangle area for refinement.
                Either list of floats for gradual refinement or float
                for single-step refinement

        Output
        ------
        Writes Triangle input and output files to triangle_out
        """
        from os.path import exists 
        from os import makedirs, environ
        from shutil import rmtree
        from shutil import which
        from subprocess import run

        # Ensure maxarea is compatible with Triangle calls
        if maxarea is not None:
            if isinstance(maxarea, (float, int)):
                maxarea = [maxarea]
            elif not isinstance(maxarea, list):
                raise TypeError("maxarea type '{}' not recognized".format(type(maxarea)))

        nnodes, nholes = 0, 0
        for region, points in regiondict['filled'].items():
            nnodes += len(points)
        for region, points in regiondict['holes'].items():
            nholes += len(points)

        if exists('triangle_out'):
            rmtree('triangle_out') 
        makedirs('triangle_out')
        # First, write data to a triangle file
        with open(f'triangle_out/{fname}.poly', 'w') as f:
            f.write("# Triangle input file from UEDGE \n")
            f.write(f"# NODES\n") 
            f.write("{} 2 0 1\n".format(nnodes))
            ii = 0
            for region, points in regiondict['filled'].items():
                f.write(f"# {region}\n") 
                for i in range(len(points)):
                    f.write("{} {:.5f} {:.5f}   {}\n".format(ii+1, *points[i], 1))
                    ii += 1
                f.write('\n')
            f.write('\n')
            f.write(f"# ELEMENTS\n") 
            f.write("{} 1\n".format(nnodes))
            ii = 0
            for region, points in regiondict['filled'].items():
                istart = ii
                f.write(f"# {region}\n")  
                for i in range(len(points)-1):
                    f.write("{} {} {} {}\n".format(ii+1, ii+1, ii+2, 1))
                    ii += 1
                f.write("{} {} {} {}\n".format(ii+1, ii+1, istart+1, 1))
                ii+=1
            f.write("{}\n".format(nholes))
            # TODO: Implement hole write routine

        if which('triangle') is None:
            triacmd = "{}/triangle".format(environ.get("TRIA_PATH")) 
        else:
            triacmd = "triangle"

        run([triacmd, triaoptions, f"triangle_out/{fname}.poly"])

        output = 1
        if maxarea is not None:
            for area in maxarea:
                run([
                    triacmd, 
                    triaoptions.replace('-', '-r')+"a{}".format(float(area)), 
                    f"triangle_out/{fname}.{output}"]
                )
                output += 1
        
    def read_triangle(self, fname, iteration=1):
        """ Reads grids produced by Triangle

        Variables
        ---------
        fname - full path to files to be read without extensions
        iteration - optional (default: 1)
            int of grid iteration to be read

        Returns
        -------
        griddata - dict containing the following subdicts
            nodes: dict of node indices and coordinates
            trias: dict of triangle indices and node coordinates
            elements: dict of triangle indices and node indices
            neighbors: dict of traingle indices and neighboring 
                        triangle indices
        """
        from numpy import array
        grid = {'nodes': {}, 'elements': {}, 'neighbors': {}, 'trias': {}}
        with open(f"{fname}.{iteration}.node") as f:
            nnodes = int(f.readline().split()[0])
            for i in range(nnodes):
                nodedata = f.readline().split()
                grid['nodes'][int(nodedata[0])] = (float(nodedata[1]), float(nodedata[2]))
        with open(f"{fname}.{iteration}.ele") as f:
            nele = int(f.readline().split()[0])
            for i in range(nele):
                eledata = f.readline().split()
                grid['trias'][int(eledata[0])] = array((
                    grid['nodes'][int(eledata[1])],
                    grid['nodes'][int(eledata[2])],
                    grid['nodes'][int(eledata[3])],
                    grid['nodes'][int(eledata[1])],
                ))
                grid['elements'][int(eledata[0])] = \
                    [int(x) for x in eledata[1:]]
        with open(f"{fname}.{iteration}.neigh") as f:
            nneigh = int(f.readline().split()[0])
            for i in range(nneigh):
                neighdata = f.readline().split()
                grid['neighbors'][int(neighdata[0])] = \
                    [int(x) for x in neighdata[1:]]

        return grid


    def plot_triangle(self, grid, ax=None):
        """ Plots the triangular grid
        
        Arguments
        ---------
        grid - griddata dict as created by read_triangle
        ax - optinal (default: None)
            axis to plot onto. If None, creates new figure
        """
        from matplotlib.pyplot import subplots
        if ax is None:
            f, ax = subplots()
        for idx, nodes in grid['trias'].items():
            ax.plot(nodes[:,0], nodes[:,1], 'k-', linewidth=0.3)
        

    def define_snull_vacuum_region(self, gridname="tria",   
        maxlength = 0.05, plot=True, limiter=None, **kwargs):
        """ Creates arrays enclosing the main-SOL and PF vacuum regions
        
        Calls Triangle if executable found in PATH or TRIA_PATH. Else,
        applies Scipy Delaunay routines.

        Arguments
        ---------
        gridname - optional (default: tria)
            name for Triangle-related files if Triangle available
        maxlength - optional (default: 0.05)
            maximum vessel segment length in meters. Necessary to 
            find approproiate vessel neighboring points for 
            automated vacuum region detection
        plot - optional (default: True)
            switch whether to plot results or not
        limiter - optional (default: None)
            Nx2 array of points ordered clockwise (snull) or 
            counterclockwise (usn) to contain the nodes
            specifying the vessel geometry.
        **kwargs - passed to run_triangle
    
        Returns
        -------
        main_vacuum - Nx2 array points enclosing main vacuum region
        pf_vacuum - Mx2 array of points enclosing PFR vacuum region
        """
        from shapely import LinearRing, LineString, plotting, Polygon
        from shapely.ops import linemerge
        from matplotlib.pyplot import subplots, triplot
        from numpy import array, roll, cross
        from copy import deepcopy
        from scipy.spatial import Delaunay
        from os import environ
        from shutil import which
                
        # Get shifts and switches for upper-single null geometries
        if (self.get("rmagx") + self.get("zmagx") == 0):
            disp = 2.8
        else:
            disp = self.get('zmagx')
        usn = self.get('geometry')[0].decode('UTF-8').strip() == 'uppersn'

        # Specify the vessel geometry: introduce copy for plotting purposes only
        if limiter is None:
            vessel_orig = LinearRing( zip(self.get('xlim'), self.get('ylim')))
        else:
            vessel_orig = LinearRing( limiter )
        if usn:
            limiter = array(list(zip(
                vessel_orig.coords.xy[0],
                vessel_orig.coords.xy[1]
            )))
            vessel_orig_plot = LinearRing( zip(
                limiter[:,0],
                (-1)**usn*limiter[:,1] + disp*usn
            ))
        else:
            vessel_orig_plot = vessel_orig
                
        # Get the nodes for the north and south plasma boundaries
        northpoints = array(list(zip(
                self.get('rm')[:-1,-1,2],
                self.get('zm')[:-1,-1,2]
        )))
        southpoints = array(tuple(zip(
                    self.get('rm')[:self.get('ixpt1')[0]+1,0,4],
                    self.get('zm')[:self.get('ixpt1')[0]+1,0,4]
        )) + tuple(zip(
                    self.get('rm')[self.get('ixpt2')[0]+1:,0,3],
                    self.get('zm')[self.get('ixpt2')[0]+1:,0,3]
        )))
        # Create a LineString for plotting
        north = LineString( northpoints ) 
        south = LineString( southpoints )
        # Segmentize the vessel to enforce max segments lengths
        vessel=vessel_orig.segmentize(maxlength)
        vessel_plot=vessel_orig_plot.segmentize(maxlength)
        # Get the segmentized vessel coordinates
        vesselpoints = array( list( zip(
            vessel.coords.xy[0],
            vessel.coords.xy[1]
        )))
        # Ensure right-handedness of final vessel points
        if usn:
            vesselpoints = vesselpoints[::-1]
        # Get vacuum region points and information
        main_points, nmain_points, nintersects = self.get_vacuum(northpoints, vesselpoints, True)
        pf_points, npf_points, sintersects = self.get_vacuum(southpoints, vesselpoints, False)
        # Account for USN flip
        main_points[:,1] = (-1)**usn*main_points[:,1] + usn*disp
        pf_points[:,1] = (-1)**usn*pf_points[:,1] + usn*disp
        # Set up triangulation dict
        triadata = {
            'filled': {
                'main_vacuum': main_points,
                'pfr_vacuum': pf_points,
            },
            'holes': {
            }
        }
        # Call triangulation routines
        if (which('triangle', path = environ.get("TRIA_PATH")) is None) \
            and (which('triangle') is None):
            print("Package Triangle not found. Add path to Triangle to\n"
                +"PATH or TRIA_PATH and re-run. Using Scipy Delaunay...")
            triangle = False
            mainTri = Delaunay(main_points)        
            pfTri = Delaunay(pf_points)        
        else:
            self.run_triangle(gridname, triadata, **kwargs)
            triangle = True
        # Read triangle output
        iteration = 1
        if ('maxarea' in kwargs):
            if kwargs['maxarea'] is not None:
                try:
                    iteration += len(kwargs['maxarea'])
                except:
                    iteration += 1



        if plot:
            # Initialize figure
            f, ax = subplots(1,5,figsize=(25,10))
            for a in ax:
                a.set_aspect('equal', adjustable='box') 
            # Plot UEDGE grid and vessel
            f = self.plot.mesh(
                rm=self.get('rm'),
                zm=self.get('zm')*(-1)**(usn), 
                ax=ax[0], lcfs=False
            )
            ax[0].set_title("Original grid")
            # Plot Original nodes and boundaries in first pane 
            plotting.plot_line(vessel_orig_plot, ax=ax[1], color='r')
            ax[1].plot(
                northpoints[:,0], 
                (-1)**usn*northpoints[:,1] + usn*disp,
                'ko-'
            )
            ax[1].plot(
                southpoints[:,0], 
                (-1)**usn*southpoints[:,1] + usn*disp,
                'ko-'
            )
            ax[1].set_title("Region outlines")
            # Plot refined vessel nodes and detected start/end points in    
            # second pane
            plotting.plot_line(vessel_plot, ax=ax[2], color='r')
            ax[2].plot(
                northpoints[:,0], 
                (-1)**usn*northpoints[:,1] + usn*disp,
                'ko-'
            )
            ax[2].plot(
                southpoints[:,0], 
                (-1)**usn*southpoints[:,1] + usn*disp,
                'ko-'
            )
            colors = ['g', 'b']
            for i in range(len(nintersects)):
                ax[2].plot(
                    nintersects[i][0], 
                    (-1)**usn*nintersects[i][1]+usn*disp, 
                    color=colors[i], marker='s', markersize=8
                )
            for i in range(len(sintersects)):
                ax[2].plot(
                    sintersects[i][0], 
                    (-1)**usn*sintersects[i][1]+usn*disp, 
                    color=colors[i], marker='d', markersize=8
                )
            ax[2].set_title(f"Refined vessel: {maxlength} cm") 
            # Plot detected vacuum regions in third pane
            plotting.plot_line( LinearRing(main_points), ax=ax[3], color='k')
            plotting.plot_line( LinearRing(pf_points), ax=ax[3], color='k')
            colors = ('g', 'b')
            ax[3].set_title("Vacuum regions")
            # Plot triangles 
            if triangle is False:
                ax[4].triplot(
                    main_points[:,0], 
                    main_points[:,1], 
                    mainTri.simplices, 
                    linewidth=0.5, color='k'
                )
                ax[4].triplot(
                    pf_points[:,0], 
                    pf_points[:,1], 
                    pfTri.simplices, 
                    linewidth=0.5, color='k'
                )
            else:
                # Read triangle input into dict
                triangledata = self.read_triangle(
                        f'triangle_out/{gridname}', iteration
                )
                # Plot triangles
                self.plot_triangle(triangledata, ax=ax[4])
            ax[4].set_title("Triangulation")
        return main_points, pf_points
       
