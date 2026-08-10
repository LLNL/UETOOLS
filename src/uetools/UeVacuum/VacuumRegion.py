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
    def __init__(self, case):
        self.coupling = case.coupling
        self.set = case.setue
        self.info = case.info
        self.tools = case.tools
        self.getue = case.getue

    def puffing_array_calc(self, region, point, current, puffing_matrix): # helper
        """Returns input puffing array based on location as an (R, Z) coordinate, input current, and the puffing transport matrix."""
        from shapely import Point, intersects
        import numpy
        print('enter puffing array calc')
        puff_point = Point(point)

        puffing_location = 0
        min_dist = float('inf')
        for i in range(region.numSurfaces):
            if i >= region.P:
                seg = region.surfaces[i]
                seg_x = seg.normalEndX
                seg_y = seg.normalEndY

                point_x = puff_point.x
                point_y = puff_point.y

                dist = numpy.sqrt((seg_x - point_x)**2 + (seg_y - point_y)**2)
                if dist < min_dist:
                    min_dist = dist
                    puffing_location = i

        region.matrices()
        puffing_array = region.getPuffingArray(puffing_matrix, puffing_location, current)
        print('error here 1')
        puffing_array = numpy.transpose(puffing_array)
        print('error here 2')
        return puffing_array, puffing_location
    
    def save_matrices(self, matrix_list, matrix_name_list, save_file_name, open_file=None): # helper
            """Saves matrices to hdf5 file (save file)."""
            import h5py

            def save_helper(f):
                subgroup = f.require_group('vnm/bbb')
                for name, matrix in zip(matrix_name_list, matrix_list):
                    if hasattr(matrix, 'toarray'):
                        matrix = matrix.toarray()
                    if name in subgroup:
                        if subgroup[name].shape == matrix.shape:
                            subgroup[name][...] = matrix
                        else:
                            del subgroup[name]
                            subgroup.create_dataset(name, data=matrix)
                    else:
                        subgroup.create_dataset(name, data=matrix)

            if open_file is not None:
                save_helper(open_file)
            else:
                with h5py.File(save_file_name, 'a') as f:
                    save_helper(f)

    def restore(self, save_file=None, sol_puff_dict=None):
        """Restores existing telematrices from the provided save file, and calculates puffing input arrays if desired.
        Uses the current save file if none is provided.
        - sol_puff_dict = {'point': __, 'current': __, 'region':__}
         'region' is optional-- pass in a pkl save of a VacuumRegion or a previously loaded VacuumRegion
        """
        import h5py
        import warnings
        import numpy
        from uedge import com, bbb

        if save_file == None:
            save_file = self.info['savefile']

        bbb.cftelematrix[:, :, :, :] = 0

        with h5py.File(save_file, 'a') as f:
            if 'vnm/bbb/cftelematrix' not in f: # Check for main telematrix
                warnings.warn('SOL telematrix not found in save file. Call self.vnm.generate()')
            else: # re-save main telematrix
                if len(f['vnm/bbb/cftelematrix'].shape) == 2:
                    dimension = len(f['vnm/bbb/cftelematrix'][:]) + 2
                    cftelematrix_full = numpy.zeros((2, dimension, dimension, 6))
                    cftelematrix_full[1, 1:-1, 1:-1, 0] = f['vnm/bbb/cftelematrix'][:]

                    bbb.cftelematrix[1, 1:-1 , 1:-1 , 0] = f['vnm/bbb/cftelematrix'][:]
                elif len(f['vnm/bbb/cftelematrix'].shape) == 2:
                    bbb.cftelematrix[1, 1:-1, 1:-1, 0] = f['vnm/bbb/cftelematrix'][1, 1:-1, 1:-1, 0]

                self.save_matrices([f['vnm/bbb/cftelematrix']], ['cftelematrix'], save_file, open_file=f)

            if 'vnm/bbb/cftelematrix_pf' not in f: # check for pfr telematrix
                warnings.warn('Private flux region telematrix not found in save file. Call self.vnm.generate.')
            else: # re-save pfr telematrix
                small_matrix = f['vnm/bbb/cftelematrix_pf'][:]

                if 'vnm/bbb/cftelematrix' not in f:
                    warnings.warn('Generate main SOL telematrix before generating PFR telematrix.')
                else:
                    dimension = len(f['vnm/bbb/cftelematrix'][:]) + 2
                    print(dimension)
                    print(len(f['vnm/bbb/cftelematrix'].shape))
                cftelematrix_full_pf = numpy.zeros((2, dimension, dimension, 6))
                cftelematrix_full_pf[0, 1:com.ixpt1[0]+1, 1:com.ixpt1[0]+1, 0] = small_matrix[:com.ixpt1[0], :com.ixpt1[0]] # 1
                cftelematrix_full_pf[0, 1:com.ixpt1[0]+1, com.ixpt2[0]+1:com.nx+1, 0] = small_matrix[:com.ixpt1[0], com.ixpt1[0]:] # 2
                cftelematrix_full_pf[0, com.ixpt2[0]+1:com.nx+1, 1:com.ixpt1[0]+1, 0] = small_matrix[com.ixpt1[0]:, :com.ixpt1[0]] # 3 
                cftelematrix_full_pf[0, com.ixpt2[0]+1:com.nx+1, com.ixpt2[0]+1:com.nx+1, 0] = small_matrix[com.ixpt1[0]:, com.ixpt1[0]:] # 4

                bbb.cftelematrix[0, 1:com.ixpt1[0]+1, 1:com.ixpt1[0]+1, 0] = small_matrix[:com.ixpt1[0], :com.ixpt1[0]] # 1
                bbb.cftelematrix[0, 1:com.ixpt1[0]+1, com.ixpt2[0]+1:com.nx+1, 0] = small_matrix[:com.ixpt1[0], com.ixpt1[0]:] # 2
                bbb.cftelematrix[0, com.ixpt2[0]+1:com.nx+1, 1:com.ixpt1[0]+1, 0] = small_matrix[com.ixpt1[0]:, :com.ixpt1[0]] # 3 
                bbb.cftelematrix[0, com.ixpt2[0]+1:com.nx+1, com.ixpt2[0]+1:com.nx+1, 0] = small_matrix[com.ixpt1[0]:, com.ixpt1[0]:] # 4

                self.save_matrices([small_matrix], ['cftelematrix_pf'], save_file, open_file=f)

            if 'vnm/bbb/puffing_matrix' not in f: # check for main puffing matrix
                warnings.warn('SOL puffing matrix not found in save file. Call self.vnm.generate.')
            else: # re-save main puffing matrix
                self.save_matrices([f['vnm/bbb/puffing_matrix']], ['puffing_matrix'], save_file, open_file=f)

            if 'vnm/bbb/puffing_array' not in f: # check for main puffing array
                warnings.warn('SOL puffing array not found in save file. Call self.vnm.restore or self.vnm.generate.')
            else: # re-save main puffing array
                self.save_matrices([f['vnm/bbb/puffing_array']], ['puffing_array'], save_file, open_file=f)
                self.getue('fngyo_use', cp=False)[1:-1, 0] = f['vnm/bbb/puffing_array'][: , 0]

            if sol_puff_dict != None and 'vnm/bbb/puffing_matrix' in f: # creating a new puffing array
                if 'region' not in sol_puff_dict:
                    (main, pf) = self.coupling.get_snull_vacuum_regions(maxlength=0.0087)
                    sol = VacuumRegion(main[0], P=main[1])
                else:
                    print('error here')
                    sol = VacuumRegion(sol_puff_dict['region'])

                print('made it here')
                puffing_array, puff_loc = self.puffing_array_calc(sol, sol_puff_dict['point'], sol_puff_dict['current'], 
                                                    f['vnm/bbb/puffing_matrix'][:])
                print('puffing_array:', puffing_array)
                print('save error')
                self.main_puffing_array = numpy.transpose(puffing_array)
                self.save_matrices([self.main_puffing_array], ['puffing_array'], save_file, open_file=f)
                self.getue('fngyo_use', cp=False)[1:-1, 0] = f['vnm/bbb/puffing_array'][: , 0]
                print('set error')
            elif sol_puff_dict != None and 'vnm/bbb/puffing_matrix' not in f: # error if can't make new puffing array
                warnings.warn('Cannot generate SOL puffing input. Call self.vnm.restore or self.vnm.generate.')

            if 'vnm/bbb/pumping_matrix_pf' not in f:
                warnings.warn('PFR pumping matrix not found. Call self.vnm.generate.')
            else:
                pump = f['vnm/bbb/pumping_matrix_pf'][:]
                self.save_matrices([pump], ['cftelematrix_pf'], save_file_name=save_file, open_file=f)
                bbb.cftelematrix[0, 1:com.ixpt1[0]+1, 1:com.ixpt1[0]+1, 0] = pump[:com.ixpt1[0], :com.ixpt1[0]] # 1
                bbb.cftelematrix[0, 1:com.ixpt1[0]+1, com.ixpt2[0]+1:com.nx+1, 0] = pump[:com.ixpt1[0], com.ixpt1[0]:] # 2
                bbb.cftelematrix[0, com.ixpt2[0]+1:com.nx+1, 1:com.ixpt1[0]+1, 0] = pump[com.ixpt1[0]:, :com.ixpt1[0]] # 3 
                bbb.cftelematrix[0, com.ixpt2[0]+1:com.nx+1, com.ixpt2[0]+1:com.nx+1, 0] = pump[com.ixpt1[0]:, com.ixpt1[0]:] # 4

            if 'vnm/bbb/pumping_matrix_sol' not in f:
                warnings.warn('SOL pumping matrix not found. Call self.vnm.generate.')
            else:
                pump = f['vnm/bbb/pumping_matrix_sol'][:]
                self.save_matrices([pump], ['cftelematrix'], save_file_name=save_file, open_file=f)

                if len(f['vnm/bbb/pumping_matrix_sol'].shape) == 2:
                    bbb.cftelematrix[1, 1:-1, 1:-1, 0] = pump
                elif len(f['vnm/bbb/pumping_matrix_sol'].shape) == 4:
                    bbb.cftelematrix[1, 1:-1, 1:-1, 0] = pump[1, 1:-1, 1:-1, 0]

        print('error with one of these')
        self.getue('isvacuummodel', cp=False)[0] = 1
        # bbb.isvacuummodel[0] = 1
        self.set('cfteleout', 1.0)
        print('finished execution')
    
    def generate(self, sol=None, pfr=None, save_file=None, sol_puff_dict=None, 
                 sol_pump_dict=None, pf_pump_dict=None, save=False, pump_plot=True, pfr_user=None,
                 overwrite=False):
        """Generates telematrices for given VacuumRegions.
        
        Keyword arguments:
        - sol -- VacuumRegion save of the main SOL; can be a pkl file name (default None).
         If none, a new VacuumRegion will be loaded. If a pkl file name is passed in, that pkl save VacuumRegion will be used.
        - pfr -- VacuumRegion save of the private flux region; can be a pkl file (default None).
         If none, a new VacuumRegion will be loaded. If a pkl file name is passed in, that pkl save VacuumRegion will be used.
        - save_file -- save file used to generate the Case. Uses the current save file if none are provided (default None).
        - sol_puff_dict -- {'point': __, 'current': __} (default None).
         'point' is a coordinate that identifies the location of the gas puff.
         'current' is a value that determines the strength of the gas puff.
        - sol_pump_dict -- {'box_coords':__, 'albedo':__} (default None).
         'box_coords' should be a list of four or more coordinates that create a box around the pumping portion.
           Note: even if the full surfaces is not within the box, it will be counted as a pumping surface. Additionally, if not enough points are provided for a box (3 or less),
           all wall surfaces in the region will be taken as pumping with the specified albedo.
         'albedo' is a value between 0 and 1 that determines the strength of the pumping
        - pf_pump_dict -- {'box_coords': __, 'albedo': __} (default None).
         Same logic as sol_pump_dict.
        - save -- if True, will save new pkl file with names 'SOL_Vacuum.pkl' and 'PFR_Vacuum.pkl' (default False).
        - pump_plot -- if True, plots the pumping surfaces (default True).
        - pfr_user -- set of user supplied points to define the private flux region (default None).
        - overwrite -- decides whether or not to write generated matrices into save file (default False).
        """

        import h5py
        import numpy
        from uedge import com, bbb

        (main, pf) = self.coupling.get_snull_vacuum_regions(maxlength=0.0087)
        self.pfrPoints = pfr_user ########
        self.box_loop = True
        pf_plasma_number = pf[1] ########
        if pfr_user is not None: ########
            pf = [pfr_user, pf_plasma_number] ########

        if save_file == None:
            save_file = self.info['savefile']
        
        def matrix_calculate(main_vac, pf_vac):
            """Generates the telematrices."""
            t_main = main_vac
            t_pf = pf_vac

            t_main.matrices()
            t_pf.matrices()

            cftelematrix_ = t_main.getOutputMatrix(t_main.AB_matrix, 1000000) # transport matrices
            cftelematrix_ = numpy.transpose(cftelematrix_)
            cftelematrix_pf = t_pf.getOutputMatrix(t_pf.AB_matrix, 1000000)
            cftelematrix_pf = numpy.transpose(cftelematrix_pf)

            puffing_matrix = t_main.getPuffingMatrix(1000000)

            self.cftelematrix = cftelematrix_
            self.cftelematrix_pf = cftelematrix_pf

            return [cftelematrix_, cftelematrix_pf, puffing_matrix]

        if sol == None or pfr == None: # generate new VacuumRegions
            # print('generate new Vacuum Regions')
            sol = VacuumRegion(main[0], P=main[1])
            pfr = VacuumRegion(pf[0], P=pf[1] - 1)
            matrices = matrix_calculate(sol, pfr)
            if save: # save the new VacuumRegions as pkl files
                if pfr_user is None:
                    # print('normal save')
                    sol.saveVacuumRegion('SOL_Vacuum.pkl')
                    pfr.saveVacuumRegion('PFR_Vacuum.pkl')
                else:
                    # print('user save')
                    sol.saveVacuumRegion('SOL_Vacuum.pkl')
                    pfr.saveVacuumRegion('PFR_userVacuum.pkl')
        elif type(sol) != VacuumRegion or type(pfr) != VacuumRegion: # provide a specific pkl file to read from
            # print('load from pkls')
            sol = VacuumRegion(sol)
            pfr = VacuumRegion(pfr)
            matrices = matrix_calculate(sol, pfr)
        else: # passing in a VacuumRegion object as sol and pfr arguments
            matrices = matrix_calculate(sol, pfr)
        
        # print('error here')
        self.sol = sol
        self.pfr = pfr
        print(self.pfr.numSurfaces)
        if sol_puff_dict is not None:
            sol_puff_dict['region'] = sol
        # print('matrix error')
        
        cftelematrix = matrices[0]
        cftelematrix_pf = matrices[1]
        puffing_matrix = matrices[2]

        dimension = len(cftelematrix) + 2
        cftelematrix_full = numpy.zeros((2, dimension, dimension, 6))
        cftelematrix_full[1, 1:-1, 1:-1, 0] = cftelematrix
        cftelematrix_full[1, 1:-1, 1:-1, 1] = cftelematrix

        if overwrite:

            bbb.cftelematrix[1, 1:-1 , 1:-1 , 0] = cftelematrix

            bbb.cftelematrix[0, 1:com.ixpt1[0]+1, 1:com.ixpt1[0]+1, 0] = cftelematrix_pf[:com.ixpt1[0], :com.ixpt1[0]] # 1
            bbb.cftelematrix[0, 1:com.ixpt1[0]+1, com.ixpt2[0]+1:com.nx+1, 0] = cftelematrix_pf[:com.ixpt1[0], com.ixpt1[0]:] # 2
            bbb.cftelematrix[0, com.ixpt2[0]+1:com.nx+1, 1:com.ixpt1[0]+1, 0] = cftelematrix_pf[com.ixpt1[0]:, :com.ixpt1[0]] # 3 
            bbb.cftelematrix[0, com.ixpt2[0]+1:com.nx+1, com.ixpt2[0]+1:com.nx+1, 0] = cftelematrix_pf[com.ixpt1[0]:, com.ixpt1[0]:] # 4
            
            self.save_matrices([cftelematrix, cftelematrix_pf, puffing_matrix], 
                            ['cftelematrix', 'cftelematrix_pf', 'puffing_matrix'], 
                            save_file)

            self.set('cftelematrix', cftelematrix_full)

        self.getue('isvacuummodel', cp=False)[0] = 1
        # bbb.isvacuummodel[0] = 1
        self.set('cfteleout', 1.0)

        if sol_puff_dict != None:
            # print('enter puff sol')
            main_puffing_array, main_puffing_location = self.puffing_array_calc(sol, sol_puff_dict['point'], sol_puff_dict['current'], puffing_matrix)
            self.main_puffing_array = numpy.transpose(main_puffing_array)
            if overwrite:
                self.save_matrices([self.main_puffing_array], ['puffing_array'], save_file)
            # a = self.plot_grid(sol, pfr, sol_test_surf=main_puffing_location)
            self.puffing_matrix = puffing_matrix
        
        def pump_helper(region, pump_dictionary, plot):
            from shapely import Polygon, intersects, Point
            # print('enter pump helper')
            pumping_surf = []
            points = []
            xbox = []
            ybox = []
            original_R = region.R_dictionary.copy()
            # print('Reflection coefficients:', region.R_dictionary)
            if region==self.pfr and (self.pfrPoints is not None) or ('box_coords' in pump_dictionary and len(pump_dictionary['box_coords']) > 3):
                # print('box error')
                for coord in pump_dictionary['box_coords']:
                    points.append(Point(coord))
                    xbox.append(coord[0])
                    ybox.append(coord[1])
                xbox.append(xbox[0])
                ybox.append(ybox[0])
                box = Polygon(points)

                for i in range(region.numSurfaces):
                    if i >= region.P:
                        seg = region.surfaces[i].segment
                        if intersects(seg, box):
                            region.R_dictionary[i] = pump_dictionary['albedo']
                            pumping_surf.append(i)
            else:
                print('pumping all walls')
                for i in range(region.numSurfaces):
                    if i >= region.P: # Pump on all wall surfaces
                        region.R_dictionary[i] = pump_dictionary['albedo']
                        pumping_surf.append(i)

            region.matrices()
            pumping_matrix_general = region.getOutputMatrix(region.AB_matrix, 1000000)
            pumping_matrix_general = numpy.transpose(pumping_matrix_general)

            if plot and len(pump_dictionary['box_coords']) > 3:
                # print('pump plot error')
                surfx = []
                surfy = []
                for i in pumping_surf:
                    surfx.append(region.surfaces[i].start.x)
                    surfy.append(region.surfaces[i].start.y)

                    surfx.append(region.surfaces[i].end.x)
                    surfy.append(region.surfaces[i].end.y)
                
                a = self.plot_grid(self.sol, self.pfr, label=False)
                import matplotlib.pyplot as plt
                ax = plt.gca()
                ax.plot(surfx, surfy, color='blue')
                ax.plot(xbox, ybox, color='gray')

            region.R_dictionary.update(original_R)
            
            return pumping_matrix_general

        if pf_pump_dict != None:
            # print('enter pump loop')
            self.box_loop = False
            pumping_matrix_pf = pump_helper(self.pfr, pf_pump_dict, pump_plot)
            self.pumping_matrix_pf = pumping_matrix_pf
            self.cftelematrix_pf = pumping_matrix_pf
            print(self.pfr.numSurfaces)
            print('pfr albedo:', pf_pump_dict['albedo'])
            a = self.plot_grid(sol, pfr, pf_test_surf=[], label=False)
            if overwrite:
                self.save_matrices([self.pumping_matrix_pf], ['cftelematrix_pf'], save_file_name=save_file)
                # print('bbb error')
                bbb.cftelematrix[0, 1:com.ixpt1[0]+1, 1:com.ixpt1[0]+1, 0] = self.pumping_matrix_pf[:com.ixpt1[0], :com.ixpt1[0]] # 1
                bbb.cftelematrix[0, 1:com.ixpt1[0]+1, com.ixpt2[0]+1:com.nx+1, 0] = self.pumping_matrix_pf[:com.ixpt1[0], com.ixpt1[0]:] # 2
                bbb.cftelematrix[0, com.ixpt2[0]+1:com.nx+1, 1:com.ixpt1[0]+1, 0] = self.pumping_matrix_pf[com.ixpt1[0]:, :com.ixpt1[0]] # 3 
                bbb.cftelematrix[0, com.ixpt2[0]+1:com.nx+1, com.ixpt2[0]+1:com.nx+1, 0] = self.pumping_matrix_pf[com.ixpt1[0]:, com.ixpt1[0]:] # 4

        if sol_pump_dict != None:

            pumping_matrix_sol = pump_helper(self.sol, sol_pump_dict, plot=False)
            self.pumping_matrix_sol = pumping_matrix_sol
            self.cftelematrix = pumping_matrix_sol
            print('sol albedo:', sol_pump_dict['albedo'])
            if overwrite:
                self.save_matrices([self.pumping_matrix_sol], ['cftelematrix'], save_file_name=save_file)
                bbb.cftelematrix[1, 1:-1 , 1:-1 , 0] = self.pumping_matrix_sol

        # print('completed')

    def plot_grid(self, sol, pfr, sol_plot=True, pfr_plot=True, sol_test_surf=[], pf_test_surf=[], label=False):
        if sol_plot and pfr_plot:
            m = sol.plotGeometry(labels=label, testsurf=sol_test_surf, showCircle=True)
            p = pfr.plotGeometry(labels=label, ax=m.get_axes()[0], testsurf=pf_test_surf, showCircle=True)
        elif sol_plot:
            m = sol.plotGeometry(labels=label, testsurf=sol_test_surf, showCircle=True)
        elif pfr_plot:
            p = pfr.plotGeometry(labels=label, testsurf=pf_test_surf, showCircle=True)

        # input("Press 'Enter' to close plots.")   

class VacuumRegion:
    def __init__(self, nodeList, P=0, variation=True, multiprocess=True, ncores=None, verbose=True):
        from shapely import Point, Polygon
        from tqdm import tqdm
        from pickle import load
        from multiprocessing import Process, Pipe, Manager
        from os import cpu_count, environ
        from itertools import islice
        from numpy import array, array_split
        from time import time
        from copy import deepcopy
    
        self.surfaces = {}

        starttime = time()
        if isinstance(nodeList, str):
            with open(nodeList, 'rb') as f:
                save = load(f)
                self.surfaces = save['surfaces']
                self.P = save['P']
        
        else:
            # Set up surfaces of geometry and the polygon object 
            self.P = P
            for i in range(len(nodeList)):
                startNode = Point(nodeList[i])
                if i == len(nodeList) - 1:
                    endNode = Point(nodeList[0])
                else:
                    endNode = Point(nodeList[i + 1])

                # Have plasma surfaces use a uniform dist., while wall surfaces use a cosine dist.
                if variation:
                    if i >= self.P: # non-plasma surfaces
                        offset = 1
                    else:
                        offset = 0
                    self.surfaces[i] = Surface((startNode.x, startNode.y), (endNode.x, endNode.y), i, r_offset=offset)
                else:
                    self.surfaces[i] = Surface((startNode.x, startNode.y), (endNode.x, endNode.y), i)
 
            # Create Polygon of Vacuum region for intersect checks
            self.geometry = Polygon(nodeList) 

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
        self.R_dictionary = {}
        for i in range(self.numSurfaces):
            if i >= self.P: # non-plasma surfaces
                self.R_dictionary[i] = 1
            else:
                self.R_dictionary[i] = 0

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

    def matrices(self):
        import numpy
        from numpy import zeros, identity, percentile, log
        from scipy.sparse import csr_array, block_array
        import seaborn as sns
        import matplotlib.pyplot as plt

        '''Creates R (self.R_matrix), C (self.C_matrix), A (self.A_matrix), 
            B (self.B_matrix), and AB (self.AB_matrix) matrices.'''

        # Array representations of R and C
        self.R_array = zeros((self.numSurfaces, self.numSurfaces))
        self.C_array = zeros((self.numSurfaces, self.numSurfaces))

        # Populate R array
        for surfaceID, rVal in self.R_dictionary.items():
            self.R_array[surfaceID][surfaceID] = rVal

        # Populate C array and take transpose
        for surfaceID, surface in self.surfaces.items(): # self.surfaces.items()
            for outputID in surface.neighbors.keys():
                self.C_array[surfaceID][outputID] = surface.neighbors[outputID]['flux']
        self.C_array = self.C_array.transpose()

        # R and C into sparse matrices
        self.R_matrix = csr_array(self.R_array)
        self.C_matrix = csr_array(self.C_array)

        # Zero and identity sparse matrices
        Zero_matrix = csr_array(zeros((self.numSurfaces, self.numSurfaces)))
        Identity_matrix = csr_array(identity(self.numSurfaces))

        # Create A, B, and AB sparse matrices
        self.A_matrix = block_array([[self.C_matrix, Zero_matrix], [Zero_matrix, Identity_matrix]])
        self.B_matrix = block_array([[self.R_matrix, Zero_matrix], [Identity_matrix - self.R_matrix, Identity_matrix]])

        # A * B
        self.AB_matrix = self.A_matrix @ self.B_matrix 


    def heatmapPlot(self):
        import numpy
        from numpy import zeros, identity, percentile, log
        from scipy.sparse import csr_array, block_array
        import seaborn as sns
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm

        '''Creates a heatmap of C, R, and Transport matrices.'''

        # Get A, B, AB
        self.matrices()

        # Generate output (transport) matrix
        self.getOutputMatrix(self.AB_matrix, 1000000)

        '''Print statements to check for unity of transport matrix.'''
        # print(f"sum(sum(output)): {sum(sum(self.output))}")
        # print(f"P: {self.P}")

        # Plotting heatmaps of C, R, and Output (Transport)
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        fig.suptitle("Cosine Distribution", fontsize=16)

        matricesToPlot = [self.C_array, self.R_array, self.output]
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

        return self.output

    def matrixPower(self, matrix, power):
        from numpy import zeros, identity
        from scipy.sparse import csr_array, block_array, linalg

        '''Raises a given matrix to the specified power.'''

        resultMatrix = linalg.matrix_power(matrix, power)

        return resultMatrix

    def getOutputMatrix(self, AB, power):
        from numpy import zeros, identity, transpose
        from scipy.sparse import csr_array, block_array

        AB_power_A = self.matrixPower(self.AB_matrix, power) @ self.A_matrix # (AB)^M * A

        # Final transport matrix
        self.output = zeros((self.P, self.P))

        gamma_array = zeros((self.numSurfaces * 2, 1))
        for i in range(0, self.P):
            gamma_array[i, 0] = 1

            rowCalculation = AB_power_A @ gamma_array

            gamma_array[i, 0] = 0

            gammaOut = rowCalculation[self.numSurfaces:]
            gammaFinal = gammaOut[0:self.P].flatten()

            for j in range(self.P):
                self.output[i, j] = gammaFinal[j]

        self.output = transpose(self.output)
        return self.output

    def getPuffingMatrix(self, reflections):
        from numpy import zeros, identity, transpose
        from scipy.sparse import csr_array, block_array
        '''Set-up matrix for puffing, gives surface to surface transport.'''

        # Control number of reflections
        puffingMatrix = self.matrixPower(self.AB_matrix, reflections) @ self.A_matrix

        return puffingMatrix
    
    def getPuffingArray(self, puffingMatrix, surfaceIndex, current):
        from numpy import zeros, identity, transpose
        from scipy.sparse import csr_array, block_array
        '''1-D array (geometric output flux vector) when puffing at one source surface given a
            source strength.'''

        gamma_array = zeros((self.numSurfaces * 2, 1))

        gamma_array[surfaceIndex, 0] = current

        rowCalculation = puffingMatrix @ gamma_array
        gammaOut = rowCalculation[self.numSurfaces:]
        puffingArray = gammaOut[0:self.P]

        return puffingArray
        

    def saveVacuumRegion(self, savename):
        from pickle import dump
        '''Use to save a Vacuum Region to avoid having to generate a new one every time. 
            Be sure to set pf/main (tokamakPlot), r_offset (Surface constructor), and 
            variation (Vacuum Region constructor).'''

        save = {
            'surfaces': self.surfaces,
            'P': self.P
        }
        with open(savename, 'wb') as f:
            dump(save, f)

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
            
        for _, surface in self.surfaces.items():
            color = 'k'
            if (surface.ID <  self.P):
                color = 'red'
            surface.plotSelf(color=color, ax=ax, label=labels, showCircle=showCircle)

        for itest in testsurf:
            self.surfaces[itest].plotConnections(
                            ax=ax, 
                            linewidth=connectionLineWidth,
                            **kwargs
            )
    
        for line in ax.lines:
            line.set_marker(".")

        ax.set_aspect('equal')
        ax.grid(False)
        plt.xlabel("R [m]")
        plt.ylabel("Z [m]")
    
        # plt.savefig('fullGeometry62.svg', dpi=300)
        
        plt.show(block=False)
    

        return ax.get_figure()
   
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
 


        







