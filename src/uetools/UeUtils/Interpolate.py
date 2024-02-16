# NOTE
# The Grid objects set up here could be expanded upon to eventually
# completely control and drive/handle any UEDGE grids and/or gridding.
# The objects could read the grid data from an HDF5 file, calll the 
# necessary routines, etc. to handle grid modifications and interpolations

class Interpolate():

    def interpolate_snull(
        self, oldgrid, newgrid, oldsave=None, ishdf5=None, radtranspfile=None,
        newsavename=None, **kwargs
):
        """ Interpolates new solution based on previous state and new grid """
        from h5py import File

#       TODO: Make oldgrid option. If not supplied, use variables in memory.
#       TODO: Read and pass magnetic and physical data: calculate poloidal/
#           parallel distance between end-points of each flux-tube.
#           Start by using flux-tube averaged values? Will let one use 
#           RegularGridInterpolator. Then, try to use unstructured grid
#           interpolators instead.

        if ishdf5 is None:
            ishdf5 = self.get('isgriduehdf5')
        # Start by reading in new and old grid dimensions
        grid_old = {}
        grid_new = {}
        # Check if gridue is HDF5 file or not: read accordingly
        proplist = ['nxm', 'nym', 'ixpt1', 'ixpt2', 'iysptrx1', 'b', 'bpol', 'rm', 'zm']
        if ishdf5==True:
            with File(oldgrid) as f_oldgrid:
                for prop in proplist:
                    grid_old[prop] = f_oldgrid['grid/com'][prop][()]
            with File(newgrid) as f_newgrid:
                for prop in proplist:
                    grid_new[prop] = f_newgrid['grid/com'][prop][()]
        else:
            with open(oldgrid) as f_oldgrid:
                for i in range(len(proplist)):
                    gridprops = [int(x) for x in f_oldgrid.readline().split()]
                    grid_old[proplist[i]] = gridprops[i]
            with open(newgrid) as f_newgrid:
                    gridprops = [int(x) for x in f_newgrid.readline().split()]
                    grid_old[proplist[i]] = gridprops[i]


        # Record solution
        savedata = {}
        if oldsave is None:
            for var in ['nis', 'ngs', 'ups', 'tes', 'tis', 'tgs', 'phis']:
                savedata[var] = self.getue(var)
        else:
            with File(oldsave) as f_save:
                # Get handle to save group: UeCase save or not
                try:
                    f_save['restore']
                    savegroup = f_save['restore/bbb']
                except:
                    savegroup = f_save['bbb']
                for var in ['nis', 'ngs', 'ups', 'tes', 'tis', 'tgs', 'phis']:
                    savedata[var] = savegroup[var][()]

        # Record radial transport coefficients
        radtransp = {}
        varlist = [
            ['2D', ['dif_use', 'kye_use', 'kyi_use', 'tray_use']],
            ['1D', ['difniv', 'kyev', 'kyiv', 'travisv']],
            ['0D', ['difni', 'kye', 'kyi', 'travis']]
        ]
        if radtranspfile is None:
            for vartype, vars_ in varlist:
                radtransp[vartype] = {}
                for var in vars_:
                    radtransp[vartype][var] = self.getue(var)
        else:
            with File(radtranspfile) as f_radtransp:
                for vartype, vars_ in varlist:
                    radtransp[vartype] = {}
                    for var in vars_:
                        radtransp[vartype][var] = \
                            f_radtransp['diffusivities/bbb'][var][()]

        (save_nx, save_ny) = savedata['tis'].shape
        # Assert solution matches old grid
        if (save_nx != grid_old['nxm']+2) or (save_ny != grid_old['nym']+2):
            raise ValueError('Save file dimensions do not match the original'
                'grid dimensions. Grid (nx, ny) = ({}, {}), save (nx, ny) = '
                '({}, {})'.format(nx_old, ny_old, save_nx-2, save_ny-2))
    
        self.oldgrid = GridSnull(grid_old, savedata, radtransp, **kwargs)
        newgrid = self.oldgrid.interpolate_grid(grid_new, **kwargs)
        newsave = newgrid.savedata
        newradtransp = newgrid.radtransp
        if newsavename is None:
            newsavename = 'interpolated_{}x{}-{}x{}'.format(grid_old['nxm'],
                grid_old['nym'], grid_new['nxm'], grid_new['nym'])

        with File(newsavename, 'w') as f_save:
            f_save.create_group('bbb')
            for var in ['nis', 'ngs', 'ups', 'tes', 'tis', 'tgs', 'phis']:
                f_save['bbb'].create_dataset(var, data=newsave[var])
            f_save.create_group('diffusivities')
            f_save['diffusivities'].create_group('bbb')
            for _, varlist in newradtransp.items():
                for var, data in varlist.items():
                    f_save['diffusivities/bbb'].create_dataset(var, \
                        data=data)

        
        return newgrid

    # TODO
    # Try using an unstructured grid for interpolation:
    # The y-dimension is radial distance from the core in each 'patch', 
    # the x-dimension is the poloidal distance from the target, 
    # centered at 'cuts'. Scheme can be upgraded to include magnetics,
    # and calculating parallel distances rather than poloidal

    def gridmorph(self, oldgrid, newgrid, **kwargs):
        """ Uses the continuation solver to morph the UEDGE grids """
        from copy import deepcopy

        manualgrid = deepcopy(self.get('manualgrid'))
        self.setue('gridmorph', 0)
        self.setue('manualgrid', 1)
        self.continuation_solve(
                "gridmorph", 
                1, 
                commands=[f"self.morphed_mesh('{oldgrid}','{newgrid}',"+\
                        "self.getue('gridmorph'), standalone=False)",
                ],
                newgeo=True, 
                **kwargs
        )
        self.setue('manualgrid', manualgrid)
        
        

    def morphed_mesh(self, oldgrid, newgrid, fraction, 
                standalone=True, reread=False
        ):
        """ Morphs the physical and magnetic mesh between old and new grids 

        Grids must have same number of poloidal and radial grid cells.
        The number of cells in each macro-region should be identical.
        First use interpolation routines to interpolate solution 
        to new grid.
        """
        from h5py import File
        from Forthon import packageobject
        from copy import deepcopy
        
        
        # List of variables to morph
        variables = ["rm", "zm", "psi", "br", "bz", "bpol", "bphi", "b"]
        # "nlim", "xlim", "ylim", "nplate1", "nplate2", "rplate1", "rplate2",
        # "zplate1", "zplate2"
        if ("griddata"  not in dir(self)) or reread:
            self.griddata = {'old': {}, 'new': {}, 'delta': 0}
            for var in variables:
                self.griddata['old'][var] = self.hdf5search(oldgrid, var)
                self.griddata['new'][var] = self.hdf5search(newgrid, var)
        for var in variables:
            self.setue(var, (1-fraction)*self.griddata['old'][var] \
                    + fraction*self.griddata['new'][var])
        self.griddata['delta'] = fraction
        packageobject('bbb').__getattribute__('guardc')()
        if standalone:
            iprint = deepcopy(self.get('iprint'))
            self.setue("iprint", 0)
            self.populate(silent=True, verbose=False)
            self.setue("iprint", iprint)
        



class GridSnull:
    """ Object containing single null grid data for interpolation """
    def __init__(self, griddata, savedata, radtransp, interpolation='index'):
        self.nx = griddata['nxm']
        self.ny = griddata['nym']
        self.ixpt1 = griddata['ixpt1']
        self.ixpt2 = griddata['ixpt2']
        self.iysptrx = griddata['iysptrx1']
        self.savedata = savedata
        self.radtransp = radtransp
        # Define the snull topology:
        #       ____________________________________________
        #      |            |                   |           |
        # SOL  |            |                   |           |
        #      |____________|___________________|___________|
        #      |            |                   |           |
        # CORE |            |                   |           |
        #      |____________|___________________|___________|
        #           ILEG             CORE            OLEG

        if interpolation.lower() == 'index':
            PatchType = IndexGridPatch
            self.connlen = None
        elif interpolation.lower() == 'parallel':
            PatchType = ParallelGridPatch
            self.connlen = self.calc_connlen(griddata)


        else:
            raise ValueError('Interpolation type "{interpolation}" not recognized!')
            
        self.patches = {
            'sol': {
                'ileg': PatchType(0, self.ixpt1+1, self.iysptrx+1, 
                            self.ny+2, savedata, radtransp, self.connlen),
                'core': PatchType(self.ixpt1+1, self.ixpt2+1, self.iysptrx+1, 
                            self.ny+2, savedata, radtransp, self.connlen),
                'oleg': PatchType(self.ixpt2+1, self.nx+2, self.iysptrx+1, 
                            self.ny+2, savedata, radtransp, self.connlen),
            },
            'core': {
                'ileg': PatchType(0, self.ixpt1+1, 0, 
                            self.iysptrx+1, savedata, radtransp, self.connlen),
                'core': PatchType(self.ixpt1+1, self.ixpt2+1, 0, 
                            self.iysptrx+1, savedata, radtransp, self.connlen),
                'oleg': PatchType(self.ixpt2+1, self.nx+2, 0, 
                            self.iysptrx+1, savedata, radtransp, self.connlen),
            }
        }

    def interpolate_grid(self, griddata, **kwargs):
        """ Interpolates the current grid to a new GridSnull object """ 
        from numpy import concatenate
        nx_new = griddata['nxm']
        ny_new = griddata['nym']
        ixpt1_new = griddata['ixpt1']
        ixpt2_new = griddata['ixpt2']
        iysptrx_new = griddata['iysptrx1']
        connlen_new = self.calc_connlen(griddata)

        
        
        # Store the new dimensions
        nxny = {
            'sol':{
#                'common': [nx_new+2, ny_new-iysptrx_new+1]
                'ileg': [0, ixpt1_new+1, iysptrx_new+1, ny_new+2, connlen_new],
                'core': [ixpt1_new+1, ixpt2_new+1, iysptrx_new+1, ny_new+2, connlen_new],
                'oleg': [ixpt2_new+1, nx_new+2, iysptrx_new+1, ny_new+2, connlen_new]
            },
            'core':{
                'ileg': [0, ixpt1_new+1, 0, iysptrx_new+1, connlen_new],
                'core': [ixpt1_new+1, ixpt2_new+1, 0, iysptrx_new+1, connlen_new],
                'oleg': [ixpt2_new+1, nx_new+2, 0, iysptrx_new+1, connlen_new]
            }
        }

        patches_new = {}
        for radialkey, poloidalpatches in self.patches.items():
            try: 
                patches_new[radialkey]
            except:
                patches_new[radialkey] = {}
            for poloidalkey, patch in poloidalpatches.items():
                try:
                    patches_new[radialkey][poloidalkey]
                except:
                    patches_new[radialkey][poloidalkey] = \
                        patch.interpolate_solution(*nxny[\
                        radialkey][poloidalkey])

        
        savedata_new = {}
        for variable in self.savedata.keys():
            savedata_new[variable] = concatenate((
                    concatenate((
                        patches_new['core']['ileg'][variable], 
                        patches_new['core']['core'][variable],
                        patches_new['core']['oleg'][variable]
                    )), 
                    concatenate((
                        patches_new['sol']['ileg'][variable], 
                        patches_new['sol']['core'][variable],
                        patches_new['sol']['oleg'][variable]
                    ))
                ), axis = 1)
         
        radtransp_new = {}
        for vartype, varlist in self.radtransp.items():
            radtransp_new[vartype] = {}
            for variable in varlist.keys():
                try:
                    radtransp_new[vartype][variable] = concatenate((
                            concatenate((
                                patches_new['core']['ileg'][variable], 
                                patches_new['core']['core'][variable],
                                patches_new['core']['oleg'][variable]
                            )), 
                            concatenate((
                                patches_new['sol']['ileg'][variable], 
                                patches_new['sol']['core'][variable],
                                patches_new['sol']['oleg'][variable]
                            ))
                        ), axis = 1)   
                except:
                    radtransp_new[vartype][variable] = \
                        self.radtransp[vartype][variable]
        
        # TODO:
        # Add interpolation in 1D
                    
        return GridSnull(griddata, savedata_new, radtransp_new, **kwargs)

    def calc_connlen(self, griddata):
        ''' Calculates the parallel connection length '''
        from copy import deepcopy
        from numpy import zeros
        # Calculate the parallel length of each cell
        # Pad with negligible value to avoid zero-divisor
        rbfbt = (deepcopy(griddata['b'])/(deepcopy(griddata['bpol'])+1e-50))[:,:,0]
        # And again for future use
        rbfbt[rbfbt==0] = 1
        rm = griddata['rm']
        zm = griddata['zm']
        (nxg, nyg, _) = rm.shape
        midpoints = zeros((2, nxg+1, nyg))
        # (R,Z)-coordinates of the midpoints
        midpoints[0,:-1] = 0.5*(rm[:,:,1]+rm[:,:,3])
        midpoints[1,:-1] = 0.5*(zm[:,:,1]+zm[:,:,3])
        # Do right boundary
        midpoints[0,-1] = 0.5*(rm[-1,:,2]+rm[-1,:,4])
        midpoints[1,-1] = 0.5*(zm[-1,:,2]+zm[-1,:,4])
        # Calculate poloidal distance for each cell
        polloc = midpoints[:,1:] - midpoints[:,:-1]
        pollength = (polloc[0]**2 + polloc[1]**2)**0.5
        # Caclulate the parallel length of each cell
        parlength = pollength * rbfbt
        # Calculate parallel distance from left boundary
        connlen = zeros((nxg, nyg))
        for i in range(1,nxg):
            connlen[i] += connlen[i-1] + parlength[i]
        # Set guard cells
        connlen[:,0] = connlen[:,1]
        connlen[:,-1] = connlen[:,-2]
        return connlen

            
        

class IndexGridPatch:
    """ Object containing data for topological patch """
    # TODO: Radial interpolation according to midplane/PFR PSIN?
    def __init__(self, nxl, nxu, nyl, nyu, savedata, radtransp, *args, **kwargs):
        """ Set up the required interpolators """
        from numpy import linspace
        from scipy.interpolate import RegularGridInterpolator
        from copy import deepcopy
        self.nxl = nxl
        self.nxu = nxu
        self.nx = nxu - nxl # Number of nodes in X-direction of patch
        self.nyl = nyl
        self.nyu = nyu
        self.ny = nyu - nyl # Number of nodes in Y-direction of patch
        self.savedata = deepcopy(savedata)
        self.radtransp = deepcopy(radtransp)
        # Create linearly distributed points
        self.x = linspace(0, 1, self.nx) 
        self.y = linspace(0, 1, self.ny) 
        # Break data into patch and create interpolation functions
        self.interp = {}
        for variable, data in self.savedata.items():
            self.savedata[variable] = data[nxl:nxu, nyl:nyu]
            if len(data.shape) == 3:
                nz = data.shape[2]
                # Do interpolation with multiple species
                z = linspace(0, 10*(nz-1), nz)
                self.interp[variable] = [RegularGridInterpolator( \
                    (self.x, self.y, z), self.savedata[variable]), nz]
            else:
                # Do interpolation with single species
                self.interp[variable] = RegularGridInterpolator( \
                    (self.x, self.y), self.savedata[variable])

        for variable, data in self.radtransp['2D'].items():
            self.radtransp[variable] = data[nxl:nxu, nyl:nyu]
            if len(data.shape) == 3:
                nz = data.shape[2]
                # Do interpolation with multiple species
                z = linspace(0, 10*(nz-1), nz)
                self.interp[variable] = [RegularGridInterpolator( \
                    (self.x, self.y, z), self.radtransp[variable]), nz]
            else:
                # Do interpolation with single species
                self.interp[variable] = RegularGridInterpolator( \
                    (self.x, self.y), self.radtransp[variable])

           

    def interpolate_solution(self, nxl, nxu, nyl, nyu, *args):
        """ Returns an interpolation of solution in index space """
        from numpy import linspace, meshgrid
        new_solution = {}
        nx = nxu - nxl
        ny = nyu - nyl
        x = linspace(0, 1, nx)
        y = linspace(0, 1, ny)
        X, Y = meshgrid(x, y, indexing='ij')
        for variable, interpolator in self.interp.items():
            # Single-dimensioned variable
            try:
                new_solution[variable] = interpolator((X, Y))
            # Multi-species variable
            except:
                z = linspace(0, 10*(interpolator[1]-1), interpolator[1])
                Xz, Yz, Z = meshgrid(x, y, z, indexing='ij')
                new_solution[variable] = interpolator[0]((Xz, Yz, Z))
        return new_solution
        

class ParallelGridPatch:
    """ Object containing data for topological patch """
    def __init__(self, nxl, nxu, nyl, nyu, savedata, radtransp, connlen):
        """ Set up the required interpolators """
        from numpy import linspace, zeros
        from scipy.interpolate import griddata
        from copy import deepcopy
        self.nxl = nxl
        self.nxu = nxu
        self.nx = nxu - nxl # Number of nodes in X-direction of patch
        self.nyl = nyl
        self.nyu = nyu
        self.ny = nyu - nyl # Number of nodes in Y-direction of patch
        self.savedata = deepcopy(savedata)
        self.radtransp = deepcopy(radtransp)
        self.connlen = deepcopy(connlen)

        self.xy = self.get_xy(nxl, nxu, nyl, nyu, 
            self.nx, self.ny, self.connlen)
        # Break data into patch and create interpolation functions
        self.interp = {}
        for variable, data in self.savedata.items():
            self.savedata[variable] = data[nxl:nxu, nyl:nyu]
            self.interp[variable] = []
            if len(data.shape) == 3:
                nz = data.shape[2]
                for iz in range(nz):
                # Append list [ (x,y), z ] containing data necessary for
                # interpolation 
                    self.interp[variable].append([
                            self.xy, self.savedata[variable][:,:,iz].ravel()
                        ]
                    )
            else:
                # Do interpolation with single species
                self.interp[variable].append([
                        self.xy, self.savedata[variable].ravel()
                    ]
                )

        for variable, data in self.radtransp['2D'].items():
            self.radtransp[variable] = data[nxl:nxu, nyl:nyu]
            self.interp[variable] = []
            if len(data.shape) == 3:
                nz = data.shape[2]
                # Do interpolation with multiple species
                for iz in range(nz):
                    self.interp[variable].append([
                            self.xy, self.radtransp[variable][:,:,iz].ravel()
                        ]
                    )
            else:
                # Do interpolation with single species
                self.interp[variable].append([
                        self.xy, self.radtransp[variable].ravel()
                    ]
                )
           

    def interpolate_solution(self, nxl, nxu, nyl, nyu, connlen):
        """ Returns an interpolation of solution in index space """
        from numpy import linspace, meshgrid, array, mean
        from scipy.interpolate import griddata
        from copy import deepcopy
        new_solution = {}
        nx = nxu - nxl
        ny = nyu - nyl
        connlen = deepcopy(connlen)
        newpoints = self.get_xy(nxl, nxu, nyl, nyu, nx, ny, connlen)
        for variable, interpolator in self.interp.items():
            new_solution[variable] = []
            for i in range(len(interpolator)):
                new_solution[variable].append(
                    griddata(
                        interpolator[i][0], interpolator[i][1], newpoints,
                    ).reshape((nx, ny))
                )
        for variable, solution in new_solution.items():
            if variable in ['tes', 'tis', 'phis']:
                new_solution[variable] = new_solution[variable][0]
            else:
                new_solution[variable] = array(new_solution[variable]).\
                    transpose((1,2,0))
            # Pop into arrays according to shape
        return new_solution
        

    def get_xy(self, nxl, nxu, nyl, nyu, nx, ny, connlen):
        from numpy import linspace, zeros, concatenate
        from copy import deepcopy
        # Truncate connection length to correspond to patch
        x = connlen[nxl:nxu, nyl:nyu]
        # Check whether we are starting from the right
        xzero = deepcopy(x[0])
        for i in range(x.shape[0]):
            x[i] -= xzero # Make sure left bound starts at zero
            x[i] /= x[-1] # Make sure right boundary ends at 1
        y = zeros((nx, ny))
        for i in range(y.shape[0]):
            y[i] = linspace(0, 1, ny) 
        return x.ravel(), y.ravel()

