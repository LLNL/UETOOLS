# NOTE
# The Grid objects set up here could be expanded upon to eventually
# completely control and drive/handle any UEDGE grids and/or gridding.
# The objects could read the grid data from an HDF5 file, calll the 
# necessary routines, etc. to handle grid modifications and interpolations

class Interpolate():

    def interpolate_snull(
        self, oldgrid, newgrid, oldsave=None, hdf5=None, newsavename=None, **kwargs
):
        """ Interpolates new solution based on previous state and new grid """
        from h5py import File

        if hdf5 is None:
            hdf5 = self.get('isgriduehdf5')
        # Start by reading in new and old grid dimensions
        grid_old = {}
        grid_new = {}
        # Check if gridue is HDF5 file or not: read accordingly
        proplist = ['nxm', 'nym', 'ixpt1', 'ixpt2', 'iysptrx1']
        if hdf5==True:
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

        # Read solution
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

        (save_nx, save_ny) = savedata['tis'].shape
        # Assert solution matches old grid
        if (save_nx != grid_old['nxm']+2) or (save_ny != grid_old['nym']+2):
            raise ValueError('Save file dimensions do not match the original'
                'grid dimensions. Grid (nx, ny) = ({}, {}), save (nx, ny) = '
                '({}, {})'.format(nx_old, ny_old, save_nx-2, save_ny-2))
    
        newgrid = GridSnull(grid_old, savedata).interpolate_grid(grid_new)
        newsave = newgrid.savedata
        if newsavename is None:
            newsavename = 'interpolated_{}x{}-{}x{}'.format(grid_old['nxm'],
                grid_old['nym'], grid_new['nxm'], grid_new['nym'])

        with File(newsavename, 'w') as f_save:
            f_save.create_group('bbb')
            for var in ['nis', 'ngs', 'ups', 'tes', 'tis', 'tgs', 'phis']:
                f_save['bbb'].create_dataset(var, data=newsave[var])
        
        return newgrid



    # TODO
    # Try using an unstructured grid for interpolation:
    # The y-dimension is radial distance from the core in each 'patch', 
    # the x-dimension is the poloidal distance from the target, 
    # centered at 'cuts'. Scheme can be upgraded to include magnetics,
    # and calculating parallel distances rather than poloidal

class GridSnull:
    """ Object containing single null grid data for interpolation """
    def __init__(self, dimensions, savedata):
        self.nx = dimensions['nxm']
        self.ny = dimensions['nym']
        self.ixpt1 = dimensions['ixpt1']
        self.ixpt2 = dimensions['ixpt2']
        self.iysptrx = dimensions['iysptrx1']
        self.savedata = savedata
        # Define the snull topology:
        #       ____________________________________________
        #      |                                            |
        # SOL  |                    COMMON                  |
        #      |____________________________________________|
        #      |            |                   |           |
        # CORE |    ILEG    |        CORE       |   OLEG    |
        #      |____________|___________________|___________|
        # 
        self.patches = {
            'sol': {
                'common': IndexGridPatch(0, self.nx+2, self.iysptrx+1, 
                            self.ny+2, savedata)
            },
            'core': {
                'ileg': IndexGridPatch(0, self.ixpt1+1, 0, 
                            self.iysptrx+1, savedata),
                'core': IndexGridPatch(self.ixpt1+1, self.ixpt2+1, 0, 
                            self.iysptrx+1, savedata),
                'oleg': IndexGridPatch(self.ixpt2+1, self.nx+2, 0, 
                            self.iysptrx+1, savedata),
            }
        }

    def interpolate_grid(self, dimensions):
        """ Interpolates the current grid to a new GridSnull object """ 
        from numpy import concatenate
        nx_new = dimensions['nxm']
        ny_new = dimensions['nym']
        ixpt1_new = dimensions['ixpt1']
        ixpt2_new = dimensions['ixpt2']
        iysptrx_new = dimensions['iysptrx1']
        
        # Store the new dimensions
        nxny = {
            'sol':{
                'common': [nx_new+2, ny_new-iysptrx_new+1]
            },
            'core':{
                'ileg': [ixpt1_new+1, iysptrx_new+1],
                'core': [ixpt2_new-ixpt1_new, iysptrx_new+1],
                'oleg': [nx_new-ixpt2_new+1 , iysptrx_new+1]
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
                    patches_new['sol']['common'][variable]
                ), axis = 1)
            
        
                    
        return GridSnull(dimensions, savedata_new)
            
        

class IndexGridPatch:
    """ Object containing data for topological patch """
    def __init__(self, nxl, nxu, nyl, nyu, savedata):
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

    def interpolate_solution(self, nx, ny):
        """ Returns an interpolation of solution in index space """
        from numpy import linspace, meshgrid
        new_solution = {}
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
        

