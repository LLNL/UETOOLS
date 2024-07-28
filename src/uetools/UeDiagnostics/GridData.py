
class Grid():
    def __init__(self, case, flip=True, variables = ['te', 'ne', 'ni', 'ng']):
        from shapely import coverage_union_all, Polygon
        self.case = case
        self.geometry = self.case.get("geometry")[0].strip().lower().decode("UTF-8")
        self.cells = []
        self.case.about.set_speciesarrays()
        rm = self.case.get('rm')
        zm = self.case.get('zm')
        self.vars = {}
        for var in variables:
            self.vars[var] = self.case.get(var)

        if (self.geometry == "uppersn") and (flip is True):
            zm = self.case.disp - zm



        self.densities = {'e': {0: self.vars['ne']}}
        for i in range(len(case.about.ionarray)):
            ion = case.about.ionarray[i]
            species = (ion.split("+")[0]).lower()
            if species in ["d", "t"]:
                species = "h"
            try: 
                charge = int(ion.split("+")[1])
            except:
                if '+' not in ion:
                    # No charge: inertial neutral species
                    continue
                else:
                    charge = 1
            if species not in self.densities:
                self.densities[species] = {}
            self.densities[species][charge] = self.vars["ni"][:,:,i]
        for g in range(len(case.about.gasarray)):
            gas = case.about.gasarray[g]
            species = (gas.replace("0", "").replace("_2","")).lower()
            if species in ["d", "t"]:
                species = "h"
            mols = "_2" in gas
            if not mols:
                self.densities[species][0] = self.vars["ng"][:,:,g]
            else:
                self.densities[species]['mol'] = self.vars["ng"][:,:,g]




        mapdict = {}
        geometries = []
        # Create polygons for all cells
        (nx, ny, _) = rm.shape
        for ix in range(1,nx-1):
            if ix not in mapdict:
                mapdict[ix] ={}
            for iy in range(1,ny-1):
                verts = []
                for iv in [1, 2, 4, 3]:
                    verts.append([rm[ix, iy, iv], zm[ix, iy, iv]])
                self.cells.append(Cell(verts, (ix,iy)))#, self.vars)),
#                    case.ionarray, case.gasarray
#                ))
                # Create list of plygons for union
                geometries.append(self.cells[-1].polygon)
                # Map the coordinates to the appropriate Cell object
                mapdict[ix][iy] = self.cells[-1]
        # Create union cells, consisting of the plasma grid
        self.area = coverage_union_all(geometries)
        # Create and store Polygon for core
        corepts = []
        for nxpt in range(case.get('nxpt')):
            for ix in range(case.get('ixpt1')[nxpt]+1, case.get('ixpt2')[nxpt]+1):
                corepts.append([
                            case.get('rm')[ix, 0, 3], 
                            case.get('zm')[ix, 0, 3],
                ])
        self.core = Polygon(corepts)
        self.map = mapdict


    def plot_intensity(self, interval, ax=None,crm=None,zrange=(None,None), 
        mol=True, cbar=True,zscale=1
    ):
        from matplotlib.pyplot import subplots
        from numpy import array, nonzero, transpose, log10, floor
        from tqdm import tqdm
        from matplotlib.colors import Normalize,LogNorm
        from matplotlib.cm import ScalarMappable
        from matplotlib.pyplot import get_cmap,colorbar,figure

        if ax is None:
            f, ax = subplots()
            ret = True
        else:
            ret = False

        # Store all cell intensities in the requested interval to a list
        intensity = []
        for i in tqdm(range(len(self.cells))):
            cell = self.cells[i]
            if cell.crmeval is False:
                cell.calc_emissivity(crm)

            ci = [0,0]
            for i in nonzero(cell.H2_emission[1])[0]:
                ene = cell.H2_emission[0][i]
                if 10*ene > interval[0] and 10*ene < interval[1]:
                    ci[0] += cell.H2_emission[1][i]
            for i in nonzero(cell.H_emission[1])[0]:
                ene = cell.H_emission[0][i]
                if 10*ene > interval[0] and 10*ene < interval[1]:
                    ci[1] += cell.H_emission[1][i]
            intensity.append(ci)
        intensity = zscale*transpose(array(intensity))

        ind = 0
        if mol is False:
            ind = 1
        zmin, zmax = intensity[ind].min(), intensity[ind].max()
        if zrange[0] is not None:
            zmin=zrange[0]
        if zrange[1] is not None:
            zmax=zrange[1]
        Zcol=((log10(intensity[ind])-floor(log10(zmin)))/(floor(log10(zmax))-floor(log10(zmin))))
        cmap=get_cmap('magma')

        for i in range(len(self.cells)):
            cell = self.cells[i]
            xs, ys = cell.polygon.exterior.xy    
            col=cmap(Zcol[i])
            ax.fill(xs, ys, fc=col, ec='none')

        if cbar is True:
            norm = Normalize(vmin=floor(log10(zmin)),vmax=floor(log10(zmax)))
            norm = LogNorm(vmin=zmin,vmax=zmax)
            sm = ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar=colorbar(sm,ax=ax)
        ax.set_aspect('equal')
        return ax.get_figure()

    def plot_grid(self, ax=None):
        ''' Plots a polygrid of Patches '''
        from matplotlib.pyplot import subplots, Figure
        if ax is None:
            f, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]

        for cell in self.cells:
            cell.plot_cell(ax)

        self.case.plot.vessel(f)
        self.case.plot.plates(f)


        ax.set_aspect('equal')
        return ax.get_figure()
  





class Cell():
    ''' Containter object for grid cell info '''
    
    # TODO: create function calculating and storing the spectra
    def __init__(self, vertices, indices):
#, variables, ionarray=None,
#        gasarray=None
#    ):
        from shapely.geometry import Polygon
        self.vertices = vertices
        self.indices = indices
        self.polygon = Polygon(vertices)

                    
    def plot_cell(self, ax=None, linewidth = 0.05):
        if ax is None:
            f, ax = subplots()
        ax.plot(*self.polygon.exterior.xy, 'k-', linewidth=linewidth) 


