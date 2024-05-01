
class Grid():
    def __init__(self, case, flip=True, variables = ['te', 'ne', 'ni', 'ng']):
        self.case = case
        self.geometry = self.case.get("geometry")[0].strip().lower().decode("UTF-8")
        self.cells = []
        self.case.set_speciesarrays()
        rm = self.case.get('rm')
        zm = self.case.get('zm')
        self.vars = {}
        for var in variables:
            self.vars[var] = self.case.get(var)

        if (self.geometry == "uppersn") and (flip is True):
            zm = self.case.disp - zm
        # Create polygons for all cells
        (nx, ny, _) = rm.shape
        for ix in range(1,nx-1):
            for iy in range(1,ny-1):
                verts = []
                for iv in [1, 2, 4, 3]:
                    verts.append([rm[ix, iy, iv], zm[ix, iy, iv]])
                self.cells.append(Cell(verts, (ix,iy), self.vars,
                    case.ionarray, case.gasarray
                ))


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

        self.case.plotvessel(f)
        self.case.plotplates(f)


        ax.set_aspect('equal')
        return ax.get_figure()
  





class Cell():
    ''' Containter object for grid cell info '''
    
    # TODO: create function calculating and storing the spectra
    def __init__(self, vertices, indices, variables, ionarray=None,
        gasarray=None
    ):
        from shapely.geometry import Polygon
        self.vertices = vertices
        self.indices = indices
        self.ionarray = ionarray
        self.gasarray = gasarray
        self.polygon = Polygon(vertices)
        for var, value in variables.items():
            self.__setattr__(var, variables[var][*indices])
        self.emission = {}

        self.nh = self.ng[0]
        for species in self.gasarray:
            species = species.split("0")[0].lower()
            if '2' in species:
                continue
            self.__setattr__(f"n{species}", [])
            for i in range(len(self.gasarray)):
                if (species in self.gasarray[i].lower()) and \
                    ('2' not in self.gasarray[i]) and \
                    (species not in ['h','d','t']):
                    self.__getattribute__(f"n{species}").append(self.ng[i])
            for i in range(len(self.ionarray)):
                if species in self.ionarray[i].lower():
                    self.__getattribute__(f"n{species}").append(self.ni[i])
                    
    def plot_cell(self, ax=None, linewidth = 0.05):
        if ax is None:
            f, ax = subplots()
        ax.plot(*self.polygon.exterior.xy, 'k-', linewidth=linewidth) 

    def set_emission(self, adasspecies, lam=None, chargestate=None,
        rtype = ['excit', 'recom', 'chexc']):
        species = adasspecies.species.lower()
        self.emission[species] = adasspecies.calc_emission(
            self.ne, self.te, self.__getattribute__(f"n{species}"), self.nh, 
            lam, chargestate, rtype)
        
        


