class ChordLine():
    """ Creates a single pencil-beam chord object """
    def __init__(self, p0, p1, grid, *args, res = 500, **kwargs):
        """ Creates a generalized pencil-beam

        Pencil-beams are used for 1D line-integration and can accommodate
        tangential view-chords.

        Arguments:
        ==========
        p0 - start point of chord, Shapely Point object
        p1 - end point of chord, Shapely Point object
        grid - UETOOLS Grid object used to link chord to geometry

        Keyword arguments:
        =================
        res [500] - Number of segments to split the chord into
        """
        from shapely.geometry import Point, LineString
        from numpy import array, pi, linspace
        from scipy.interpolate import interp1d

        self.p0 = p0 # Chord starting point
        self.p1 = p1 # Chord end point
        # Catch edge-case of purely horizontal tangential view-chord
        if abs(self.p0.z - self.p1.z) < 1e-5:
            raise Exception("Tangential lines in the horizontal plane not yet implemented")
        # Get the total vector length 
        self.length = p0.distance(p1)
        # Parametrize the lines as a function of path-length along chord 
        # f(L) = x
        self.fx = interp1d([0, self.length], [self.p0.x, self.p1.x])
        # f(L) = y
        self.fy = interp1d([0, self.length], [self.p0.y, self.p1.y])
        # f(L) = z
        self.fz = interp1d([0, self.length], [self.p0.z, self.p1.z])
        # f(z) = L
        self.fL = interp1d([self.p0.z, self.p1.z], [0, self.length])
        # TODO: handle exception when chord is horizontal!
        def generator(L):
            return ( 
                ( self.fx(L)**2 + self.fy(L)**2)**0.5,
                self.fz(L)
            )
        # Create a res long set of poloidal coordinates along the chord 
        points = [generator(0)]
        for chunk in linspace(0, self.length, res):
            points.append(generator(chunk))
        # Ceate Shapely LineString object of polidal projection of chord 
        self.chord = LineString(points)
        # Get start- and end-points in poloidal plane
        p0rz = Point(points[0])
        p1rz = Point(points[-1])
        # Find the cells intersected by the chord and store 
        # their path-lengths and indices
        self.cells = {}
        # Loop through all cells to find the cells intersected by the chord
        for ix, row in grid.map.items():
            for iy, cell in row.items():
                # Check whether the chord starts/ends within the grid and raise errors
                # TODO: implement handling of chords starting/terminating within grid?
                # TODO: raise Exception if chord starts in core?
                if cell.polygon.contains(p0rz):
                    raise ValueError("Chord start point inside grid")
                elif cell.polygon.contains(p1rz):
                    raise ValueError("Chord end point inside grid")
                buff = []
                # See whether chord intersect (ix, iy) grid-cell
                if grid.map[ix][iy].polygon.intersects(self.chord):
                    # Store the start- and end-points in paramtrized chord-distance
                    for intersect in grid.map[ix][iy].polygon.exterior.intersection(self.chord).geoms:
                        # Use the polidal Z-coordinate for back-projection
                        # NOTE: this back-projection is what limits the routine
                        # to non-horizontal tangential chords only
                        buff.append(self.fL(intersect.y))
                    # Set up objects if this is the first chord intersection 
                    # with the cell
                    if ix not in self.cells:
                        self.cells[ix] = {}
                    if iy not in self.cells[ix]:
                        self.cells[ix][iy] = {
                            'dL': 0
                        }
                    # Store the path-length through the cell
                    buff.sort()
                    self.cells[ix][iy]['dL'] += buff[1] - buff[0]

    def plot_chord(self, ax=None, color='r', linewidth=0.5):
        """ Plots poloidal projection of the chord object
        
        Keyword arguments:
        ==================
        ax [None] - axis to plot chord on. If none, creates a new figure
        color ['r'] - color of chord
        linewidth [0.5] - width of chord

        Returns:
        ========
        Figure object
        
        """
        from matplotlib.pyplot import Figure, subplots
        # Check whether axis is specified, create a new fig if not 
        if ax is None:
            f, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        # Plot the polidal projection
        ax.plot(*self.chord.xy, '-', color=color, linewidth=linewidth)
        # Ensure Figure object is returned
        return ax.get_figure()
        
    def integrate_field(self, field):
        """ Returns the path-integral of chord through field

        Arguments:
        ==========
        field - 2D array of shape compatible with UEDGE setup containing 
                field values

        Returns:
        ========
        Float: line-integral of chord through field in [field]*m 
        """
        integral = 0
        # Loop through all cells intersected by chord
        for ix, row in self.cells.items():
            for iy, L in row.items():
                # Access the correct cell data from field and perform
                # numerical line-integral
                integral += field[ix,iy]*L['dL']
        return integral


class ChordPoly():
    ''' Creates a 2D poloidal chord object '''
    def __init__(self, p0, p1, grid, width=None, omega=None):
        """ Creates a generalized pencil-beam

        Pencil-beams are used for 1D line-integration and can accommodate
        tangential view-chords.

        Arguments:
        ==========
        p0 - start point of chord, Shapely Point object
        p1 - end point of chord, Shapely Point object
        grid - UETOOLS Grid object used to link chord to geometry

        Keyword arguments:
        =================
        width [None] - Width of beam at end-point in m. Calculated if 
                       omega is set. Assumed to be 1e-9m if both width
                       and omega are None (approximate pencil-beam).
        omega [None] - Solid-angle of detector at start-point. Calculated
                       if width is set. 

        """
        from shapely.geometry import Point
        from numpy import array, pi
        self.p0 = p0 # Chord starting point
        self.p1 = p1 # Chord end point
        # Get total chord length
        self.length = self.p0.distance(self.p1)
        # Set up the 2D polygon triangle approximating the chord 
        # view geomtery
        if (width is None) and (omega is None):
            # Approximate pencil beam
            self.width = 1e-9 # Pencil beam
            self.omega = self.calc_solidangle(self.width)
        elif (width is None) and (omega is not None):
            self.omega = omega
            self.width = self.calc_width(self.omega)
        elif (width is not None) and (omega is None):
            self.width = width
            self.omega = self.calc_solidangle(self.width)
        else:
            raise ValueError("Set either width or omega, not both")

        self.create_poly()
        self.cells = {}
        self.cell_polys = []

        for ix, row in grid.map.items():
            for iy, cell in row.items():
                if cell.polygon.contains(self.p0):
                    raise ValueError("Chord start point inside grid")
                elif cell.polygon.contains(self.p1):
                    raise ValueError("Chord end point inside grid")
                elif self.poly.intersects(cell.polygon):
                    if ix not in self.cells:
                        self.cells[ix] = {}
                    if iy not in self.cells[ix]:
                        self.cells[ix][iy] = {
                            'dL': 0
                        }
                    self.cells[ix][iy]['dL'] += self.intersected_dL(cell.polygon)
                    if self.poly.contains(cell.polygon):
                        self.cells[ix][iy]['dL'] = self.contained_dL(cell.polygon)
                    
                    

    def calc_width(self, omega):
        ''' Calculates the 2D 2*r at the end-point for a given solid angle '''
        # Calculates 
        from numpy import pi
        return self.length*(2/(2*pi-omega))*(4*pi*omega-omega**2)**0.5

    def calc_solidangle(self, w):
        ''' Calculates the solid angle for a give 2D width 2*r at the end-point '''
        from numpy import pi
        return 2*pi*(1-1/(1+w**2/4/self.length**2)**0.5)


    def contained_dL(self, cell):
        # Add the entire cell polygon to the list of grid cells
        # inside the LOS polygon:

        # Calculate dL for the integration by postulating a
        # rectangular cell orthogonal to the line-of-sight with
        # area of the original cell and width determined by the
        # distance of the cell from the origin of the LOS and
        # the width of the LOS at its end:
        L_cell = cell.centroid.distance(self.p0)
        w_orthog = L_cell/self.length*self.width
        self.cell_polys.append(cell)
        return cell.area/w_orthog

    def intersected_dL(self, cell):
        # Determine the part of the grid cell polygon that
        # intersects with the line-of-sight polygon:
        clipped = self.poly.intersection(cell)

        # Calculate dL for the integration by postulating a
        # rectangular cell orthogonal to the line-of-sight with
        # area of the original cell and width determined by the
        # distance of the cell from the origin of the LOS and
        # the width of the LOS at its end: 
        L_cell = clipped.centroid.distance(self.p0)
        w_orthog = L_cell/self.length*self.width
        self.cell_polys.append(clipped)
        return clipped.area/w_orthog
        
    
    def create_poly(self):
        from numpy import sqrt, cos, sin, arctan
        from shapely.geometry import Polygon

        denom = self.p1.y - self.p0.y
        # Avoid zero-divisors
        if denom == 0:
            denom = 1e-10
        theta = arctan((self.p1.x - self.p0.x)/denom)
        
        # Add the elongated line-of-sight to the list of LOS polygons: 
        self.poly = Polygon([
            (self.p0.x, self.p0.y), 
            (self.p1.x - 0.5*self.width*cos(theta), 
                                self.p1.y + 0.5*self.width*sin(theta)), 
            (self.p1.x + 0.5*self.width*cos(theta), 
                                self.p1.y - 0.5*self.width*sin(theta))
            ]
        )
 
    def integrate_field(self, field):
        integral = 0
        for ix, row in self.cells.items():
            for iy, L in row.items():
                integral += field[ix,iy]*L['dL']
        return integral


    def plot_chord(self, ax=None, poly=True, line=True, color='r', alpha=0.2):
        from matplotlib.pyplot import Figure, subplots
        if ax is None:
            f, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        if line is True:
            ax.plot([self.p0.x, self.p1.x], [self.p0.y, self.p1.y], '-', 
                linewidth=0.5,color=color)
        if poly is True:
             xs, ys = self.poly.exterior.xy    
             ax.fill(xs, ys, alpha=alpha, fc=color, ec='none')

    def plot_cell_polys(self, ax=None):
        from matplotlib.pyplot import Figure
        if ax is None:
            f, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        for poly in self.cell_polys:
             xs, ys = poly.exterior.xy    
             ax.fill(xs, ys, fc='g', ec='none')





class Spectrometer():

    def __init__(self, case, chords=None, width=None, omega=1e-6, flip=True):
#, displ=0, width = 0.017226946, norm_zmag=True, crm=None):
        ''' Creates polygons for the chords from a data file '''

        # TODO: Add capability of specifiyng coords according to arrays
        from .GridData import Grid
        from numpy import ndarray
        uevars = [
            "ne",
            "te",
            "ni",
            "prad",
            "ng",
            "pradcff",
            "pradhyd",
        ]

        if ((width is None) or (width is False)) and \
            ((omega is None) or (omega is False)):
            self.line = True
        else:
            self.line = False


        self.case = case
        self.grid = Grid(case, flip=flip, variables=uevars)
        self.chords = []
        self.width = width
        self.omega = omega
        if isinstance(chords, str):
            self.read_chordfile(chords)
        elif isinstance(chords, (ndarray, list)):
            self.read_chordarray(chords)
            
    def read_chordarray(self, chords):
        for chord in chords:
            self.add_chord(chord, self.width)

    def read_chordfile(self, file):
        with open(file, 'r') as f:
            for line in f:
                if line.strip()[0] == "#":
                    continue
                parsed = line.strip().split('#')[0].split()
                [xs, ys, xe, ye] = [float(z) for z in parsed]
                self.add_chord(((xs, ys), (xe, ye)), self.width)
       
    def add_chord(self, points, width=None, omega=None):
        from shapely.geometry import Point
        if (width is None) and (omega is None):
            width = self.width
            omega = self.omega

        if (len(points[0]) == 3) and (self.line is False):
            print("Tangential chords treated as pencil-beams!")
            self.line = True

        pointcoords = []
        for point in points:
            pointcoords.append(Point(tuple(point)))

        if self.line:
            Chord = ChordLine
        else:
            Chord = ChordPoly

        self.chords.append(
            Chord(
                *pointcoords,
                self.grid,
                width = width,
                omega = omega
            )
        )

    def add_rates(self, path, species, ratetype, **kwargs):
        from uetools.UePostproc.ADASclass import ADASSpecies
        self.rates = ADASSpecies(path, species, ratetype, **kwargs)

    def plot_spectrometer(self, ax=None, **kwargs):
        ''' Plots a polygrid of Patches '''
        from matplotlib.pyplot import subplots, Figure
        if ax is None:
            f, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        for chord in self.chords:
            chord.plot_chord(ax, **kwargs)
        return ax.get_figure()

    def plot_setup(self):
        f = self.grid.plot_grid()
        self.plot_spectrometer(f)


    def get_nph(self, lam, chargestate, rates=None, 
        rtype = ['excit', 'recom', 'chexc']):
        if rates is None:   
            rates = self.rates
        species = rates.species.lower()
        return  rates.calc_emission(
                    self.grid.densities['e'][0], 
                    self.grid.vars['te'], 
                    self.grid.densities[species][chargestate], 
                    self.grid.densities[species][chargestate+1], 
                    self.grid.densities['h'][0],
                    chargestate, 
                    lam=lam, 
                    rtype=rtype
                )[lam]
 
    def calc_chord_emission(self, chargestate, rates=None, lam=None,
        rtype = ['excit', 'recom', 'chexc']):
        from numpy import array, zeros
        # Reset emission dictionary to avoid double-counting
        self.nph = {}
        self.emission = {}
        if rates is None:   
            rates = self.rates
        species = rates.species.lower()
        if lam is None:
            lam = list(rates.lines[chargestate].keys())
        elif isinstance(lam, (int, float)):
            if lam not in rates.lines[chargestate]:
                lam = rates.get_closest_line(lam, rates.linelist[chargestate])
            lam = [lam]
        elif isinstance(lam, (list, tuple)):
            _lam = []
            for l in lam:
                if l not in rates.lines[chargestate]:
                    l = rates.get_closest_line(l, rates.linelist[chargestate])
                _lam.append(l)
            lam = _lam
        for l in lam:
            self.nph[l] = self.get_nph(l, chargestate, rates=rates, rtype=rtype)
            emission_chord = []
            for chord in self.chords:
                emission_chord.append(chord.integrate_field(self.nph[l]))
            self.emission[l] = emission_chord
                
    def plot_chord_integral(self, field, ax=None, linestyle='-', marker='o',
        color='k', x = None):
        from matplotlib.pyplot import subplots, Figure
        if ax is None:
            f, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        if x is None:
            x = range(1, len(self.chords)+1)
        # TODO: Figure out what to use as X-axis
        integrals = []
        for chord in self.chords:
            integrals.append(chord.integrate_field(field))

        ax.plot(x, integrals, linestyle=linestyle, marker=marker, color=color)
        
        return ax.get_figure()


    def plot_chord_spectra(self, chord, chargestate, ax=None, linestyle='-',
        color='k', rates=None, **kwargs):
        from matplotlib.pyplot import subplots, Figure
        if ax is None:
            f, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        self.calc_chord_emission(chargestate)
        for line, chords in self.emission.items():
            ax.semilogy([line, line], [0, chords[chord]], linestyle=linestyle, color=color)
        return ax.get_figure()
        

    def plot_chord_emission(self, lam, chargestate, ax=None, linestyle='-',
        marker='o', color='k', rates=None, x=None, **kwargs):
        from matplotlib.pyplot import subplots, Figure
        if ax is None:
            f, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        if not isinstance(lam, (float, int)):
            raise Exception("Please select a single spectroscopic line to plot")
        self.calc_chord_emission(chargestate, lam=lam, rates=rates, **kwargs)
        if x is None:
            x = range(1, len(self.chords)+1)
        # TODO: Figure out what to use as X-axis
        ax.plot(x, self.emission[lam], linestyle=linestyle, marker=marker, color=color)
        
        return ax.get_figure()

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


