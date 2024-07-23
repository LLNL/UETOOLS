class ChordLine():
    """ Creates a single 0D chord object """
    def __init__(self, p0, p1, grid, *args, res = 500, **kwargs):
        from shapely.geometry import Point, LineString
        from numpy import array, pi, linspace
        from scipy.interpolate import interp1d

        self.p0 = p0 # Chord starting point
        self.p1 = p1 # Chord end point
        if abs(self.p0.z - self.p1.z) < 1e-5:
            raise Exception("Tangential lines in the horizontal plane not yet implemented")


        # Get the total vector length 
        self.length = p0.distance(p1)
        
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
    
        self.chunkL = self.length/res
        points = [generator(0)]
        segments = {
            'outside': [],
            'inside': [],
        }
        for chunk in linspace(0, self.length, res):
            points.append(generator(chunk))
    
        self.chord = LineString(points)

        p0rz = Point(points[0])
        p1rz = Point(points[-1])
        # Find the cells intersected by the chord and store 
        # their path-lengths and indices
        self.intersects = []
        self.cells = {}
        for ix, row in grid.map.items():
            for iy, cell in row.items():
                if cell.polygon.contains(p0rz):
                    raise ValueError("Chord start point inside grid")
                elif cell.polygon.contains(p1rz):
                    raise ValueError("Chord end point inside grid")
                buff = []
                if grid.map[ix][iy].polygon.intersects(self.chord):
                    for intersect in grid.map[ix][iy].polygon.exterior.intersection(self.chord).geoms:
                        buff.append(self.fL(intersect.y))
                    if ix not in self.cells:
                        self.cells[ix] = {}
                    if iy not in self.cells[ix]:
                        self.cells[ix][iy] = {
                            'dL': 0
                        }
                    buff.sort()
                    self.cells[ix][iy]['dL'] += buff[1] - buff[0]
                    self.intersects.append(buff)


    def plot_chord(self, ax=None, color='r', alpha=0.2):
        from matplotlib.pyplot import Figure, subplots
        if ax is None:
            f, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        ax.plot(*self.chord.xy, '-', color=color)
        return ax.get_figure()


        
    def integrate_field(self, field):
        integral = 0
        for ix, row in self.cells.items():
            for iy, L in row.items():
                integral += field[ix,iy]*L['dL']
        return integral



        

class ChordPoly():
    ''' Creates a single chord object '''
    def __init__(self, p0, p1, grid, width=None, omega=None):
        from shapely.geometry import Point
        from numpy import array, pi
        self.p0 = p0 # Chord starting point
        self.p1 = p1 # Chord end point

        self.length = self.p0.distance(self.p1)

        if (width is None) and (omega is None):
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





    # TODO: Implement using integrate_field


    def add_emission(self, emission, dL):
        """ Populates self.emission with emissivities from volumetric emission
        
            Arguments:
            emission - dictionary with lines and corresponding volum. emission 
            in 1/(s cm**-3)
            dL - path length through cell in m
            
            Populates self.emission with emissivities in units
            1/(s * sr * cm**-2)


        """
        from numpy import pi
        for lam, e in emission.items():
            try:
                self.emission[lam] += e*dL
            except:
                self.emission[lam] = e*dL
        for lam, emission in self.emission.items():
            emission *= 1/(4*pi*1e-2)

    def calc_emission(self, rates, chargestate, lam=None, rerun=True,
        rtype = ['excit', 'recom', 'chexc']):
        from numpy import array
        species = rates.species.lower()
        # Reset emission dictionary to avoid double-counting
        self.emission = {}
        if isinstance(lam, (int, float)):
            if lam not in rates.lines[chargestate]:
                lam = rates.get_closest_line(lam, rates.linelist[chargestate])
        elif isinstance(lam, (list, tuple)):
            _lam = []
            for l in lam:
                if l not in rates.lines[chargestate]:
                    l = rates.get_closest_line(l, rates.linelist[chargestate])
                _lam.append(l)
            lam = _lam
        for _, obj in self.dL.items():
            species = rates.species.lower()
            o = obj['cell']
            obj['cell'].emission[species] = rates.calc_emission(
                    obj['cell'].ne, 
                    obj['cell'].te, 
                    obj['cell'].densities[species][chargestate], 
                    obj['cell'].densities[species][chargestate+1], 
                    obj['cell'].densities['h'][0],
                    chargestate, 
                    lam=lam, 
                    rtype=rtype
            )
            self.add_emission(obj['cell'].emission[species], obj['dL'])
        
             
    def plot_emission(self, rates, chargestate,
        rtype = ['excit', 'recom', 'chexc']):
        from matplotlib.pyplot import subplots
        species = rates.species.lower()
        # Reset emission dictionary to avoid double-counting
        self.emission = {}
        for _, obj in self.dL.items():
            species = rates.species.lower()
            o = obj['cell']
            obj['cell'].emission[species] = rates.calc_emission(
                    obj['cell'].ne, 
                    obj['cell'].te, 
                    obj['cell'].densities[species][chargestate], 
                    obj['cell'].densities[species][chargestate+1], 
                    obj['cell'].densities['h'][0],
                    chargestate, 
                    lam=None, 
                    rtype=rtype
            )
            self.add_emission(obj['cell'].emission[species], obj['dL'])


        f, ax = subplots()
        for x, y in self.emission.items():
            ax.semilogy([x,x], [0,y], 'r-')
        return f











    def calc_L(self, point1, point2):
        return ( (point1.x - point2.x)**2 + (point1.y - point2.y)**2 )**0.5

    def add_contained_cell(self, cell):
        from numpy import sqrt, pi
        # Add the entire cell polygon to the list of grid cells
        # inside the LOS polygon:

        # Calculate dL for the integration by postulating a
        # rectangular cell orthogonal to the line-of-sight with
        # area of the original cell and width determined by the
        # distance of the cell from the origin of the LOS and
        # the width of the LOS at its end:
        L_cell = self.calc_L(cell.polygon.centroid, self.p0)
        w_orthog = L_cell/self.length*self.width
        dL = cell.polygon.area/w_orthog
        
        self.dL[str(cell.indices)] = {
            'cell': cell,
            'dL': dL,
            'L': L_cell,
            'A': cell.polygon.area,
        }
        self.cell_polys.append(cell.polygon)

    def add_intersected_cell(self, cell):
        from numpy import sqrt, pi
        # Determine the part of the grid cell polygon that
        # intersects with the line-of-sight polygon:
        clipped = self.poly.intersection(cell.polygon)

        # Calculate dL for the integration by postulating a
        # rectangular cell orthogonal to the line-of-sight with
        # area of the original cell and width determined by the
        # distance of the cell from the origin of the LOS and
        # the width of the LOS at its end: 
        L_cell = self.calc_L(clipped.centroid, self.p0)
        w_orthog = L_cell/self.length*self.width
        dL = clipped.area/w_orthog
            
        self.dL[str(cell.indices)] = {
            'cell': cell,
            'dL': dL,
            'L': L_cell,
            'A': cell.polygon.area,
        }
        self.cell_polys.append(clipped)


    def calc_LOS_spectra(self,lower, upper):
        from numpy import nonzero
        ret = [0,0]
        for i in nonzero(self.LOS_H2_emission_integral)[0]:
            ene = self.LOS_H2_energies[i]
            if 10*ene > lower and 10*ene < upper:
                ret[0] += self.LOS_H2_emission_integral[i]
        for i in nonzero(self.LOS_H_emission_integral)[0]:
            ene = self.LOS_H_energies[i]
            if 10*ene > lower and 10*ene < upper:
                ret[1] += self.LOS_H_emission_integral[i]
        return ret     
            
    def plot_LOS_spectra(self, ax=None, xlim=(500,8500), ylim=(None, None),yaxis='lin',yscale=1):
        from matplotlib.pyplot import subplots
        from numpy import nonzero
        ret = False
        
        if ax is None:
            f, ax = subplots()
            ret = True

        if yaxis == 'lin':
            pl = ax.plot
        elif yaxis == 'log':
            pl = ax.semilogy
        
    
        for i in nonzero(self.LOS_H2_emission_integral)[0]:
            ene = 10*self.LOS_H2_energies[i]
            if (ene > xlim[0]) and (ene < xlim[1]):
                pl([ene, ene], [0, yscale*self.LOS_H2_emission_integral[i]], 'r-')

        for i in nonzero(self.LOS_H_emission_integral)[0]:
            ene = 10*self.LOS_H_energies[i]
            if (ene > xlim[0]) and (ene < xlim[1]):
                pl([ene, ene], [0, yscale*self.LOS_H_emission_integral[i]], 'b-')
        

        ax.set_xlabel(r'Wavelength [Å]')
        ax.set_ylabel(r'Intensity [ph/sr/s/$\rm cm^3$]')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        if ret is True:
            return f

    def LOS_integral(self, grid, reevaluate=False, nH=None):
        from numpy import sqrt, sum, array
        from tqdm import tqdm
        self.LOS_cells = [] 
        self.LOS_area = []
        self.LOS_H2_emission = []
        self.LOS_H_emission = []
        self.LOS_dL = []
        self.LOS_L = []
        
        
        # Loop through each inversion grid cell:
        for cell in tqdm(grid.cells):
            
            # Check if the grid cell is completely inside the line-of-
            # sight polygon:
            if self.poly.contains(cell.polygon):
                
                # Add the entire cell polygon to the list of grid cells
                # inside the LOS polygon:
                self.LOS_cells.append(cell)
                
                # Calculate dL for the integration by postulating a
                # rectangular cell orthogonal to the line-of-sight with
                # area of the original cell and width determined by the
                # distance of the cell from the origin of the LOS and
                # the width of the LOS at its end:
                L_cell = sqrt((cell.polygon.centroid.x - self.points[0].x)**2 + 
                                (cell.polygon.centroid.y - self.points[0].y)**2)
                w_orthog = L_cell/sqrt((self.points[1].x - self.points[0].x)**2 + 
                                (self.points[1].y - self.points[0].y)**2)*self.width
                dL = cell.polygon.area/w_orthog
                
                # Store the emission in the grid cell multiplied by dL
                # in the list of emission inside the LOS polygon:
                # TODO: do CRM if it is not done
                if (cell.crmeval is False) or (reevaluate is True):
                    cell.calc_emissivity(self.crm, nH=nH)
                self.LOS_H_emission.append(cell.H_emission[1]*dL)
                self.LOS_H2_emission.append(cell.H2_emission[1]*dL)
                self.LOS_L.append(L_cell)
                self.LOS_dL.append(dL)
                
                self.LOS_area.append(cell.polygon.area)
                self.LOS_H_energies = cell.H_emission[0]
                self.LOS_H2_energies = cell.H2_emission[0]
            
            # Alternatively, check if part of the grid cell is inside
            # the line-of-sight polygon:
            elif self.poly.intersects(cell.polygon):
                
                # Determine the part of the grid cell polygon that
                # intersects with the line-of-sight polygon:
                grid_poly_clipped = self.poly.intersection(cell.polygon)
                
                # Add the intersecting part of the cell polygon to the
                # list of grid cells inside the LOS polygon:
                # TODO: How to do this using objects?
                self.LOS_cells.append(cell)
                
                # Calculate dL for the integration by postulating a
                # rectangular cell orthogonal to the line-of-sight with
                # area of the original cell and width determined by the
                # distance of the cell from the origin of the LOS and
                # the width of the LOS at its end: 
                L_cell = sqrt((grid_poly_clipped.centroid.x - self.points[0].x)**2 + 
                                (grid_poly_clipped.centroid.y - self.points[0].y)**2)
                w_orthog = L_cell/sqrt((self.points[1].x - self.points[0].x)**2 + 
                                (self.points[1].y - self.points[0].y)**2)*self.width
                dL = grid_poly_clipped.area/w_orthog
                
                # Store the emission in the grid cell multiplied by dL
                # in the list of emission inside the LOS polygon:
                if (cell.crmeval is False) or (reevaluate is True):
                    cell.calc_emissivity(self.crm, nH=nH)
                self.LOS_H_emission.append(cell.H_emission[1]*dL)
                self.LOS_H2_emission.append(cell.H2_emission[1]*dL)
                self.LOS_dL.append(dL)
                self.LOS_L.append(L_cell)
                
                self.LOS_area.append(cell.polygon.area)
                self.LOS_H_energies = cell.H_emission[0]
                self.LOS_H2_energies = cell.H2_emission[0]
        
        self.LOS_cells = array(self.LOS_cells)
        self.LOS_area = array(self.LOS_area)
        self.LOS_H_emission = array(self.LOS_H_emission)
        self.LOS_H2_emission = array(self.LOS_H2_emission)
        self.LOS_dL = array(self.LOS_dL)
        self.LOS_L = array(self.LOS_L)
        # Calculate the LOS-integrated emission as a sum of each
        # emission*dL element stored for the current LOS:   
        self.LOS_H_emission_integral = sum(self.LOS_H_emission, axis=0)
        self.LOS_H2_emission_integral = sum(self.LOS_H2_emission, axis=0)
        


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
            print("Tangential chords treated as bencil-beams!")
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

