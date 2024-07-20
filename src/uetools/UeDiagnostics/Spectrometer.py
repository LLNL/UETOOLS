
class Chord():
    ''' Creates a single chord object '''
    def __init__(self, p0, p1, grid, width=None, omega=None):
        from shapely.geometry import Point
        from numpy import pi
        self.p0 = p0 # Chord starting point
        self.p1 = p1 # Chord end point

        self.length = ( (p0.x - p1.x)**2 + (p0.y - p1.y)**2)**0.5

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
        self.contained = {}
        self.intersected = {}
        self.dL = {}
        self.cell_polys = []
        self.emission = {}
        # TODO: Implement check on end-point in grid cell!
        for cell in grid.cells:
            if cell.polygon.contains(self.p0):
                raise ValueError("Chord start point inside grid")
            elif cell.polygon.contains(self.p1):
                raise ValueError("Chord end point inside grid")
            elif self.poly.contains(cell.polygon):
                self.add_contained_cell(cell)
            elif self.poly.intersects(cell.polygon):
                self.add_intersected_cell(cell)


    def calc_width(self, omega):
        ''' Calculates the 2D 2*r at the end-point for a given solid angle '''
        # Calculates 
        from numpy import pi
        return self.length*(2/(2*pi-omega))*(4*pi*omega-omega**2)**0.5

    def calc_solidangle(self, w):
        ''' Calculates the solid angle for a give 2D width 2*r at the end-point '''
        from numpy import pi
        return 2*pi*(1-1/(1+w**2/4/self.length**2)**0.5)

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

    def integrate(self, variable=None):
        ret = 0
        if variable is None:
            for key, obj in self.dL.items():
                ret += obj['dL']
        else:
            for key, obj in self.dL.items():
                ret += obj['dL']*obj['dL'].__getattribute__(variable)
        return ret






    def plot_chord(self, ax=None, poly=True, line=True, color='r', alpha=0.2):
        from matplotlib.pyplot import Figure
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

    def __init__(self, case, chordfile=None, width=None, omega=1e-6, flip=True):
#, displ=0, width = 0.017226946, norm_zmag=True, crm=None):
        ''' Creates polygons for the chords from a data file '''

        # TODO: Add capability of specifiyng coords according to arrays
        from .GridData import Grid
        uevars = [
            "ne",
            "te",
            "ni",
            "prad",
            "ng",
            "pradcff",
            "pradhyd",
        ]

        self.case = case
        self.grid = Grid(case, flip=flip, variables=uevars)
        self.chords = []
        self.width = width
        self.omega = omega
        if chordfile is not None:
            self.read_chordfile(chordfile)


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

        p0 = Point((points[0][0], points[0][1]))
        p1 = Point((points[1][0], points[1][1]))

        self.chords.append(
            Chord(
                p0,
                p1,
                self.grid,
                width = width,
                omega = omega
            )
        )

    def add_rates(self, path, species, ratetype, **kwargs):
        from uetools.UePostproc.ADASclass import ADASSpecies
        self.rates = ADASSpecies(path, species, ratetype, **kwargs)

    
    def calc_chord_emission(self, lam, chargestate, rerun=True,
        rtype = ['excit', 'recom', 'chexc'], rates=None, **kwargs):
        if rates is None:
            rates = self.rates
        for chord in self.chords:
            chord.calc_emission(rates, chargestate, lam, rerun, rtype)
        return [x.emission[lam] for x in self.chords]


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


    def plot_chord_intensity(self, lam, chargestate, ax=None, linestyle='-',
        marker='o', color='k', rates=None, **kwargs):
        from matplotlib.pyplot import subplots, Figure
        if ax is None:
            f, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        emission = self.calc_chord_emission(lam, chargestate, rates=rates, **kwargs)

        # TODO: Figure out what to use as X-axis
        ax.plot(range(1, len(self.chords)+1),
            emission,
            linestyle=linestyle, marker=marker, color=color)
        
        return ax.get_figure()












    def calculate_LOS_integral(self, grid, **kwargs):
        ''' Calculates the LOS integrals from a Grid object '''
        for chord in self.chords:
            if chord.integrated is False:
                chord.LOS_integral(grid, **kwargs)
    
    def calculate_LOS_integral_data(self, data, grid, **kwargs):
        ret = []
        for chord in self.chords:
            chord.LOS_integral_value(data, grid, **kwargs)
            ret.append(chord.LOS_value_integral)

        return ret

    def plot_chord_angle(self, grid, interval = (5980,6330), mol=True, ax=None, style='ko-',yscale=1,dtheta=0, **kwargs):
        from matplotlib.pyplot import subplots
        from numpy import array, mean, array#, arccos, dot, degrees
        from math import degrees, atan2
        from numpy.linalg import norm
        from shapely.geometry import Point
        if ax is None:
            f, ax = subplots()
        self.calculate_LOS_integral(grid)
        chord_ends = []
        detector = []
        angles = []
        intensities = []
        for chord in self.chords:
            detector.append(array(chord.points[0].xy))
            chord_ends.append(array(chord.points[1].xy))
            # Integrate chord intensity in interval
            intensities.append(chord.calc_LOS_spectra(interval[0], interval[1]))

        # Mean upper point
        detector = mean(array(detector), axis=0)
        for point in chord_ends:
#        for i in range(len(chord_ends)):
#            angles.append(degrees(atan2(chord_ends[i][1]-detector[i][1], chord_ends[i][0]-detector[i][0]) - 
#                    atan2(detector[i][1]-detector[i][1], detector[i][0]+1-detector[0][0])))
            angles.append(degrees(atan2(point[1]-detector[1], point[0]-detector[0]) - 
                    atan2(self.osep[1]-detector[1], self.osep[0]-detector[0])))
#            chord = point - detector
#            sep = self.osep - detector
#            angles.append(degrees(arccos(dot(chord[:,0], sep[:,0]) / (norm(chord) * norm(sep)))))
        # Calculate angle between line and sep
        ax.plot(array(angles)+dtheta, yscale*array(intensities)[:,0**mol],style,**kwargs)
        return ax.get_figure()
                
 

    def plot_chord_angle_data(self,data,  grid, ax=None, style='ko-',yscale=1, dtheta=0, **kwargs):
        from matplotlib.pyplot import subplots
        from numpy import array, mean, array#, arccos, dot, degrees
        from math import degrees, atan2
        from numpy.linalg import norm
        from shapely.geometry import Point
        if ax is None:
            f, ax = subplots()
        chord_ends = []
        detector = []
        angles = []
        for chord in self.chords:
            detector.append(array(chord.points[0].xy))
            chord_ends.append(array(chord.points[1].xy))

        # Mean upper point
        detector = mean(array(detector), axis=0)
#        for i in range(len(chord_ends)):
        for point in chord_ends:
#            angles.append(degrees(atan2(chord_ends[i][1]-detector[i][1], chord_ends[i][0]-1-detector[i][0]) - 
#                    atan2(detector[i][1]-detector[i][1], detector[i][0]+1-detector[0][0])))
            angles.append(degrees(atan2(point[1]-detector[1], point[0]-detector[0]) - 
                    atan2(self.osep[1]-detector[1], self.osep[0]-detector[0])))
#            chord = point - detector
#            sep = self.osep - detector
#            angles.append(degrees(arccos(dot(chord[:,0], sep[:,0]) / (norm(chord) * norm(sep)))))
        # Calculate angle between line and sep
        ax.plot(array(angles)+dtheta, yscale*array(data),style,**kwargs)
        return ax.get_figure()
                
 


