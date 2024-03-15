
class Chord():
    ''' Creates a single chord object '''
    def __init__(self, points, width, grid):
        from shapely.geometry import Point
        self.points = (Point((points[0][0], points[0][1])), 
                        Point((points[1][0], points[1][1])))
        self.width = width
        self.create_poly()
        self.contained = {}
        self.intersected = {}
        self.dL = {}
        self.cell_polys = []
        for cell in grid.cells:
            if self.poly.contains(cell.polygon):
                self.add_contained_cell(cell)
            elif self.poly.intersects(cell.polygon):
                self.add_intersected_cell(cell)

    def add_contained_cell(self, cell):
        from numpy import sqrt, pi
        # Add the entire cell polygon to the list of grid cells
        # inside the LOS polygon:

        # Calculate dL for the integration by postulating a
        # rectangular cell orthogonal to the line-of-sight with
        # area of the original cell and width determined by the
        # distance of the cell from the origin of the LOS and
        # the width of the LOS at its end:
        L_cell = (
                    (cell.polygon.centroid.x - self.points[0].x)**2 + \
                    (cell.polygon.centroid.y - self.points[0].y)**2
                )**0.5
        w_orthog = L_cell/(
                        (self.points[1].x - self.points[0].x)**2 + \
                        (self.points[1].y - self.points[0].y)**2
                )**0.5*self.width
        dL = cell.polygon.area/w_orthog
        
#        self.contained[str(cell.indices)] = {
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
        L_cell = (
                    (clipped.centroid.x - self.points[0].x)**2 + \
                    (clipped.centroid.y - self.points[0].y)**2
                )**0.5
        w_orthog = L_cell/(
                        (self.points[1].x - self.points[0].x)**2 + \
                        (self.points[1].y - self.points[0].y)**2
                )**0.5*self.width
        dL = clipped.area/w_orthog
            
#        self.intersected[str(cell.indices)] = {
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

        L = sqrt((self.points[1].x - self.points[0].x)**2 + 
                        (self.points[1].y - self.points[0].y)**2)
        theta = arctan((self.points[1].x - self.points[0].x)/(self.points[1].y - 
                        self.points[0].y))
        
        # Elongate the LOS arbitrarily to make sure that it crosses the
        # wall boundary and calculate the new end coordinates, length
        # and end width of the LOS:
        r_elong = self.points[1].x - 1.0 * sin(theta)
        z_elong = self.points[1].y - 1.0 * cos(theta)
        L_elong = sqrt(
                    (r_elong - self.points[0].x)**2 + \
                    (z_elong - self.points[0].y)**2
        )
        w_elong = L_elong/L*self.width
        
        # Add the elongated line-of-sight to the list of LOS polygons: 
        self.poly = Polygon([(self.points[0].x, self.points[0].y), 
                                (r_elong - 0.5*w_elong*cos(theta), 
                                    z_elong + 0.5*w_elong*sin(theta)), 
                                (r_elong + 0.5*w_elong*cos(theta), 
                                    z_elong - 0.5*w_elong*sin(theta))])

    def integrate(self, variable=None):
        if variable is None:
            ret = 0
            for key, obj in self.dL.items():
                ret += obj['dL']
        return ret






    def plot_chord(self, ax=None, poly=True, line=True, color='r', alpha=0.2):
        from matplotlib.pyplot import Figure
        if ax is None:
            f, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        if line is True:
            (p1, p2) = self.points
            ax.plot([p1.x, p2.x], [p1.y, p2.y], '-', linewidth=0.5,color=color)
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

    def __init__(self, case, fname=None, width=0.001, flip=True):
#, displ=0, width = 0.017226946, norm_zmag=True, crm=None):
        ''' Creates polygons for the chords from a data file '''
        from .GridData import Grid
        self.case = case
        self.grid = Grid(case, flip=flip)
        self.chords = []
        self.width = width
        if fname is not None:
            self.read_chordfile(fname)


    def read_chordfile(self, file):
        with open(file, 'r') as f:
            for line in f:
                try:
                    parsed = line.strip().split('#')[0].split()
                    [xs, ys, xe, ye] = [float(z) for z in parsed]
                    self.add_chord(((xs, ys), (xe, ye)), self.width)
                except:
                    print(f"Could not read line: '{line}'")

    def add_chord(self, points, width):
        self.chords.append(
            Chord(
                points,
                width,
                self.grid
            )
        )








    def plot_spectrometer(self, ax=None, **kwargs):
        ''' Plots a polygrid of Patches '''
        from matplotlib.pyplot import subplots, Figure
        if ax is None:
            f, ax = subplots()
            ret = True
        else:
            ret = False

        for chord in self.chords:
            chord.plot_chord(ax, **kwargs)

        if ret is True:
            return f

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
                
 


