class ChordLine:
    """Creates a single pencil-beam chord object

    Atrributes
    ----------
    p0: shapely Point object of start point in (x,y,z) space
    p1: shapely Point object of end point in (x,y,z) space
    length: total length of chord in (x,y,z) space
    fx: parametrized function f(L) = x
    fy: paremetrized function f(L) = y
    fz: parametrized function f(L) = z
    fL: paramtrized function f(z) = L
    chord: shapely LineString object of chord in (r,z) space
    p0rz: shapely Point object of start point in (r,z) space
    p1rz: shapely Point object of end point in (r,z) space
    cells: nested dict containing optical path-length through cell
    in (x,y,z) space, cells[ix][iy] = dL

    Methods
    -------
    plot_chord(ax=None, color='r', linewidth=0.5)
        plots the (r,z) representation of the chord
    integrate_field(field)
        returns the LOS integral through field in (x,y,z) space
    """

    def __init__(self, p0, p1, grid, *args, res=500, **kwargs):
        """Creates a generalized pe ncil-beam

        Pencil-beams are used for 1D line-integration and can accommodate
        tangential view-chords.

        Arguments
        ---------
        p0 - start point of chord, Shapely Point object
        p1 - end point of chord, Shapely Point object
        grid - UETOOLS Grid object used to link chord to geometry

        Keyword arguments
        -----------------
        res : int (default = 500)
            number of segments to split the chord into

        Returns
        -------
        None
        """
        from shapely.geometry import Point, LineString
        from numpy import array, pi, linspace
        from scipy.interpolate import interp1d

        self.p0 = p0  # Chord starting point
        self.p1 = p1  # Chord end point
        # Catch edge-case of purely horizontal tangential view-chord
        if abs(self.p0.z - self.p1.z) < 1e-5:
            raise Exception(
                "Tangential lines in the horizontal plane not yet implemented"
            )
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
            return ((self.fx(L) ** 2 + self.fy(L) ** 2) ** 0.5, self.fz(L))

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
                    for intersect in (
                        grid.map[ix][iy].polygon.exterior.intersection(self.chord).geoms
                    ):
                        # Use the polidal Z-coordinate for back-projection
                        # NOTE: this back-projection is what limits the routine
                        # to non-horizontal tangential chords only
                        buff.append(self.fL(intersect.y))
                    # Set up objects if this is the first chord intersection
                    # with the cell
                    if ix not in self.cells:
                        self.cells[ix] = {}
                    if iy not in self.cells[ix]:
                        self.cells[ix][iy] = {"dL": 0}
                    # Store the path-length through the cell
                    buff.sort()
                    self.cells[ix][iy]["dL"] += buff[1] - buff[0]

    def plot_chord(self, ax=None, color="r", linewidth=0.5):
        """Plots poloidal projection of the chord object

        Keyword arguments
        -----------------
        ax : matplotlib.Axes or matplotlib.Figure (default = None)
            axis to plot chord on. If none, creates a new figure
        color : str (default = 'r')
            color of chord
        linewidth : float (default = 0.5)
            width of line representing chord

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
        ax.plot(*self.chord.xy, "-", color=color, linewidth=linewidth)
        # Ensure Figure object is returned
        return ax.get_figure()

    def integrate_field(self, field):
        """Returns the path-integral of chord through field

        Arguments
        ---------
        field: 2D array of shape compatible with UEDGE grid containing
                field values to integrate through

        Returns
        -------
        Float: line-integral of chord through field in [field]*m
        """
        integral = 0
        # Loop through all cells intersected by chord
        for ix, row in self.cells.items():
            for iy, L in row.items():
                # Access the correct cell data from field and perform
                # numerical line-integral
                integral += field[ix, iy] * L["dL"]
        return integral


class ChordPoly:
    """Creates a 2D poloidal chord object

    Attributes
    ----------
    p0: shapely Point object of start point in (x,y,z) space
    p1: shapely Point object of end point in (x,y,z) space
    length: total length of chord in (x,y,z) space
    width: 2D fan width at end point in meters
    omega: solid angle (in 3D) of chord
    cells: nested dict containing optical path-length through cell
        in (x,y,z) space, cells[ix][iy]['dL'] = dL
    cell_polys: list of shapely polygons constiuting the chord throguh
        the plasma grid

    Methods
    -------
    calc_width(omega)
        Calculates 2D width at end-point based on given 3D solid angle
    calc_solidangle(w)
        Calculates 3D solid angle based on given 2D width at end-point
    get_dL(cell)
        Calculates path-length representation through the cell
    create_poly()
        Creates the 2D chord fan shapely Polygon
    integrate_field(field)
        returns the LOS integral through field in (x,y,z) space
    plot_chord(ax=None, poly=True, line=True, color='r', alpha=0.2)
        plots the chord
    plot_cell_polys(ax=None)
        plots the compounded chord from intersected/contained cell
        polygons
    """

    def __init__(self, p0, p1, grid, width=None, omega=None):
        """Creates a generalized pencil-beam

        Pencil-beams are used for 1D line-integration and can accommodate
        tangential view-chords. Representation considers a 2D fan
        representation of the chords for small solid angles. For large
        solid angles, the chords are strictly 3D objects and this
        representation will yield inaccurate results

        Arguments
        ---------
        p0 - shapely Point object of start point in (x,y,z) space
        p1 - shapely Point object of end point in (x,y,z) space
        grid - uetools.UeDiagnostics Grid object

        Keyword arguments
        -----------------
        width : float (default = None)
            Width of beam at end-point in m. Calculated if omega is set.
            Assumed to be 1e-9m if both width and omega are None
            (approximate pencil-beam).
        omega : floar (default = None)
            Solid-angle of detector at start-point. Calculated if width
            is set.

        """
        from shapely.geometry import Point
        from numpy import array, pi

        self.p0 = p0  # Chord starting point
        self.p1 = p1  # Chord end point
        # Get total chord length
        self.length = self.p0.distance(self.p1)
        # Set up the 2D polygon triangle approximating the chord
        # view geomtery
        if (width is None) and (omega is None):
            # Approximate pencil beam
            self.width = 1e-9  # Pencil beam
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
                        self.cells[ix][iy] = {"dL": 0}
                    self.cells[ix][iy]["dL"] += self.get_dL(cell.polygon)

    def calc_width(self, omega):
        """Calculates the 2D 2*r at the end-point for a given solid angle"""
        # Calculates
        from numpy import pi

        return (
            self.length * (2 / (2 * pi - omega)) * (4 * pi * omega - omega**2) ** 0.5
        )

    def calc_solidangle(self, w):
        """Calculates the solid angle for a give 2D width 2*r at the end-point"""
        from numpy import pi

        return 2 * pi * (1 - 1 / (1 + w**2 / 4 / self.length**2) ** 0.5)

    def get_dL(self, cell):
        """Adds the path-length of contained cell to ChordPoly.cells


        Calculate the line-integration path-length through the cell by
        postulating a rectangular cell orthogonal to the line-of-sight
        with area of the original cell, or union of the chord beam and
        cell if cell only intersected, and width determined by the
        distance of the cell from the origin of the LOS and the width of
        the LOS at its end. Calculations based on B. Lomanowski's
        PYPROC/PEST scripts.

        Arguments
        ---------
        cell - shapely Polygeon representation of the cell

        Returns
        -------
        float : path-length representation through cell
        """
        # If the cell is not contained, get the intersecting polygon
        if not self.poly.contains(cell):
            cell = self.poly.intersection(cell)
        L_cell = cell.centroid.distance(self.p0)
        w_orthog = L_cell / self.length * self.width
        return cell.area / w_orthog

    def create_poly(self):
        """Creates a shapely.Polygon object representing the chord fan"""
        from numpy import sqrt, cos, sin, arctan
        from shapely.geometry import Polygon

        denom = self.p1.y - self.p0.y
        # Avoid zero-divisors
        if denom == 0:
            denom = 1e-10
        theta = arctan((self.p1.x - self.p0.x) / denom)

        # Add the elongated line-of-sight to the list of LOS polygons:
        self.poly = Polygon(
            [
                (self.p0.x, self.p0.y),
                (
                    self.p1.x - 0.5 * self.width * cos(theta),
                    self.p1.y + 0.5 * self.width * sin(theta),
                ),
                (
                    self.p1.x + 0.5 * self.width * cos(theta),
                    self.p1.y - 0.5 * self.width * sin(theta),
                ),
            ]
        )

    def integrate_field(self, field):
        """Returns the path-integral of chord through field

        Arguments
        ---------
        field: 2D array of shape compatible with UEDGE grid containing
                field values to integrate through

        Returns
        -------
        Float: line-integral of chord through field in [field]*m
        """
        integral = 0
        for ix, row in self.cells.items():
            for iy, L in row.items():
                integral += field[ix, iy] * L["dL"]
        return integral

    def plot_chord(self, ax=None, poly=True, line=True, color="r", alpha=0.2):
        """Plots poloidal projection of the chord object

        Keyword arguments
        -----------------
        ax : matplotlib.Axes or matplotlib.Figure (default = None)
            axis to plot chord on. If none, creates a new figure
        poly : bool (default = True)
            draws the polygon representation of the chord if True
        line : bool (default = True)
            draws the center-line of the polygon if True
        color : str (default = 'r')
            color of chord
        alpha : float (default = 0.2)
            alpha of the chord polygon representation

        Returns:
        ========
        Figure object
        """
        from matplotlib.pyplot import Figure, subplots

        if ax is None:
            f, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        if line is True:
            ax.plot(
                [self.p0.x, self.p1.x],
                [self.p0.y, self.p1.y],
                "-",
                linewidth=0.5,
                color=color,
            )
        if poly is True:
            xs, ys = self.poly.exterior.xy
            ax.fill(xs, ys, alpha=alpha, fc=color, ec="none")
        return ax.get_figure()

    def plot_cell_polys(self, ax=None):
        """Plots compound chord from cell intersections on ax

        Keyword arguments
        -----------------
        ax : matplotlib.Axes or matplotlib.Figure (default = None)
            axis to plot chord on. If none, creates a new figure

        Returns:
        ========
        Figure object
        """
        from matplotlib.pyplot import Figure

        if ax is None:
            f, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        for poly in self.cell_polys:
            xs, ys = poly.exterior.xy
            ax.fill(xs, ys, fc="g", ec="none")
        return ax.get_figure()


class Spectrometer:
    """Creates a spectrometer Object

    Attributes
    ----------
    case: uetools.Case object linking spectrometer to the Case
    grid: uetools.UeDiagnostics.Grid object for resolving geometry
    chords: list of chords in the spectrometer
    width: width of chords in spectrometer
    omega: solid angle of chords in spectrometer

    Created attributes
    ------------------
    rates: uetools.UeUtils.ADASSpecies object containing rate data

    Methods
    -------
    read_chordarray(chords)
        Creates chords from a nested list start-end point coord pairs
    read_chordfile(file)
        Reads nested list start-end point coordinate pairs from file
    add_chord(points, width=None, omega=None)
        Adds a chord to the list from a len-2 list of r,z coord pairs
    add_rates(path, species, ratetype, **kwargs)
        Adds a ADASSpecies objects to Spectrometer.rates
    plot_spectrometer(ax=None, **kwargs)
        Plots the spectrometer chords on ax
    plot_setup()
        Plots the spectrometer and setup on a new figure
    get_nph(lam, chargestate, rates=None,
            rtype=['excit', 'recom', 'chexc'])
        Gets the photon density on the UEDGE grid
    calc_chord_emission(chargestate, rates=None, lam=None
            rtype=['excit', 'recom', 'chexc'])
        Calculates the emission for each chord
    plot_chord_integral(field, ax=None, linestyle='-', marker='o',
            color='k', x=None)
        Plots the line-integrated values of field along each chord
    plot_chord_spectra(chord, chargestate, ax=None, linestyle='-',
            color='k', rates=None, **kwargs)
        Plots the spectral line distribution integrated along chord
    plot_chord_emission(lam, chargestate, ax=None, linestyle='-',
            marker='o', color='k', rates=None, x=None, **kwargs)
        Plots the emission along each chord


    """

    def __init__(self, case, chords=None, width=None, omega=1e-6, flip=True):
        """Creates a Spectrometer object

        Arguments
        ---------
        case - uetools.Case obect linking Spectrometer to a case

        Keyword arguments
        -----------------
        chords : list or str (default = None)
            Nested list of [[r_start, z_start], [r_end, z_end]] points
            to create data points or path to file containing the chord
            start and end points (see read_chordfile for details). If
            None, initializes Spectrometer without chords
        width : float (default = None)
            Width of chord at end-points. Calculated from omega if None.
            Assumes pencil-beams if both omega and width are None. If
            both width and omega None, uses pencil-beam approximation.
        omega : float (default = 1.e-6)
            Solid angle of chord start point. Results in uniform chords.
            Assumes pencil-beams if both omega and width are None. If
            both width and omega None, uses pencil-beam approximation.
        flip : bool (default = True)
            Represents USN as USN on grid if True

        Returns
        -------
        None
        """

        # TODO: Add capability of specifiyng coords according to arrays
        from numpy import ndarray

        # Switch whether to use fan or pencil-beam
        if ((width is None) or (width is False)) and (
            (omega is None) or (omega is False)
        ):
            self.line = True
        else:
            self.line = False
        # Initialize local attributes
        self.case = case
        self.grid = Grid(case, flip=flip)
        self.chords = []
        self.width = width
        self.omega = omega
        # If chords defined, read them
        if isinstance(chords, str):
            self.read_chordfile(chords)
        elif isinstance(chords, (ndarray, list)):
            self.read_chordarray(chords)

    def read_chordarray(self, chords):
        """ " Reads chord data points from a nested list

        len-2 nested list with either 2 (R,Z) or 3 coordinates (x,y,z)
        for start and end points. If 3 coordinates, chords are assumed
        to be pencil-beam projections onto the poloidal plane.

        Returns
        -------
        None
        """
        for chord in chords:
            self.add_chord(chord, self.width)

    def read_chordfile(self, file):
        """Reads chords from file

        Reads the data from file line-by-line, where each coordinate
        is separated by blanks. The lines can contain 4 (2 coordinates
        for start and end points in R,Z-space) or 6 (3 coordinates for
        start and end points in z,y,z-space). If 6 coordinates, chords
        are assumed to be pencil-beam projections onto the poloidal
        plane.

        Returns
        -------
        None
        """
        with open(file, "r") as f:
            for line in f:
                if line.strip()[0] == "#":
                    continue
                parsed = line.strip().split("#")[0].split()
                self.add_chord(
                    [
                        [float(x) for x in parsed[: len(parsed) // 2]],
                        [float(x) for x in parsed[len(parsed) // 2 :]],
                    ],
                    self.width,
                    self.omega,
                )

    def add_chord(self, points, width=None, omega=None):
        """ " Adds a chord object to Spectrometer.chords

        Arguments
        ---------
        points - len-2 list of point (R,Z) or (x,y,z) coordinates.
            if three coordinates given for the points, the chord is
            assumed to be a pencil-beam projection onto the poloidal
            plane.

        Keyword arguments
        -----------------
        width : float (default = None)
            Width of chord fan at end-point. If both width and omega is
            None, values specified on Spectrometer creation are used.
        omega : float (default = None)
            Solid angle of chord. If both width and omega is None,
            values specified on Spectrometer creation are used.

        Returns
        -------
        None
        """
        from shapely.geometry import Point

        # Default to Spectrometer values
        if (width is None) and (omega is None):
            width = self.width
            omega = self.omega
        # Check whether to use pencil beam or fan
        if (len(points[0]) == 3) and (self.line is False):
            print("Tangential chords treated as pencil-beams!")
            self.line = True
        # Create shapely.Point objects from coords
        pointcoords = []
        for point in points:
            pointcoords.append(Point(tuple(point)))
        # Select the appropriate chord object to create
        if self.line:
            Chord = ChordLine
        else:
            Chord = ChordPoly
        # Append created chord to Spectrometer.chords
        self.chords.append(Chord(*pointcoords, self.grid, width=width, omega=omega))

    def add_rates(self, path, species, ratetype, **kwargs):
        """Creates ADASSpecies object as Spectrometer.rates

        Please see uetools.UeUtils.ADASSpecies for more information

        Arguments
        ---------
        path - str to adas directory containing adf* folders
        species - species to gather data for
        ratetype - rate types to include in analysis

        Returns
        -------
        None
        """
        from uetools.UeUtils import ADASSpecies

        self.adas = ADASSpecies(path, species, ratetype, **kwargs)

    def plot_spectrometer(self, ax=None, **kwargs):
        """Plots the spectometer chords on ax

        Passes **kwargs to Chord.plot_chord

        Returns
        -------
        matplotlib.pyploy.Figure
        """
        from matplotlib.pyplot import subplots, Figure, Axes

        # Validate ax input
        if ax is None:
            f, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        elif isinstance(ax, Axes):
            pass
        else:
            raise TypeError("ax type {} not recognized!".format(type(ax)))
        # Plot each chord individually
        for chord in self.chords:
            chord.plot_chord(ax, **kwargs)
        return ax.get_figure()

    def plot_setup(self):
        """Plots the spectrometer chords overlaid on UEDGE grid"""
        f = self.grid.plot_grid()
        self.plot_spectrometer(f)

    def get_nph(self, lam, chargestate, rates=None, rtype=["excit", "recom", "chexc"]):
        """Returns the photon density on the UEDGE grid

        Species is determined by the associated ADASSpecies object

        Arguments
        ---------
        lam - float or photon wavelength associated to chargestate to
            return
        chargestate - int of species charge state that causes emission
        rate - ADASSpecies object associated with species. If none,
            uses Spectrometer.rates

        Keyword arguments
        -----------------
        rtype : list of str (default = ['excit', 'recom', 'chexc']
            Rate types to include in the calculation if they contribute
            to the requested line

        Returns
        -------
        2D array of floats corresponding to ph/s/m**3
        """
        if rates is None:
            rates = self.adas
        species = rates.species.lower()
        return rates.calc_emission(
            self.grid.densities["e"][0],
            self.grid.vars["te"],
            self.grid.densities[species][chargestate],
            self.grid.densities[species][chargestate + 1],
            self.grid.densities["h"][0],
            chargestate,
            lam=lam,
            rtype=rtype,
        )[lam]

    def calc_chord_emission(
        self, chargestate, rates=None, lam=None, rtype=["excit", "recom", "chexc"]
    ):
        """Calcualtes the line-integrated emission along each chord

        Species is determined by the associated ADASSpecies object

        Arguments
        ---------
        chargestate - int of species charge state that causes emission
        rate - ADASSpecies object associated with species. If none,
            uses Spectrometer.rates

        Keyword arguments
        -----------------
        lam : float or list of floats (default = None)
            photon wavelength associated to chargestate to return. If
            None, returns all lines associated with chargestate of
            species
        rtype : list of str (default = ['excit', 'recom', 'chexc']
            Rate types to include in the calculation if they contribute
            to the requested line

        Returns
        -------
        None

        Modidifies
        ----------
        Spectrometer.emission
            Populates the dict with lists of line-integrated intensities
            for the associated wavelengths in ph/s/m**2
        """
        from numpy import array, zeros

        # Reset emission dictionary to avoid double-counting
        self.nph = {}
        self.emission = {}
        # Ensure a rate file is used
        if rates is None:
            rates = self.adas
        # Get the species from the rates
        species = rates.species.lower()
        # Validate lambda input
        # If none, calculate all lines
        if lam is None:
            lam = list(rates.rates[chargestate].keys())
        # If a single line, wrap it in a list
        elif isinstance(lam, (int, float)):
            # If lambda is not found, use the closest line
            if lam not in rates.rates[chargestate]:
                lam = rates.get_closest_line(lam, rates.linelist[chargestate])
            lam = [lam]
        # If list, use all lambdas in list
        elif isinstance(lam, (list, tuple)):
            _lam = []
            for l in lam:
                # Check that line is present, otherwise, use closest lambda
                if l not in rates.rates[chargestate]:
                    l = rates.get_closest_line(l, rates.linelist[chargestate])
                _lam.append(l)
            lam = _lam
        # Perform line-integration for each line and chord
        for l in lam:
            # Get the photon field of the rrequested wavelength(s)
            self.nph[l] = self.get_nph(l, chargestate, rates=rates, rtype=rtype)
            emission_chord = []
            # Perform line-integration for all the chords
            for chord in self.chords:
                emission_chord.append(chord.integrate_field(self.nph[l]))
            # Associate list to dict key
            self.emission[l] = emission_chord

    def plot_chord_integral(
        self, field, ax=None, linestyle="-", marker="o", color="k", x=None
    ):
        """Plots the line-integrated chord values through field

        Arguments
        ---------
        field - 2D field with dimensions compatible with UEDGE case to
            perform line-integration through

        Keyword arguments
        -----------------
        ax : matplotlib.pyplot.Figure or Axes object (default = None)
            Axes to plot onto. If None, creates new figure
        linestyle : str (default = '-')
            style of line plotted
        marker : str (default = 'o')
            marker of data points
        color : str (default = 'k')
            color of line and points
        x : 1D list or array (default = None)
            Plot x-axis of len equals to len(Spectrometer.chords). If
            None, plots chords vs range(len(Spectrometer.chords))

        """
        from matplotlib.pyplot import subplots, Figure, Axes

        if ax is None:
            f, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        elif isinstance(ax, Axes):
            pass
        else:
            raise TypeError("ax type {} not recognized!".format(type(ax)))
        if x is None:
            x = range(1, len(self.chords) + 1)
        # TODO: Figure out what to use as X-axis
        integrals = []
        for chord in self.chords:
            integrals.append(chord.integrate_field(field))

        ax.plot(x, integrals, linestyle=linestyle, marker=marker, color=color)

        return ax.get_figure()

    def plot_chord_spectra(
        self,
        chord,
        chargestate,
        ax=None,
        linestyle="-",
        color="k",
        rates=None,
        **kwargs,
    ):
        """Plots the line-integrated emission spectra for chord

        Calls and passes **kwargs to Spectometer.calc_chord_emission

        Arguments
        ---------
        chord - int of chord in Spectrometer.chords to plot spectra for
        chargestate - int of species charge state that causes emission

        Keyword arguments
        -----------------
        ax : matplotlib.pyplot.Figure or Axes object (default = None)
            Axes to plot onto. If None, creates new figure
        linestyle : str (default = '-')
            style of line plotted
        color : str (default = 'k')
            color of line and points
        rate : ADASSpecies object (default = None)
            Rates associated with species. If none, uses
            Spectrometer.rates
        """
        from matplotlib.pyplot import subplots, Figure

        if ax is None:
            f, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        self.calc_chord_emission(chargestate, **kwargs)
        for line, chords in self.emission.items():
            ax.semilogy(
                [line, line], [0, chords[chord]], linestyle=linestyle, color=color
            )
        return ax.get_figure()

    def plot_chord_emission(
        self,
        lam,
        chargestate,
        ax=None,
        linestyle="-",
        marker="o",
        color="k",
        rates=None,
        x=None,
        **kwargs,
    ):
        """Plots the line-integrated emission spectra of line

        Calls and passes **kwargs to Spectometer.calc_chord_emission

        Arguments
        ---------
        lam - float or photon wavelength associated to chargestate to
            return
        chargestate - int of species charge state that causes emission

        Keyword arguments
        -----------------
        ax : matplotlib.pyplot.Figure or Axes object (default = None)
            Axes to plot onto. If None, creates new figure
        linestyle : str (default = '-')
            style of line plotted
        color : str (default = 'k')
            color of line and points
        marker : str (default = 'o')
            marker of data points
        rate : ADASSpecies object (default = None)
            Rates associated with species. If none, uses
            Spectrometer.rates
        x : 1D list or array (default = None)
            Plot x-axis of len equals to len(Spectrometer.chords). If
            None, plots chords vs range(len(Spectrometer.chords))
        """

        from matplotlib.pyplot import subplots, Figure

        if ax is None:
            f, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        if not isinstance(lam, (float, int)):
            raise Exception("Please select a single spectroscopic line to plot")
        self.calc_chord_emission(chargestate, lam=lam, rates=rates, **kwargs)
        if x is None:
            x = range(1, len(self.chords) + 1)
        # TODO: Figure out what to use as X-axis
        ax.plot(x, self.emission[lam], linestyle=linestyle, marker=marker, color=color)
        return ax.get_figure()




class D3DBolometers:
    """ DIII-D bolometer arrays for projection and back-projection
    
    Based on routines by W.H. Meyer supplied Aug 14th 2024

    Geometry response matrix and routines for DIII-D bolometer diagnostic.

    The website for the diagnostic is at https://diii-d.gat.com/diii-d/Bolo

    Only the first 48 channels implemented by the matrix are in use.
    They represent two fans or 24 channels each. The first 24 are 
    from an R+1 port. The channels from 0-23 sweep from
    looking almost straight down to looking up into the ceiling.
    Channels 24-47 are from a fan located in an R-1 port. The
    channels sweep from looking down to the shelf and sweeping up to 
    look at the ceiling.

    Variable channel_names match the names to channel names used by DIII-D.
    They are also the MDSPlus node names used to read the data. The data
    is located in tree BOLOM and the path is .prad_01.power:<channel_name>.
    Use function read_bolom to read mdsplus data.
    """
    def __init__(self, case, flip=True):
        from  uetools.UeDiagnostics import __path__ as diagnosticspath
        from h5py import File
        from scipy.sparse import csr_matrix
        # Get path to folder containing UeDiagnostics module
        self.matrixdata = {}
        with File('{}/d3d_bolometer_projection.hdf5'.format(
                    diagnosticspath[0])) as f:
            for var, data in f.items(): 
                self.matrixdata[var] = data[()]
        # Set up interpolation finction onto 65x65 EFIT grid
        self.interpdata = lambda x, y: case.utils.squareinterp(
                x,
                r=(self.matrixdata['rmin'], self.matrixdata['rmax']),
                z=(self.matrixdata['zmin'], self.matrixdata['zmax']),
                resolution=(
                        complex(0, self.matrixdata['grid_xpts']),
                        complex(0, self.matrixdata['grid_ypts'])
                ),
                mask = True,
                fill = 0.,
                zshift=y
        )
        # Store the 1D array size
        self.grid_pts = self.matrixdata['grid_xpts']*self.matrixdata['grid_ypts']
        # Create a scipy sparse matrix for bolometer projection matrix data
        self.bolomatrix = csr_matrix((
                self.matrixdata['gmatrix'], 
                (self.matrixdata['gpos'][:,0], self.matrixdata['gpos'][:,1])),
                shape = (self.grid_pts, self.matrixdata['matrix_cols'])
        )
        # Create list of bolometer pointnames
        self.channel_names_n = ['bol_u{:02d}_p'.format(x) for x in range(1,25)]+\
                            ['bol_l{:02d}_p'.format(x) for x in range(1,25)]


    def read_bolom(self, shot,time,server='atlas.gat.com'):
        """
        Read the time point from each of the channel nodes on the given
        shot. Time is given in msecs.

        MDSPlus tree: BOLOM
        Location: .prad_01.power
        Example MDSplus path: \BOLOM::top.prad_01.power.bol_u01_p

        Channels are read in order bol_u01_p to bol_u24_p
        and bol_l01_p to bol_l24_p where the "u" and "l" are
        for upper and lower arrays. 

        Server is defaulted but would require a direct connection. 
        Remote access will require tunneling and a server value of 
        something like 'localhost'. An ssh tunnel might be setup with 
        something like: ssh ...... -L 8000:atlas.gat.com:8000 ......
        """
        chans = np.zeros(self.matrixdata['bolom_chans'])
        try:
            import MDSplus as mds
            try:
                server = mds.Connection(server)
                try:
                    tree = server.openTree('BOLOM',shot)
                except:
                    server.disconnect()
                    print('Failed to open tree')
                    return chans
            except:
                print("Failed to connect to server")
                return chans
        except:
            print("Failed to import MDSplus")
            return chans

        for i in range(self.matrixdata['bolom_chans']):
            try:
                chans[i] = server.get('prad_01.power:{%s}[{%g}]'.format(self.channel_names[i],time))
            except:
                print('Failed to read channel {%d}: {%s}'.format(i,self.channel_names[i]))
           
        server.closeTree('BOLOM',shot)
        server.disconnect()
        return chans

    def project_uedge(self, field, zshift=-1.6, plotinterp=False, **kwargs):
        from numpy import nan_to_num
        from matplotlib.pyplot import subplots
        if plotinterp:
            f, ax = subplots()
            ax.pcolormesh(*[nan_to_num(x) for x in self.interpdata(field, zshift)])
            ax.set_aspect('equal')
        return self.project(nan_to_num(self.interpdata(field, zshift)[-1]), **kwargs)

    def plot_project_uedge(self, field, ax=None, **kwargs):
        from matplotlib.pyplot import subplots, Figure, Axes
        if ax is None:
            f, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
            f = ax.get_figure()
        elif isinstance(ax, Axes):
            f = ax.get_figure()
        else:
            raise TypeError("ax type {} not compatibel".format(type(ax)))
        ax.plot(self.project_uedge(field, **kwargs))

    def project(self, field=None, **kwargs):
        """
        This is the main synthetic diagnotic routine. The input needs to
        be a 65x65 array of radiated power values (W/m^3).  If A is the
        diagnostic response matrix then this returns A*im.

        Use flat=True to generate a flat field projection. This will
        show how diagnostic geometry can effect the channelr-to-channel
        signal variation.
        """
        from numpy import zeros, arange, ones
        if field is None:
            field = ones(self.grid_pts) 
        else:
            field = field.reshape(self.grid_pts)

        return (self.bolomatrix.transpose() * field)[:self.matrixdata['bolom_chans']]

    def backproject(self, signal=None, **kwargs):
        """
        This is the routine is a simple backprojection of the channel 
        data onto a blank image(65x65). If A is the diagnostic response 
        matrix then this return A_transpose*chans

        Use flat=True to generate a flat signal backprojection. This 
        will illuminate the chords in a 65x65 grid. It will show 
        brightness variation due to channel geometry differences.
        """
        from numpy import ones, ndarray, linspace, zeros
        if signal is None:
            lchans = ones(self.matrixdata['matrix_cols'])
            lchans[self.matrixdata['bolom_chans']:self.matrixdata['matrix_cols']] = 0
        elif isinstance(signal, (list, ndarray)):
            lchans = zeros(self.matrixdata['matrix_cols'])
            lchans[:self.matrixdata['bolom_chans']] = signal[:self.matrixdata['bolom_chans']]
        else:
            raise TypeError("Signal type {} not recognized!".format(type(signal)))

        return (self.bolomatrix * lchans).reshape(self.matrixdata['grid_xpts'], self.matrixdata['grid_ypts'])

    def plot_backprojection(self, signal=None, ax=None, res=None, method='linear', **kwargs):
        from matplotlib.pyplot import subplots, Figure, Axes
        from numpy import linspace, meshgrid
        from scipy.interpolate import RegularGridInterpolator
        if ax is None:
            f, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
            f = ax.get_figure()
        elif isinstance(ax, Axes):
            f = ax.get_figure()
        else:
            raise TypeError("ax type {} not compatibel".format(type(ax)))
        ox = linspace(self.matrixdata['rmin'], self.matrixdata['rmax'], 
                self.matrixdata['grid_xpts'])
        oy = linspace(self.matrixdata['zmin'], self.matrixdata['zmax'], 
                self.matrixdata['grid_ypts'])
        oz = self.backproject(signal, **kwargs)
        if res is None:
            ax.pcolormesh(ox, oy, oz, **kwargs)
        else:
            nx = linspace(self.matrixdata['rmin'], self.matrixdata['rmax'], res[0]) 
            ny = linspace(self.matrixdata['zmin'], self.matrixdata['zmax'], res[1]) 
            gx, gy = meshgrid(nx, ny, indexing='ij')
            ax.pcolormesh(gx, gy, RegularGridInterpolator((oy, ox), oz, method=method)((gy, gx)), **kwargs)
            
        ax.set_aspect('equal')
        return f
            





class Grid:
    """Geometric grid for diagnostics

    Attributes
    ----------
    case: uetools.Case object linked to grid
    geometry: str of com.geometry of UEDGE case
    cells: list of Cell objects representing the UEDGE grid
    vars: dict of te, ne, ni, ng from UEDGE on 2D grid
    densities:  species-specific arrays according to
        uetools.UeUtils.AboutSetup.set_speciesarrays() definitions
    area: shapely.Plygon object of the whole UEDGE grid
    core: shapely.Polygon object of the core region
    map: nested dict of cells conforming to UEGDE grid,
        map[ix][iy] = Cell

    Methods
    -------
    plot_grid(ax=None)

    """

    def __init__(self, case, flip=True):
        """Intializes Grid object

        Arguments
        ---------
        case - uetools.Case object to access grid and data from

        Keyword arguments
        -----------------
        flip : bool (default = True)
            Represents USNs in grids with divertor up if True
        """
        from shapely import coverage_union_all, Polygon

        # Link case and required data to object from Case
        self.case = case
        self.geometry = self.case.get("geometry")[0].strip().lower().decode("UTF-8")
        self.cells = []
        self.case.about.set_speciesarrays()
        rm = self.case.get("rm")
        zm = self.case.get("zm")
        # Store basic data to grid
        self.vars = {}
        for var in ["te", "ne", "ni", "ng"]:
            self.vars[var] = self.case.get(var)
        # Check wheter to use flipped grid or not
        if (self.geometry == "uppersn") and (flip is True):
            zm = self.case.disp - zm
        # Parse the data from self.variables to species-resolved
        # data for use with emission
        self.densities = {"e": {0: self.vars["ne"]}}
        for i in range(len(case.about.ionarray)):
            ion = case.about.ionarray[i]
            species = (ion.split("+")[0]).lower()
            if species in ["d", "t"]:
                species = "h"
            try:
                charge = int(ion.split("+")[1])
            except:
                if "+" not in ion:
                    # No charge: inertial neutral species
                    continue
                else:
                    charge = 1
            if species not in self.densities:
                self.densities[species] = {}
            self.densities[species][charge] = self.vars["ni"][:, :, i]
        # Now, repeat for gas arrays
        for g in range(len(case.about.gasarray)):
            gas = case.about.gasarray[g]
            species = (gas.replace("0", "").replace("_2", "")).lower()
            if species in ["d", "t"]:
                species = "h"
            mols = "_2" in gas
            if not mols:
                self.densities[species][0] = self.vars["ng"][:, :, g]
            else:
                self.densities[species]["mol"] = self.vars["ng"][:, :, g]
        # Create dict mapping ix,iy coordinates to Cell objects
        mapdict = {}
        geometries = []
        # Create polygons for all cells
        (nx, ny, _) = rm.shape
        for ix in range(1, nx - 1):
            if ix not in mapdict:
                mapdict[ix] = {}
            for iy in range(1, ny - 1):
                verts = []
                for iv in [1, 2, 4, 3]:
                    verts.append([rm[ix, iy, iv], zm[ix, iy, iv]])
                self.cells.append(Cell(verts, (ix, iy)))
                # Create list of plygons for union
                geometries.append(self.cells[-1].polygon)
                # Map the coordinates to the appropriate Cell object
                mapdict[ix][iy] = self.cells[-1]
        # Create union cells, consisting of the plasma grid
        self.area = coverage_union_all(geometries)
        # Create and store Polygon for core
        corepts = []
        for nxpt in range(case.get("nxpt")):
            for ix in range(case.get("ixpt1")[nxpt] + 1, case.get("ixpt2")[nxpt] + 1):
                corepts.append(
                    [
                        case.get("rm")[ix, 0, 3],
                        case.get("zm")[ix, 0, 3],
                    ]
                )
        self.core = Polygon(corepts)
        self.map = mapdict

    def plot_grid(self, ax=None):
        """Plots a polygrid of Patches on ax"""
        from matplotlib.pyplot import subplots, Figure, Axes

        if ax is None:
            f, ax = subplots()
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        elif isinstance(ax, Axes):
            pass
        else:
            raise TypeError("ax type {} not recognized!".format(type(ax)))
        for cell in self.cells:
            cell.plot_cell(ax)
        self.case.plot.vessel(f)
        self.case.plot.plates(f)

        ax.set_aspect("equal")
        return ax.get_figure()


class Cell:
    """Containter object for grid cell info

    Attributes
    ----------
    vertices: list of 4 (R,Z) pairs of the cell corner nodes
    indices: the (ix, iy) location of cell in UEDGE grid
    polygon: shapely.Polygon representation of the cell

    Methods
    ------
    plot_cell(ax=None, linewidth=0.05)
        plots cell on ax

    """

    # TODO: create function calculating and storing the spectra
    def __init__(self, vertices, indices):
        """Initializes the Cell object

        Arguments
        ---------
        vertices - list of 4 (R,Z) coordinate pairs defning the nodes
        indices - (ix,iy) tuple of cell location in UEDGE index-space
        """
        from shapely.geometry import Polygon

        self.vertices = vertices
        self.indices = indices
        self.polygon = Polygon(vertices)

    def plot_cell(self, ax=None, linewidth=0.05):
        """Plots the cell onto axis"""
        if ax is None:
            f, ax = subplots()
        ax.plot(*self.polygon.exterior.xy, "k-", linewidth=linewidth)
