class ADASSpecies:
    """Class for reading and accessing ADAS data

    Attributes for ratetype=='pec'
    ------------------------------
    rates: nested dict of lines[chargestate][lambda] containing ADF15
        objects
    linelist: dict of lists of lines available linelist[chargestate]

    Methods for ratetype=='pec'
    ---------------------------
    plot_emission(ne, te, n0, n1, nh, chargestate, lam=None,
            rtype=['excit', 'recom', 'chexc'])
        plots the photon emission density for the requested lines
    get_closest_line(lam, linelist)
        returns the value of the line closest to lam
    calc_emssion(ne, te, n0, n1, nh, chargestate, lam=None,
            rtype=['excit', 'recom', 'chexc'])
        calculates the photon emission density for the requested lines
    create_adf15()
        creates and populates the ADASSpecies.rates from adf15 data
    """

    def __init__(self, adaspath, species, ratetype, resolved=False, year=None):
        """Sets up the ADASSpecies class and reads the ADAS data

        Arguments
        ---------
        adaspath - path to local copy of ADAS data (contining adf*'s)
        species - the element that the rates in the object
        ratetype - the type of ADAS rates to be read
            Available options:
                'pec' - adf15 photon emission coefficients

        Keyword arguments
        -----------------
        resolved : bool (default = False)
            switch wheteher to look for metastable resolved data (True)
        year : int (default = None)
            YY int of ADAS data to access. If None, accesses newest
            set of data

        Returns
        -------
        None
        """
        self.path = adaspath
        self.species = species
        self.ratetype = ratetype
        self.year = year
        self.resolved = resolved
        self.rates = {}
        if ratetype.lower() == "pec":
            self.create_adf15()
        else:
            raise Exception(f"Rate type '{ratetype}' not recognized!")

    def plot_emission(
        self,
        ne,
        te,
        n0,
        n1,
        nh,
        chargestate,
        lam=None,
        rtype=["excit", "recom", "chexc"],
    ):
        """Plots the ADAS photon rates in ph/s/m**3

        Calls and passes args to calc_emission

        Arguments
        ---------
        ne - electron density in m**-3
        te - electron temperature in J
        n0 - ion density of radiator at chargestate in m**-3
        n1 - ion density of radiator at chargstate + 1 in m**-3
        chargestate - charge of the species

        Keyword arguments
        -----------------
        lam : float (default = None)
            line to calculate. If None, returns dict of all available
            lines
        rtype : list of str (default = ['excit', 'recom', 'chexc'])
            types of reactions to include in radiation calculation

        Returns
        -------
        matplotlib.pyplot.Figure
        """
        from matplotlib.pyplot import subplots

        f, ax = subplots()
        for x, y in self.calc_emission(
            ne, te, n0, n1, nh, lam, chargestate, rtype
        ).items():
            ax.plot([x, x], [0, y], "r-")
        f.show()
        return f

    def get_closest_line(self, lam, linelist):
        """Returns the line from linelist closeset to lam verbosely"""
        from numpy import array

        # TODO: use local linelist instead?
        print(f"Line at {lam} nm not found.")
        lam = linelist[(abs(array(linelist) - lam).argmin())]
        print(f"Using nearest line found, {lam} nm")
        return lam

    def calc_emission(
        self,
        ne,
        te,
        n0,
        n1,
        nh,
        chargestate,
        lam=None,
        rtype=["excit", "recom", "chexc"],
    ):
        """Calculates the ADAS photon rates in ph/s/m**3

        Calls and passes args to calc_emission

        Arguments
        ---------
        ne - electron density in m**-3
        te - electron temperature in J
        n0 - ion density of radiator at chargestate in m**-3
        n1 - ion density of radiator at chargstate + 1 in m**-3
        chargestate - charge of the species

        Keyword arguments
        -----------------
        lam : float (default = None)
            line to calculate. If None, returns dict of all available
            lines
        rtype : list of str (default = ['excit', 'recom', 'chexc'])
            types of reactions to include in radiation calculation

        Returns
        -------
        Dictionary with the photon rate for the corresponding wavelength
        found from ADAS rates at chargestate and reactions in rtype
        """
        from numpy import array

        if lam is None:
            lines = self.linelist[chargestate]
        elif isinstance(lam, (float, int)):
            if lam not in self.rates[chargestate]:
                lam = self.get_closest_line(lam, self.linelist[chargestate])
            lines = [lam]
        elif isinstance(lam, (list, tuple)):
            lines = []
            for l in lam:
                if l not in self.rates[chargestate]:
                    l = self.get_closest_line(l, self.linelist[chargestate])
                lines.append(l)
        elif isinstance(lam, slice):
            lines = [
                x
                for x in self.linelist[chargestate]
                if ((x > lam.start) and (x < lam.stop))
            ]
        else:
            raise Exception(f"Line specifier option {lam} not recognized!")
        if isinstance(rtype, str):
            rtype = [rtype]
        for typ in rtype:
            if typ not in ["excit", "recom", "chexc"]:
                raise Exception(f"rtype option '{typ}' not implemented!")
        output = {}
        chglines = self.rates[chargestate]
        for lam in lines:
            intensity = 0
            line = chglines[lam]
            for reaction in rtype: 
                if reaction in line.keys():
                    intensity += line[reaction](ne, te, n0, n1, nh, self.resolved)
            output[lam] = intensity
        return output

    def create_adf15(self):
        """Populates the nested dictionary ADASSpecies.rates

        Recursively looks through self.path for adf15 rates and reads
        the requested/newest available rates.
        """
        from os import listdir, path

        # Find all folders in adf15 of path
        folders = [
            f
            for f in listdir("/".join([self.path, "adf15"]))
            if path.isdir("/".join([self.path, "adf15", f]))
        ]
        # Find the folder corresponding to the selected species
        speciesfolder = {}
        for folder in folders:
            try:
                if folder.split("#")[1].lower() == self.species.lower():
                    speciesfolder[int(folder[3:5])] = "/".join(
                        [self.path, "adf15", folder]
                    )
            except:
                pass
        # Check for the requested year, use newest if not specified
        if self.year is not None:
            try:
                self.adfpath = speciesfolder[self.year]
            except:
                print(
                    "ADAS data for requested year not found. Using "
                    + "{}.".format(max(speciesfolder.keys()))
                )
                self.adfpath = speciesfolder[max(speciesfolder.keys())]
        else:
            print(
                "No year of ADAS data requested. Using ADAS{}.".format(
                    max(speciesfolder.keys())
                )
            )
            self.adfpath = speciesfolder[max(speciesfolder.keys())]
        # Find all files in the correct folder
        files = listdir(self.adfpath)
        # Loop through all files, extracting the data to get all lines
        if self.resolved:
            raise Exception("Unresolved rate handling not implemented yet.")
        else:
            for f in files:
                filepath = "/".join([self.adfpath, f])
                # Open every file
                if path.isfile(filepath) and (f.split("#")[1][-1] == "u"):
                    # Get the carge-state
                    chargestate = int(f.split("#")[-1][1])
                    # Create an ADF15 object based on the file
                    rateobj = ADF15(filepath)
                    # Access all lines read from the file via the ADF15 obj
                    for lam, obj in rateobj.adf15.items():
                        # Create a place for the rates in self.rates
                        if not chargestate in self.rates.keys():
                            self.rates[chargestate] = {}
                        rate = {}
                        # Store the different reactions for the line to a buffer
                        for rtype in ["excit", "recom", "chexc"]:
                            try:
                                rate[rtype] = rateobj.adf15[lam][rtype]["emission"]
                            except:
                                pass
                        # Append the rates to the struct
                        self.rates[chargestate][lam] = rate
        # Compile lists of lines by chargestate
        self.linelist = {}
        for chargestate, lines in self.rates.items():
            self.linelist[chargestate] = list(lines.keys())
            self.linelist[chargestate].sort()


class ADF15:
    """Object to read and calculate rates from ADF15 files


    Attributes
    ----------
    file: path to the adf15 file
    metaresolved: bool indicating whether rates are metastably resolved
    chargestate: the charge-state of the radiating elememnt
    adf15: neested dict adf15[lambda][reaction] containing the data

    Methods
    -------
    get_pec(lam, ne, te, rate)
        calculates the photon emission coefficients (ph*cm**3/s)
    get_emission(lam, ne, te, n0, n1, nh, metastate=None,
            rtype=['excit', 'recom', 'chexc'])
        calculates the photon emission rate (ph/s/m**3)
    read_adf15(file, order=1)
        reads the ADF15 file specified self.path

    Returns
    -------
    None
    """

    def __init__(self, fname):
        """Initializes the ADF15 object from file at fname"""
        self.file = fname
        # Read first line
        with open(fname) as f:
            adftype = f.readline()
        # Validate file
        if "PHOTON EMISSIVITY COEFFICIENTS" in adftype:
            self.read_adf15(fname)
        else:
            raise Exception(f"File format of {fname} not recognized!")

    def get_pec(self, lam, ne, te, rate, metastate=None):
        """Retrieves the PEC associated with lam at (ne, Te) of rate type

        Arguments
        ---------
        lam - wavelength in nm
        ne - electron density in m**-3
        te - electron temperature in J
        rate - rate type from options ['excit', 'recom', 'chexc']

        Keyword arguments
        -----------------
        metastate : bool (default = None)
            use metastable resolved rates if True

        Returns
        -------
        PEC at parameters as float in ph*cm**3/s
        """
        from numpy import log10

        if self.metaresolved:
            reaction = self.adf15[lam][rate][metastate]["pecfun"]
        else:
            reaction = self.adf15[lam][rate]["pecfun"]
        return reaction.ev(log10(ne / 1e6), log10(te / 1.602e-19))

    def get_emission(
        self, lam, ne, te, n0, n1, nh, metastate=None, rtype=["excit", "recom", "chexc"]
    ):
        """Retrieves the total photon rate of rtypes in ph/m**3/s

        Arguments
        ---------
        lam - wavelength in angstroms
        ne - electron density in m**-3
        te - electron temperature in J
        n0 - ion density of radiator at chargestate in m**-3
        n1 - ion density of radiator at chargstate + 1 in m**-3

        Keyword arguments
        -----------------
        metastate : bool (default = None)
            use metastable resolved rates if True
        rate : rtype list of strings (default : ['excit', 'recom', 'chexc'])
            reaction types to include in returend rate

        Returns
        -------
        Photon rate for given parameters as float in ph/m**3/s
        """
        if isinstance(rtype, str):
            rtype = [rtype]
        ret = 0
        for rate in rtype:
            if rate not in ["excit", "recom", "chexc"]:
                raise Exception(f"Unknown rate type {rate}!")
            # Grab the densities of the species involved in the process
            nn = (n0 * (rate == "excit") + n1 * (rate != "excit")) * (
                nh * (rate == "chexc") + ne * (rate != "chexc")
            )
            # Get the rates from ADAS, multiply by densities to get
            # volumetric rates, convert to ph/s/cm**3
            ret += 10 ** self.get_pec(lam, ne, te, rate, metastate) * nn * 1e-12
        return ret

    def read_adf15(self, file, order=1):
        """Read photon emissivity coefficients from an ADAS ADF15 file.

        Recursively reads the adf15 file and creates a dictionary at
        self.adf15 whose keys are the wavelengths of the lines
        in angstroms. The value is an interpolant that will evaluate the
        PEC at a desired density and temperature.

        Adapted from AURORA routines

        Units for interpolation: :math:`cm^{-3}` for density;
        :math:`eV` for temperature.

        This function expects the format of PEC files produced via the
        ADAS adas810 or adas218 routines.

        Arguments
        ----------
        file - string with path to file to read

        Keyword arguments
        -----------------
        order : int (default = 1)
            Parameter to control the order of interpolation. Default is
            1 (linear interpolation).

        Returns
        -------
        None

        Modifies
        --------
        Creates the local attributes metaresolved, chargestate, and
        adf15 based on the parsed adf15 file
        """
        from scipy.interpolate import RectBivariateSpline
        from numpy import array, delete, log10, ceil

        if isinstance(file, str):
            self.adf15 = {}
            with open(file, "r") as file:
                header = file.readline()
                self.metaresolved = header.split("/")[1][1].lower() == "r"
                self.chargestate = int(header.split("+")[1].strip()[0:2])
                self.read_adf15(file)
                return
        # adapted from read_adf15 from the Aurora package
        else:
            # Extract info about data block from subheader
            line = file.readline()
            while (line[0].lower() == "c") or (len(line.strip()) == 0):
                line = file.readline()
                if len(line) == 0:
                    return
            [data, filmem, rtype, indm, isel] = line.split("/")
            lam, res = data.split("A")
            lam = float(lam)
            [ndens, nte] = [int(x) for x in res.split()]
            indm = indm.split("=")[1].strip()
            isel = int(isel.split("=")[1].strip())
            rtype = rtype.split("=")[1].strip().lower()
            filmem = filmem.split("=")[1].strip()
            data = []
            while len(data) < ndens + nte + ndens * nte:
                line = file.readline()
                if ("----" not in line) and (line[0].lower() != "c"):
                    [data.append(float(x)) for x in line.split()]
            dens = data[:ndens]
            te = data[ndens : ndens + nte]
            pec = array(data[-ndens * nte :]).reshape((ndens, nte))
            if not lam in self.adf15.keys():
                self.adf15[lam] = {}
            try:
                block = self.adf15[lam][rtype]
            except:
                self.adf15[lam][rtype] = {}
                block = self.adf15[lam][rtype]
            if self.metaresolved:
                adf15[lam][rtype][f"{indm}_{isel}"] = {}
                block = adf15[lam][rtype][f"{indm}_{isel}"]
            # Store data
            for var in ["ndens", "nte", "indm", "isel", "filmem", "te", "dens", "pec"]:
                block[var] = locals()[var]
            block["pecfun"] = RectBivariateSpline(
                log10(dens), log10(te), log10(pec), kx=order, ky=order
            )

            fun = exec(
                f"lambda *args, **kwargs: self.get_emission({lam}, *args, rtype={rtype}, **kwargs)"
            )
            block["emission"] = lambda *args, **kwargs: self.get_emission(
                lam, *args, rtype=rtype, **kwargs
            )
            self.read_adf15(file)
            return
