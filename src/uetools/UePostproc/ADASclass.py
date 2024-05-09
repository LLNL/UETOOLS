
class ADASSpecies:
    def __init__(self, adaspath, species, ratetype, resolved=False, year=None):
        self.path = adaspath
        self.species = species
        self.ratetype = ratetype
        self.year = year
        self.resolved = resolved
        self.lines = {}
        self.rates = {}
        if ratetype.lower() == 'pec':
            self.create_adf15()
        else:
            raise Exception("Rate type '{ratetype}' not recognized!")
       
    def plot_emission(self, ne, te, nz, lam=None, chargestate=None, 
        rtype=['excit', 'recom', 'chexc']
        ):
        from matplotlib.pyplot import subplots
        f, ax = subplots()
        for x, y in self.calc_spectra(ne, te, nz, lam, chargestate, rtype).items():
            ax.plot([x,x], [0,y], 'r-')
        f.show()
        return f

    def calc_emission(self, ne, te, nz, nh, lam=None, chargestate=None,
        rtype = ['excit', 'recom', 'chexc']
        ):
        if lam is None:
            lines = self.linelist
        elif isinstance(lam, float):
            try:
                self.lines[lam]
                lines = lam
            except:
                raise Exception(f"Line at {lam} A not found.")
        elif isinstance(lam, (list, tuple)):
            # TODO: Will cause unexpected behavior if asking for 2 distict lines
            if len(lam) == 2:
                lines = [x for x in self.lines if ((x > lam[0]) and (x < lam[1]))]
            else:
                lines = []
                for l in lam:
                    try:
                        self.lines[l]
                        lines.append(l)
                    except:
                        raise Exception(f"Line at {lam} A not found.")
        else:
            raise Exception(f"Line specifier option {lam} not recognized!")
        if isinstance(rtype, str):
            rtype = [rtype]
        for typ in rtype:
            if typ not in ['excit', 'recom', 'chexc']:
                raise Exception(f"rtype option '{typ}' not implemented!")
        if isinstance(chargestate, int):
            chargestate = [chargestate]
        output = {}
        for lam in lines:
            intensity = 0
            line = self.lines[lam]
            for zi, rate in line.items():
                if (chargestate is not None) and (zi not in chargestate):
                    continue
                else:
                    for r in rtype:
                        try:
                            intensity += line[zi][r](ne, te, nz, self.resolved) 
                        except:
                            pass
            if intensity != 0:
                output[lam] = intensity
        return output

    def create_adf15(self):
        from os import listdir, path
        folders = [f for f in listdir('/'.join([self.path, 'adf15'])) \
                if path.isdir('/'.join([self.path, 'adf15', f]))]
        speciesfolder = {}
        for folder in folders:
            try:
                if folder.split('#')[1].lower() == self.species.lower():
                    speciesfolder[int(folder[3:5])] = '/'.join([self.path, 'adf15', folder])
            except:
                pass
        if self.year is not None:
            try:
                self.adfpath = speciesfolder[self.year]
            except:
                print("ADAS data for requested year not found. Using " + \
                    "{}.".format(max(speciesfolder.keys()))) 
                self.adfpath = speciesfolder[max(speciesfolder.keys())]
        else:
            print("No year of ADAS data requested. Using ADAS{}.".format(\
                max(speciesfolder.keys()))) 
            self.adfpath = speciesfolder[max(speciesfolder.keys())]
        
        files = listdir(self.adfpath)
        if self.resolved:
            raise Exception("Unresolved rate handling not implemented yet.")
        else:
            for f in files:
                filepath = '/'.join([self.adfpath, f])
                if path.isfile(filepath) and (f.split('#')[1][-1]=='u'):
                    chargestate = int(f.split("#")[-1][1])
                    rateobj = ADASRate(filepath)
                    try: 
                        self.rates[chargestate].append(rateobj)
                    except:
                        self.rates[chargestate] = []
                        self.rates[chargestate].append(rateobj)
                    for lam, obj in rateobj.adf15.items():
                        if not lam in self.lines.keys():
                            self.lines[lam] = {}
                        if not chargestate in self.lines[lam].keys():
                            self.lines[lam][chargestate] = []
                        rate = {}
                        for rtype in ['excit', 'recom', 'chexc']:
                            try:
                                rate[rtype] = rateobj.adf15[lam][rtype]['emission']
                            except:
                                pass
                        self.lines[lam][chargestate] = rate
        self.linelist = list(self.lines.keys())
        self.linelist.sort()


class ADASRate:
    def __init__(self, fname):
        
        self.file = fname

        with open(fname) as f:
            # TODO: verify style of ADF15 and other files
            adftype = f.readline()

        if "PHOTON EMISSIVITY COEFFICIENTS" in adftype:
            self.read_adf15(fname)
        else:
            raise Exception(f"File format of {fname} not recognized!")
      
    def get_pec(self, lam, ne, te, rate, metastate = None):
        from numpy import log10
        if self.metaresolved:
            reaction = self.adf15[lam][rate][metastate]['pecfun']
        else:
            reaction = self.adf15[lam][rate]['pecfun']
        return reaction.ev(log10(ne/1e6), log10(te/1.602e-19))

    def get_emission(self, lam, ne, te, nz, nh, metastate=None, rtype=['excit', 'recom', 'chexc']):
        if isinstance(rtype, str):
            rtype = [rtype]
        ret = 0
        for rate in rtype:
            if rate == 'chexc':
                try:
                    nn = nz[self.chargestate+1]*nh
                except: 
                    print("WARNING")
                    nn = 1
            elif rate == 'excit':
                nn = nz[self.chargestate]*ne
            elif rate == 'recom':
                try:
                    nn = nz[self.chargestate+1]*ne
                except:
                    nn = 1
                    print("WARNING")
            else:   
                raise Exception(f"Unknown rate type {rate}!")
            try:
                ret += 10**self.get_pec(lam, ne, te, rate, metastate)*nn*1e-12
            except:
                pass
        return ret


    def read_adf15(self, file, order=1):
        """Read photon emissivity coefficients from an ADAS ADF15 file.

        Returns a dictionary whose keys are the wavelengths of the lines in angstroms.
        The value is an interpolant that will evaluate the PEC at a desired density and temperature.

        Units for interpolation: :math:`cm^{-3}` for density; :math:`eV` for temperature.

        Parameters
        ----------
        path : str
            Path to adf15 file to read.
        order : int, opt
            Parameter to control the order of interpolation. Default is 1 (linear interpolation).
        Returns
        -------
        None

        Notes
        -----
        This function expects the format of PEC files produced via the ADAS adas810 or adas218 routines.

        """
        from scipy.interpolate import RectBivariateSpline
        from numpy import array, delete, log10, ceil

        if isinstance(file, str):
            self.adf15 = {}
            with open(file, "r") as file:
                header = file.readline()
                self.metaresolved = (header.split('/')[1][1].lower() == 'r')
                self.chargestate = int(header.split('+')[1].strip()[0:2])
                self.read_adf15(file)
                return
        # adapted from read_adf15 from the Aurora package
        else:
            # Extract info about data block from subheader
            line = file.readline()
            while (line[0].lower() == 'c') or (len(line.strip())==0):
                line = file.readline()
                if len(line)==0:
                    return
            [data, filmem, rtype, indm, isel ] = line.split('/')
            lam, res = data.split('A')
            lam = float(lam)
            [ndens, nte] = [int(x) for x in res.split()]
            indm = indm.split('=')[1].strip()
            isel = int(isel.split('=')[1].strip())
            rtype = rtype.split('=')[1].strip().lower()
            filmem = filmem.split('=')[1].strip()
            data = []
            while len(data) < ndens + nte + ndens*nte:
                line = file.readline()
                if ("----" not in line) and (line[0].lower() != 'c'):
                    [data.append(float(x)) for x in line.split()]
            dens = data[:ndens]
            te = data[ndens:ndens+nte]
            pec = array(data[-ndens*nte:]).reshape((ndens, nte))
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
            for var in ['ndens', 'nte', 'indm', 'isel', 'filmem', 'te', 
                'dens', 'pec'
            ]:
                block[var] = locals()[var]
            block['pecfun'] = RectBivariateSpline(
                log10(dens), log10(te), log10(pec), kx=order, ky=order
            )

            fun = exec(f'lambda *args, **kwargs: self.get_emission({lam}, *args, rtype={rtype}, **kwargs)')
            block['emission'] = (lambda *args, **kwargs: self.get_emission(lam, *args, rtype=rtype, **kwargs))
            self.read_adf15(file)
            return

#def get_emission(self, lam, ne, te, nz, nh, metastate=None, rtype=['excit', 'recom', 'chexc']):

