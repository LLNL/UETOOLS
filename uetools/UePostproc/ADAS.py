class ADAS:
    def emission_CIII(self, fname):
        from numpy import log10

        # Create self.adf15
        self.read_adf15(fname)

        for rate in ["excit", "recom", "chexc"]:
            self.__setattr__(
                "CIII_pec_{}".format(rate),
                self.adf15[4650.1][rate].ev(
                    log10(self.get("ne") / 1e6), log10(self.get("te") / 1.602e-19)
                ),
            )
            self.__setattr__(
                "CIII_emission_{}".format(rate),
                10 ** self.__getattribute__("CIII_pec_{}".format(rate))
                * (self.get("ne")*(rate != "chexc") + self.get("ng")[:,:,0]*(rate == "chexc"))
                * self.get("ni")[:, :, 4 - (rate == "excit")]
                / 1e12,
            )
        self.CIII_emission = (
            self.CIII_emission_excit
            + self.CIII_emission_recom
            + self.CIII_emission_chexc
        )
        
    def i_emiss(self, fname, i_idx, wavel, plot=True, **kwargs):
        from numpy import log10

        # Create self.adf15
        self.read_adf15(fname)
        t = self.get_tname(i_idx)
        ne = self.get('ne')
        te = self.get('te')
        ni = self.get('ni')[:,:,i_idx]
        niup = self.get('ni')[:,:,i_idx+1]
        ng = self.get('ng')[:,:,0]
        j2ev = 1.6022e-19
        
        self.__setattr__("{}_pec_excit".format(t),self.adf15[wavel]["excit"].ev(log10(ne/1e6),log10(te/j2ev)))
        self.__setattr__("{}_pec_recom".format(t),self.adf15[wavel]["recom"].ev(log10(ne/1e6),log10(te/j2ev)))
        self.__setattr__("{}_pec_chexc".format(t),self.adf15[wavel]["chexc"].ev(log10(ne/1e6),log10(te/j2ev)))
        
        self.__setattr__("{}_emission_excit".format(t),10**self.__getattribute__("{}_pec_excit".format(t))*ne*ni/1e12)
        self.__setattr__("{}_emission_recom".format(t),10**self.__getattribute__("{}_pec_recom".format(t))*ne*niup/1e12)
        self.__setattr__("{}_emission_chexc".format(t),10**self.__getattribute__("{}_pec_chexc".format(t))*ng*niup/1e12)
        
        tot = self.__getattribute__("{}_emission_excit".format(t)) + self.__getattribute__("{}_emission_recom".format(t)) + self.__getattribute__("{}_emission_chexc".format(t))
        self.__setattr__("{}_emission".format(t),tot)
        
        if plot:
            title = "{} emission".format(t)
            return self.plotmesh(tot,title=title,**kwargs)
        else:
            return tot
        
    def get_tname(self,i_idx):
        tname = str()
        i_mass = self.get('minu')[i_idx]
        if i_mass == 12:
            tname += 'C'
        elif i_mass == 4:
            tname += 'He'
        elif i_mass == 2:
            tname += 'D'
            
        i_charge = int(self.get('ziin')[i_idx])
        if i_charge > 0:
            tname += str(i_charge)+'+'
        else:
            tname += str(i_charge)
            
        return tname

    def read_adf15(self, path, order=1, raw=False):
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
        from numpy import array, delete, log10

        # adapted from read_adf15 from the Aurora package

        # find out whether file is metastable resolved
        meta_resolved = path.split("#")[-2][-1] == "r"
        if meta_resolved:
            print("Identified metastable-resolved PEC file")

        with open(path, "r") as f:
            lines = f.readlines()
        cs = path.split("#")[-1].split(".dat")[0]

        header = lines.pop(0)
        # Get the expected number of lines by reading the header:
        num_lines = int(header.split()[0])
        self.adf15 = {}

        for i in range(0, num_lines):
            if "----" in lines[0]:
                _ = lines.pop(0)  # separator may exist before each transition
            # Get the wavelength, number of densities and number of temperatures
            # from the first line of the entry:
            l = lines.pop(0)
            header = l.split()
            # sometimes the wavelength and its units are not separated:
            try:
                header = [hh.split("A")[0] for hh in header]
            except:
                # lam and A are separated. Delete 'A' unit.
                header = delete(header, 1)
            lam = float(header[0])
            if header[1] == "":
                # 2nd element was empty -- annoyingly, this happens sometimes
                num_den = int(header[2])
                num_temp = int(header[3])
            else:
                num_den = int(header[1])
                num_temp = int(header[2])
            if meta_resolved:
                # index of metastable state
                INDM = int(header[-3].split("/")[0].split("=")[-1])
            # Get the densities:
            dens = []
            while len(dens) < num_den:
                dens += [float(v) for v in lines.pop(0).split()]
            dens = array(dens)
            # Get the temperatures:
            temp = []
            while len(temp) < num_temp:
                temp += [float(v) for v in lines.pop(0).split()]
            temp = array(temp)
            # Get the PEC's:
            PEC = []
            while len(PEC) < num_den:
                PEC.append([])
                while len(PEC[-1]) < num_temp:
                    PEC[-1] += [float(v) for v in lines.pop(0).split()]
            PEC = array(PEC)
            # find what kind of rate we are dealing with
            if "recom" in l.lower():
                rate_type = "recom"
            elif "excit" in l.lower():
                rate_type = "excit"
            elif "chexc" in l.lower():
                rate_type = "chexc"
            elif "drsat" in l.lower():
                rate_type = "drsat"
            # attempt to report unknown rate type -- this should be fairly robust
            else:
                rate_type = l.replace(" ", "").lower().split("type=")[1].split("/")[0]
            # create dictionary with keys for each wavelength:

            # add a key to the pec_dict[lam] dictionary for each type of rate: recom, excit or chexc
            # interpolate PEC on log dens,temp scales
            PEC_use = log10(PEC)
            if raw is True:
                PEC_use = PEC
            pec_fun = RectBivariateSpline(
                log10(dens), log10(temp), PEC_use, kx=order, ky=order
            )

            try:
                self.adf15[lam]
            except:
                self.adf15[lam] = dict()
            if meta_resolved:
                try:
                    self.adf15[lam][rate_type]
                except:
                    self.adf15[lam][rate_type] = dict()
                self.adf15[lam][rate_type]["meta{}".format(INDM)] = pec_fun
            else:
                self.adf15[lam][rate_type] = pec_fun
                if raw is True:
                    self.adf15[lam][rate_type] = PEC
