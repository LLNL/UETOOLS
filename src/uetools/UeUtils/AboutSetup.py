class AboutSetup:
    """Class outputting the UEDGE setup in memory

    Note: this class is still a stub

    Attributes
    ----------
    ionarray:
    gasarray:

    Methods
    -------
    elements()
        returns a nested dict of elemnts and isotopes with
        symbols as values: dict[Z][A] = str
    uedge_setup()
        TBD - prints info on the current UEGDE setup
    recycling_snull()
        TBD - prints info on the current snull recucling setup
    equations()
        TBD - prints info on the current equation setup
    sputtering()
        TBD - prints into on the current sputtering setup
    set_speciesarray()
        sets self.ionarray and self.gasarray
    species_setup()
        prints information on models and species
    """

    def __init__(self, case):
        self.get = case.get

    #        self.species_setup()

    def elements(self):
        return {
            1: {1: "H", 2: "D", 2.5: "DT", 3: "T"},
            2: {3: "He3", 4: "He4"},
            3: {6: "Li6", 7: "Li7"},
            4: {7: "Be7", 8: "Be8", 9: "Be9", 10: "Be10"},
            6: {11: "C11", 12: "C12", 13: "C13", 14: "C14"},
            7: {13: "N13", 14: "N14", 15: "N15"},
            8: {15: "O15", 16: "O16", 17: "O17", 18: "O18"},
            10: {20: "Ne20", 21: "Ne21", 22: "NE22"},
            18: {
                36: "Ar36",
                37: "Ar37",
                38: "Ar38",
                39: "Ar39",
                40: "Ar40",
                41: "Ar41",
                42: "Ar42",
            },
        }

    def uedge_setup(self):
        self.species_setup()
#        self.about_equations()
#        self.about_recycling()
#        self.about_sputtering()

    def recycling_snull(self):
        from sys import modules

        raise NotImplementedError()

        for var in ["albedo_by_user"]:
            setattr(modules[__name__], var, self.get(var))
        print("Left plate")
        print("Right plate")
        print("PFR boundary")
        print("Core boundary")
        print("Outer boundary")

    def equations(self):
        raise NotImplementedError()

    def sputtering(self):
        raise NotImplementedError()

    def set_speciesarrays(self):
        """Sets ionarray and gasarray with lists of species strings"""
        from sys import modules
        from numpy import count_nonzero

        # Initialize local variables
        for var in [
            "ngsp",
            "nhsp",
            "nisp",
            "nzsp",
            "minu",
            "znuclin",
            "ishymol",
            "nhgsp",
            "zi",
            "mp",
            "mg",
        ]:
            setattr(modules[__name__], var, self.get(var))
        # Get the elemental arrays
        elements = self.elements()
        # Total number of hydrogenic species
        try:
            self.nhydrogenic = nhsp + ishymol
        except:
            self.nhydrogenic = nhgsp + 1
        # Total number of impurity species
        self.nimp = count_nonzero(nzsp)
        # Get arrays for ion species
        ionarray = []
        for i in range(nisp):
            charge = "{}".format(
                ("+" + str(int(zi[i]))) * (zi[i] > 0) + "0" * (zi[i] == 0)
            )
            sign = "".join(
                (x for x in elements[znuclin[i]][minu[i]] if not x.isdigit())
            )
            ionarray.append("{}{}".format(sign, charge))  # .center(5))
        for i in range(nhsp - 1):
            ionarray[i] = ionarray[i][:-1]
        # Get arrays for gaseous species
        gasarray = []
        for i in range(ngsp):
            mols = False
            z = zi[0]
            if i > nhgsp - 1:
                z = znuclin[nhsp + sum(nzsp[: i - nhgsp + 1]) - 1]
            if (ishymol != 0) and (i == 1):
                mols = True
                z = zi[0]
                i = 0
            sign = "".join((x for x in elements[z][int(mg[i] / mp)] if not x.isdigit()))
            if mols:
                sign = sign + "_2"
            else:
                sign = sign + "0"
            gasarray.append(sign)  # .center(5))
        self.ionarray = ionarray
        self.gasarray = gasarray

        # TODO: add option for detecting mixes

    def species_setup(self):
        """Prints the species included and models used"""
        from sys import modules

        # Initialize local variables
        for var in [
            "ngsp",
            "nhsp",
            "nisp",
            "nusp",
            "nzsp",
            "minu",
            "znuclin",
            "isngon",
            "isupgon",
            "ishymol",
            "nhgsp",
            "isimpon",
            "atn",
            "atw",
            "zi",
            "mp",
            "mg",
            "nusp_imp",
        ]:
            setattr(modules[__name__], var, self.get(var))
        self.set_speciesarrays()
        elements = self.elements()
        # Get atomic transport model
        if (isngon[0] == 1) and (isupgon[0] == 0):
            atommodel = "Diffusive "
        elif (isngon[0] == 0) and (isupgon[0] == 1):
            atommodel = "Inertial "
        else:
            atommodel = "Undefined transport model "
        # Hydrogen key
        hysign = elements[znuclin[0]][minu[0]]
        # Impurity species signs
        impsign = []
        # Get resolved impurity
        for iimp in range(self.nimp):
            impind = nhsp + sum(nzsp[: iimp + 1]) - 1
            impsign.append(elements[znuclin[impind]][minu[impind]])
        # Get impurity transport model
        imp_transp = []
        for iimp in range(self.nimp):
            if nusp_imp == 0:
                imp_transp.append("Fully force-balance")
            elif nusp_imp >= sum(nzsp[: iimp + 1]):
                imp_transp.append("Fully momentum-resolved")
            else:
                sp_nusp_imp = nusp_imp - sum(nzsp[:iimp])
                imp_transp.append(
                    ("{} momentum-resolved and {} " + "force-balance species").format(
                        sp_nusp_imp, nzsp[iimp] - sp_nusp_imp
                    )
                )
        # Get potential force-balance information
        ffimp = isimpon in [2, 7]
        if ffimp:
            self.nimp += 1
            ffsign = elements[atn][atw]
        ionarray = [x.center(5) for x in self.ionarray]
        gasarray = [x.center(5) for x in self.gasarray]
        idxarr = [str(x).center(5) for x in range(max(len(ionarray), len(gasarray)))]
        # Output species setup
        print("The UEDGE set-up contains:")
        print("  - {} hydrogenic species:".format(self.nhydrogenic))
        print("    - {} ions".format(hysign))
        if nhsp > 1:
            print("    - {}{} atoms".format(atommodel, hysign))
        if ishymol != 0:
            print("    - {}2 molecules".format(hysign))
        if self.nimp > 0:
            print("  - {} impurity species:".format(self.nimp))
            for i in range(len(impsign)):
                print("    - Charge-state resolved {}".format(impsign[i]))
                print("      - {}".format(imp_transp[i]))
            if ffimp:
                print("    - Fixed-fraction {}".format(ffsign))
        print(
            "\n{}{}".format(
                13 * " ", str(idxarr)[1:-1].replace("'", "").replace(",", " ")
            )
        )
        print("Ion array : {}".format(str(ionarray).replace("'", "")))
        print("Gas array : {}".format(str(gasarray).replace("'", "")))

        return
