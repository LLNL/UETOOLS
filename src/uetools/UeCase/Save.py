from uetools.UeUtils import Lookup


class Save:
    """Class providing UETOOLS with save and restore functionality


    Methods
    -------
    diffusivities(savefname, **kwargs)
        Saves the current UEDGE diffusivities to savefname
    dump(savefname, **kwargs)
        Dumps all UEDGE variables to savefname HDF5. Calls Save.save
    grid(savefname, **kwargs)
        Saves the current UEDGE grid to savefname
    group(savename, group, append=True, **kwargs)
        Saves the group specified by group to savefile
    load_state(savefname=None, **kwargs)
        Restores UEDGE state variables from savefname, or Case if None
    metadata(savefile, **kwargs)
        Saves the Case metadata to savefile
    recursive(savefile, saveobj, group=[], **kwargs)
        Recursively saves saveobj to savefile
    save(   savefname, group=None, append=False,
            postprocess=True, **kwargs)
        Saves data from group, or all data if group is None, to
        savefname. Variables defined by YAML files
    setup(savefname, **kwargs)
        Saves the current UETOOLS setup to savefname
    var(savefile, groups, variable, data, **kwargs)
        Saves variable to savefile
    """

    def __init__(self, case):
        """Initializes Save object and links to Case object"""
        self.get = case.get
        self.getue = case.getue
        self.setue = case.setue
        self.info = case.info
        self.mutex = case.mutex
        self.variables = case.variables
        self.record_changes = case.tracker.record_changes
        self.getpackobj = Lookup().getpackobj

    def var(self, savefile, groups, variable, data, **kwargs):
        """Saves variable and metadata to HDF5 group and dataset

        Arguments
        ---------
        savefile : h5py File object
            h5py object where to save data
        groups : list of str
            list of subroups to create and store variable
        variable : str
            name of dataset to create. If same as UEDGE variable name,
            metadata (description, units) are stored as dataset
            attributes
        data : array/int/str/float
            data to be stored to dataset

        Returns
        -------
        None
        """
        from Forthon import packageobject

        # Make group, package, file into a list and iterate?
        output = []
        for group in groups:
            savefile.require_group(group)
            savefile = savefile[group]
            output.append(group)
        if variable not in savefile.keys():
            # Check whether variable is allocated or not
            try:
                if self.getpackobj(variable, verbose=False).allocated(variable) != 0:
                    savefile.create_dataset(variable, data=data)
                else:
                    savefile.create_dataset(variable, data=0)
                    # TODO: Omit or save False?
                    # TODO: Verbose prompt or not?
            except:
                savefile.create_dataset(variable, data=data)
        try:
            savefile.attrs[variable] = packageobject(package).getvarunit(variable)
        except:
            pass
        try:
            savefile.attrs[variable] = packageobject(package).getvardoc(variable)
        except:
            pass

    def group(self, savename, group, append=True, **kwargs):
        """Write group to HDF5 file savename

        Calls Save.recursivesave

        Arguments
        ---------
        savename : str of path to file to write
        group :  str with name of group to write

        Keyword arguments
        -----------------
        append : bool (default = True)
            Switch whether to append to file (True) or rewrite

        Returns
        -------
        None
        """
        from h5py import File

        with File(savename, "a" if append else "w") as savefile:
            if append is False:
                self.metadata(savefile)
            self.recursive(savefile, self.variables["input"][group], [group])

    def setup(self, savename, **kwargs):
        """Saves setup to savename using Save.group"""
        self.group(savename, "setup", **kwargs)

    def grid(self, savename, **kwargs):
        """Saves grid to savename using Save.group"""
        self.group(savename, "grid", **kwargs)

    def diffusivities(self, savename, **kwargs):
        """Saves diffusivities to savename using Save.group"""
        self.group(savename, "diffusivities")

    def recursive(self, savefile, saveobj, group=[], **kwargs):
        """Recursively writes data to File object savefile


        Arguments
        ---------
        savefile : open File object with write access
        saveobj: struct of data to write (dict, list, value)

        Keyword arguments
        -----------------
        group : list (default = [])
            nested list of strings defining position where to write in
            savefile

        Returns
        -------
        None
        """
        # Bottom level of structure
        if not isinstance(saveobj, dict):
            # Special setup for saving setup parameters to store
            # actual set values of any setup parameters defined in input
            if group[0] == "setup":
                variable = group.pop(-1)
                if variable in [
                    "userdifffname",
                    "radialdifffname",
                    "casename",
                    "commands",
                    "savefile",
                ]:
                    value = saveobj
                # Exception for saving objects spawned without any input files
                #                elif isinstance(saveobj, list):
                #                    for var in saveobj:
                #                        self.var(savefile, group, var, self.getue(var))
                #                    return
                else:
                    if variable in self.variables['omit']:
                        return
                    else:
                        try:
                            value = self.getue(variable)
                        except:
                            for var in saveobj:
                                for var in saveobj:
                                    self.var(savefile, group, var, self.getue(var))
                                return
                self.var(savefile, group, variable, value)
            # Store requested data
            elif isinstance(saveobj, list):
                for variable in saveobj:
                    self.var(savefile, group, variable, self.getue(variable))
        # Recursively go deeper in structure
        else:
            for key, value in saveobj.items():
                if isinstance(key, int):
                    # Save the full array once only
                    variable = group.pop(-1)
                    self.var(savefile, group, variable, self.getue(variable))
                    return
                saveobj = self.recursive(savefile, value, group + [key])
            return saveobj

    def metadata(self, savefile, **kwargs):
        """Writes metadata to File object savefile"""
        from time import time, ctime

        try:
            savefile.attrs["casename"] = self.info["casename"]
        except:
            pass
        savefile.attrs["UETOOLS_version"] = self.info["uetoolsversion"]
        savefile.attrs["time"] = time()
        savefile.attrs["ctime"] = ctime()
        savefile.attrs["code"] = "UEDGE"
        savefile.attrs["ver"] = self.info["uedge_ver"]
        savefile.attrs["pyver"] = self.info["pyver"]
        savefile.attrs["user"] = self.info["user"]
        savefile.attrs["hostname"] = self.info["hostname"]
        try:
            savefile.attrs["location"] = self.info["location"]
        except:
            pass

    def save(self, savefname, group=None, append=False, postprocess=True, **kwargs):
        """Saves HDF5 file containing UeCase data

        Arguments
        ---------
        savefname : str
            path to file to write data to

        Keyword arguments
        -----------------
        group : str (default = None)
            group identifier of group to be written to file. If None,
            all data stored in UeCase is written
        append : bool (default = False)
            switch whether to append data to file (True) or rewrite
        postprocess : bool (default = True)
            switch whether to postprocess UEDGE case for power and
            particle balance/fluxes before storing (True)

        Returns
        -------
        None
        """
        from h5py import File
        from Forthon import packageobject

        if self.info["inplace"]:
            print("Data read from file, no data to save. Aborting.")
            return
        if postprocess is True:
            bbb = packageobject("bbb")
            bbb.engbal(self.get("pcoree") + self.get("pcorei"))
            bbb.plateflux()
            bbb.wallflux()
        # Check and store any changes since case last saved/read
        self.record_changes()
        # Open file for writing
        with File(savefname, "a" if append else "w") as savefile:
            # Save metadata of case
            self.metadata(savefile)
            # Write requested data to file
            if group is None:
                self.recursive(savefile, self.variables["input"])
            else:
                self.recursive(savefile, self.variables["input"][group], [group])

    def dump(self, savefname, **kwargs):
        """Dumps all variables to savefname

        Calls self.save(**kwargs)

        Arguments
        ---------
        savefname : str
            Path to file to write

        Returns
        -------
        None
        """
        # TODO: compress?
        try:
            from Forthon import package, packageobject
        except:
            pass
        from h5py import File

        self.save(savefname, **kwargs)
        with File(savefname, "a") as savefile:
            # Create data dump class
            savefile.create_group("datadump")
            for pkg in package():
                for var in packageobject(pkg).varlist():
                    self.var(savefile, ["datadump"], var, self.getue(var))

    def load_state(self, savefname=None, **kwargs):
        """Restores a saved solution

        Keyword arguments
        -----------------
        savefname : str (default = None)
            HDF5 file to read stored solution from. If None, solution
            is read from UeCase object
        **kwargs
            passed to setgroup

        Returns
        -------
        None
        """
        from h5py import File, is_hdf5
        from os.path import exists

        if self.mutex() is False:
            return

        if savefname is None:
            savefname = "{}.hdf5".format(self.info["casename"])
        if not exists(savefname):
            raise ValueError("Save file {} does not exist!".format(savefname))
        elif not is_hdf5(savefname):
            raise ValueError("Save file {} is not in HDF5 format!".format(savefname))

        with File(savefname, "r") as savefile:
            #        # Try reading new, subdivided save file
            #        try:
            #            # Don't override user-specified name for case by reading from file
            #            if casefname is None:
            #                self.casename = savefile.attrs["casename"]
            #        except:
            #            pass
            if "restore" in savefile.keys():
                for group, variables in savefile["restore"].items():
                    for variable, value in variables.items():
                        self.setue(variable, value[()])
                        self.variables["stored"][variable] = value[()]
                # If not, try reading old-style save file
                try:
                    for group, variables in savefile["restore"].items():
                        for variable, value in variables.items():
                            self.setue(variable, value[()])
                            self.variables["stored"][variable] = value[()]
                # If not, try reading old-style save file
                except:
                    for group, variables in self.variables["input"]["restore"].items():
                        for variable in variables:
                            self.setue(variable, savefile[group][variable][()])
                            self.variables["stored"][variable] = savefile[group][
                                variable
                            ][()]
                if self.info["verbose"]:
                    prfile = savefname
                    if len(prfile.split('/'))>3:
                        prfile = ".../{}".format('/'.join(prfile.split('/')[-3:]))
                    print(
                        "UETOOLS-style save successfully restored "
                        + "from {}".format(prfile)
                    )
            elif "bbb" in savefile.keys():
                statevars = ["nis", "ngs", "tes", "tis", "ups", "phis"]
                for var in ["tipers", "tgs"]:
                    try:
                        self.getue(var)
                        statevars.append(var)
                    except:
                        pass
                for var in statevars:
                    if var in savefile["bbb"].keys():
                        try:
                            self.setue(var, savefile["bbb"][var][()])
                            self.variables["stored"][var] = savefile["bbb"][var][()]
                        except Exception as e:
                            raise Exception(f"Could not read variable {var}: {e}")
                if self.info["verbose"]:
                    print(
                        "Native UEDGE-style save successfully restored "
                        + "from {}".format(savefname)
                    )
            else:
                raise ValueError(f"Save file {savefname} format not recognized")

            # TODO
            # Implement structure to read and restore auto-detected changes
