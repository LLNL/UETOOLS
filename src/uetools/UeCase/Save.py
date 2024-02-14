class Save:
    def savevar(self, savefile, groups, variable, data, **kwargs):
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

    def savegroup(self, savename, group, append=True, **kwargs):
        # TODO: Do we want to always rewrite or try to append?
        from h5pickle import File

        with File(savename, "a" if append else "w") as savefile:
            if append is False:
                self.savemetadata(savefile)
            self.recursivesave(savefile, self.varinput[group], [group])

    def savesetup(self, savename, **kwargs):
        self.savegroup(savename, "setup", **kwargs)

    def savegrid(self, savename, **kwargs):
        self.savegroup(savename, "grid", **kwargs)

    def savediffusivities(self, savename, **kwargs):
        self.savegroup(savename, "diffusivities")

    def recursivesave(self, savefile, saveobj, group=[], **kwargs):
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
#                        self.savevar(savefile, group, var, self.getue(var))
#                    return
                else:
                    try:
                        value = self.getue(variable)
                    except:
                        for var in saveobj:
                            for var in saveobj:
                                self.savevar(savefile, group, var, self.getue(var))
                            return
                self.savevar(savefile, group, variable, value)
            # Store requested data
            elif isinstance(saveobj, list):
                for variable in saveobj:
                    self.savevar(savefile, group, variable, self.getue(variable))
        # Recursively go deeper in structure
        else:
            for key, value in saveobj.items():
                if isinstance(key, int):
                    # Save the full array once only
                    variable = group.pop(-1)
                    self.savevar(savefile, group, variable, self.getue(variable))
                    return
                saveobj = self.recursivesave(savefile, value, group + [key])
            return saveobj

    def savemetadata(self, savefile, **kwargs):
        from time import time, ctime

        try:
            savefile.attrs["casename"] = self.casename
        except:
            pass
        savefile.attrs["UETOOLS_version"] = self.uetoolsversion
        savefile.attrs["time"] = time()
        savefile.attrs["ctime"] = ctime()
        savefile.attrs["code"] = "UEDGE"
        savefile.attrs["ver"] = self.uedge_ver
        savefile.attrs["pyver"] = self.pyver
        savefile.attrs["user"] = self.user
        savefile.attrs["hostname"] = self.hostname
        try:
            savefile.attrs["location"] = self.location
        except:
            pass

    def save(self, savefname, group=None, append=False, pickle=False, 
            postprocess=True, **kwargs):
        """Saves HDF5 file containing UeCase data

        Arguments
        ---------
        savefname : str
            path to/name of file to write data to

        Keyword arguments
        -----------------
        group : str (default = None)
            group identifier of group to be written to file. If None,
            all data stored in UeCase is written
        """
        from h5py import File
        from Forthon import packageobject

        if self.inplace:
            print("Data read from file, no data to save. Aborting.")
            return
        if postprocess is True:
            bbb = packageobject('bbb')
            bbb.engbal(self.get('pcoree') + self.get('pcorei'))
            bbb.plateflux()
            bbb.wallflux()
        # Check and store any changes since case last saved/read
        self.record_changes()       
        # Open file for writing
        with File(savefname, "a" if append else "w") as savefile:
            # Save metadata of case
            self.savemetadata(savefile)
            # Write requested data to file
            if group is None:
                self.recursivesave(savefile, self.varinput)
            else:
                self.recursivesave(savefile, self.varinput[group], [group])

    def dump(self, savefname, **kwargs):
        """ Dumps all variables to a save file
        
        Calls self.save(**kwargs)
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
            savefile.create_group('datadump')
            for pkg in package():
                for var in packageobject(pkg).varlist():
                    self.savevar(savefile, ['datadump'], var, self.getue(var))
        

    def load_state(self, savefname=None, **kwargs):
        """Restores a saved solution

        Keyword arguments
        -----------------
        savefname : str (default = None)
            HDF5 file to read stored solution from. If None, solution
            is read from UeCase object
        **kwargs
            passed to setgroup
        """
        from h5py import File, is_hdf5
        from os.path import exists

        if self.mutex() is False:
            return

        if savefname is None:
            savefname = "{}.hdf5".format(self.casename)
        if not exists(savefname):
            raise ValueError('Save file {} does not exist!'.format(\
                savefname
            ))
        elif not is_hdf5(savefname):
            raise ValueError('Save file {} is not in HDF5 format!'.format(\
                savefname
            ))
            

        with File(savefname, "r") as savefile:
#        # Try reading new, subdivided save file
#        try:
#            # Don't override user-specified name for case by reading from file
#            if casefname is None:
#                self.casename = savefile.attrs["casename"]
#        except:
#            pass
            try:
                for group, variables in savefile["restore"].items():
                    for variable, value in variables.items():
                        self.setue(variable, value[()])
                        self.vars[variable] = value[()]
            # If not, try reading old-style save file
            except:
                for group, variables in self.varinput["restore"].items():
                    for variable in variables:
                        self.setue(variable, savefile[group][variable][()])
                        self.vars[variable] = savefile[group][variable][()]

            # TODO
            # Implement structure to read and restore auto-detected changes

        if self.verbose:
            print("Saved solution successfully restored from {}".format(savefname))
