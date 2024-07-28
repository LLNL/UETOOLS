from .Tools import Tools


class Tracker():
    def __init__(self, case):
        self.getue = case.getue
        self.variables = case.variables
        self.get_bottomkeys = Tools().get_bottomkeys

    def get_uevars(self):
        """ Records the attributes, packages, and hashes of all UEDGE variables

        """
        try:
            from Forthon import package, packageobject
        except:
            pass
        # Create object for storing variables
        self.variables['hashes']['undef'] = {} 
        self.variables['hashes']['input'] = {}
        # Loop through UEDGE packages
        for pkg in package():
            # Get package Object
            pkgobj = packageobject(pkg)
            # Loop through all variables in package
            for var in pkgobj.varlist():
                # Parse the variable attributes into list
                attrs = pkgobj.listvar(var).split('Attributes:')[1].split(\
                    '\n')[0].replace('  ','').split()
                try:
                    self.variables['hashes']['undef'][attrs[0]][var] = [pkg, hash(str(var))]
                except:
                    self.variables['hashes']['undef'][attrs[0]] = {}
                    self.variables['hashes']['undef'][attrs[0]][var] = [pkg, hash(str(var))]
                # All variables are assigned their group: if there is only
                # one attribute, variable is not assigned any other attributes
                # Otherwise, there are custom attributes for the variable
                if len(attrs) > 1:
                    # Loop through all assigned attributes for variable
                    for att in attrs[1:]:
                        # Assert a dict entry for every attribute
                        try:
                            self.variables['hashes'][att]
                        except:
                            self.variables['hashes'][att] = {}
                        # Add the variable to nested dict, with entry of
                        # package name and current hash
                        self.variables['hashes'][att][var] = [pkg, hash(str(\
                            packageobject(pkg).getpyobject(var)))]
            # Move full packages over: backstop for older versions
            for pkg in ['Ynorm', 'Volsrc']:
                if pkg in self.variables['hashes']['undef']:
                    for var, struct in self.variables['hashes']['undef'][pkg].items():
                        self.variables['hashes']['input'][var] = struct
#                    del( self.variables['hashes']['undef'][pkg] )

    def gather_changes(self, vardict, changes=None, key=None):
        try:
            from Forthon import packageobject
        except:
            pass
        """ Recursively checks for changes in all variables in vardict 

        
        """
        # Loop through all dictionary entries
        for key, subdict in vardict.items():
            # If there are nested dicts, traverse recursively down
            if isinstance(subdict, dict):
                changes = self.gather_changes(subdict, changes, key)
            # If dict contains list, we've hit rock bottom
            elif isinstance(subdict, list):
                # Get the current hash of the variable
                checkhash = hash(str(packageobject(\
                    subdict[0]).getpyobject(key)))
                # Check whether the hash differs
                if subdict[1] !=  checkhash:
                        # Store changed var name to list: assert list
                        try:
                            changes.append(key)
                        except:
                            changes = []
                            changes.append(key)
                        # Update the hash to current value
                        vardict[key] = [subdict[0], checkhash]
            # Catch any exceptions
            else:
                print(type(subdict), subdict)
                print(key, ' passed through loop! '
                    'You shouldnt be here, go away!')
        # Pass variables down recursively
        return changes
            
            

    # TODO: avoid purging changed variables between saves?

    def record_changes(self):
        """ Stores changed variables to self.variables['input']['setup']['detected'] """
        # NOTE
        # Conditionals could be introduced to detect whether any given
        # variable is being used given the status of another setting.
        # This could either be done by hard-coding the conditionals
        # in this function, or the settings and flags could be passed
        # as attributes/part of attributes in UEDGE and be parsed
        # by one of the tracking functions.

        # Get list of all variables written to setup block
        setupvars = self.get_bottomkeys(self.variables['input']['setup'])
        # Get a list of all input variables changed since save/read
        changedvars = self.gather_changes(self.variables['hashes']['input']) 
        # Changes detected
        if changedvars is not None:
            # Assert dictionary for detected changed variables exist
            try:
                self.variables['input']['setup']['detected']
            except:
                self.variables['input']['setup']['detected'] = {}
            for var in changedvars:
                # Values are grabbed from UEDGE, but store them
                # here regardless. Carry pointers to save memory
                self.variables['input']['setup']['detected'][var] = \
                    self.getue(var, cp=False)

    def get_detected_changes(self, savefile):
        from h5py import File
        from os.path import exists
        if not exists(savefile):
            raise ValueError(f"File {savefile} does not exist.")
        with File(savefile, 'r') as f:
            try: 
                detected = f['setup/detected']
            except:
                raise KeyError(f"File {savefile} is not a UETOOLS save")
            for key, item in detected.items():
                print(key)

