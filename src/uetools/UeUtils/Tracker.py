

class Tracker():

    def get_uevars(self):
        """ Records the attributes, packages, and hashes of all UEDGE variables

        """
        try:
            from Forthon import package, packageobject
        except:
            pass
        # Create object for storing variables
        self.uevars = {'undef':{}, 'input': {}}
        # Loop through UEDGE packages
        for pkg in package():
            # Get package Object
            pkgobj = packageobject(pkg)
            # Loop through all variables in package
            for var in pkgobj.varlist():
                # Parse the variable attributes into list
                attrs = pkgobj.listvar(var).split('Attributes:')[1].split(\
                    '\n')[0].replace('  ','').split()
                # All variables are assigned their group: if there is only
                # one attribute, variable is not assigned any other attributes
                if len(attrs) == 1:
                    self.uevars['undef'][var] = [pkg, hash(str(var))]
                # Otherwise, there are custom attributes for the variable
                else:
                    # Loop through all assigned attributes for variable
                    for att in attrs[1:]:
                        # Assert a dict entry for every attribute
                        try:
                            self.uevars[att]
                        except:
                            self.uevars[att] = {}
                        # Add the variable to nested dict, with entry of
                        # package name and current hash
                        self.uevars[att][var] = [pkg, hash(str(\
                            packageobject(pkg).getpyobject(var)))]

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
        """ Stores changed variables to self.varinput['setup']['detected'] """
        # NOTE
        # Conditionals could be introduced to detect whether any given
        # variable is being used given the status of another setting.
        # This could either be done by hard-coding the conditionals
        # in this function, or the settings and flags could be passed
        # as attributes/part of attributes in UEDGE and be parsed
        # by one of the tracking functions.

        # Get list of all variables written to setup block
        setupvars = self.get_bottomkeys(self.varinput['setup'])
        # Get a list of all input variables changed since save/read
        changedvars = self.gather_changes(self.uevars['input']) 
        # Changes detected
        if changedvars is not None:
            # Assert dictionary for detected changed variables exist
            try:
                self.varinput['setup']['detected']
            except:
                self.varinput['setup']['detected'] = {}
            for var in changedvars:
                # Values are grabbed from UEDGE, but store them
                # here regardless. Carry pointers to save memory
                self.varinput['setup']['detected'][var] = \
                    self.getue(var, cp=False)

