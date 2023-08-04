

class Tracker():
    def __init__(self, **kwargs): 
        print('PING')
        # NOTE: I cannot get this thing to initialize using super.......

    def get_uevars(self):
        from Forthon import package, packageobject
        self.uevars = {'undef':{}}
        inputs, maybes, other, unclassified= [], [], [], []
        for pkg in package():
            pkgobj = packageobject(pkg)
            for var in pkgobj.varlist():
                attrs = pkgobj.listvar(var).split('Attributes:')[1].split(\
                    '\n')[0].replace('  ','').split()
                if len(attrs) == 1:
                    self.uevars['undef'][var] = hash(str(var))
                else:
                    for att in attrs[1:]:
                        try:
                            self.uevars[att]
                        except:
                            self.uevars[att] = {}
                        self.uevars[att][var] = [pkg, hash(str(\
                            packageobject(pkg).getpyobject(var)))]

    def gather_changes(self, vardict, changes=None, key=None):
        from Forthon import packageobject
        """ Recursively checks for changes in all variables in vardict """
        # NOTE
        # A scheme where old hashes are stored could be implemented.
        # In this case, the comparison could append new hashes to 
        # a list. Then, the current hash could be checked against the
        # first hash, in order to determine whether there has been 
        # any changes from the original. Upon a save, the list of
        # hashes could be reset to their current values
        for key, subdict in vardict.items():
            if isinstance(subdict, dict):
                changes = check_changes(subdict, changes, key)
            elif isinstance(subdict, list):
                checkhash = hash(str(packageobject(\
                    subdict[0]).getpyobject(key)))
                if subdict[1] !=  checkhash:
                        try:
                            changes.append(key)
                        except:
                            changes = []
                            changes.append(key)
                        vardict[key] = [subdict[0], checkhash]
            else:
                print(subdict, ' passed through loop! '
                    'You shouldnt be here, go away!')
            
        return changes
            
            

    def record_changes(self):
        # Get list of all variables written to setup block
        setupvars = self.get_bottomkeys(self.varinput['setup'])
        # Get a list of all input variables changed since save/read
        changedvars = self.gather_changes(self.uevars['input']) 
        # Changes detected
        if changedvars is not None:
            try:
                self.varinput['setup']['detected']
            except:
                self.varinput['setup']['detected'] = {}
            for var in changedvars:
                # Values are grabbed from UEDGE, but store them
                # here regardless. Carry pointers to save memory
                self.varinput['setup']['detected'][var] = \
                    self.getue(var, cp=False)

