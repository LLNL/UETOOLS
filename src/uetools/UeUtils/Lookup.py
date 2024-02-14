class Lookup:
    def getpackage(self, var, verbose=True, **kwargs):
        """Returns the package name of variable"""
        try:
            from Forthon import package, packageobject
        except:
            pass

        ret = []
        for packagestr in package():
            if var in packageobject(packagestr).varlist():
                #                return packagestr
                ret.append(packagestr)
        if len(ret) > 1:
            if verbose is True:
                print("WARNING! {} found in {}!".format(var, ret))
                print("Using {}!".format(ret[0]))
        elif len(ret) == 0:
            return None
        return ret[0]

    def getpackobj(self, var, verbose=True, **kwargs):
        """Returns the package object of variable"""
        try:
            from Forthon import package, packageobject
        except:
            pass

        ret = []
        for packagestr in package():
            packobj = packageobject(packagestr)
            if var in packobj.varlist():
                #                return packobj
                ret.append(packobj)
        if len(ret) > 1:
            if verbose is True:
                print("WARNING! {} found in {}!".format(var, ret))
                print("Using {}!".format(ret[0]))
        elif len(ret) == 0:
            return None
        return ret[0]

    def about(self, variable):
        """Prints *.v info available for variable"""
        abt = self.infostring(variable)
        if abt is not False:
            print(self.getpackobj(variable).listvar(variable))
        else:
            print('Variable "{}" not found.'.format(variable))

    def infostring(self, variable):
        """Returns a string containing *.v contents of variable"""
        try:
            return self.getpackobj(variable, False).listvar(variable)
        except:
            return False

    def aboutdict(self, variable):
        """Creates dictionary contining information about variable

        Parses the infostring into parts based on the different
        attributes. Used by derived functions
        """
        abt = self.infostring(variable)
        attrs = [
            "Package",
            "Group",
            "Attributes",
            "Dimension",
            "Type",
            "Address",
            "Pyaddress",
            "Unit",
            "Comment",
        ]
        vals = []
        keys = []
        for parameter in attrs:
            try:
                [val, abt] = abt.strip().split("{}:".format(parameter))
                vals.append(val.strip())
                keys.append(parameter)
            except:
                pass
        vals.append(abt.strip())
        ret = {}
        for i in range(len(keys)):
            ret[keys[i]] = vals[i + 1].replace("  ", "")
        try:
            ret['Attributes'] = ret['Attributes'].split()
        except:
            raise 
        return ret

    def aboutparameter(self, variable, parameter):
        "Returns the string of variable corresponding to parameter " ""
        return self.aboutdict(variable)[parameter]

    def package(self, variable):
        "Returns the Package string of variable " ""
        return self.aboutparameter(variable, "Package")

    def group(self, variable):
        "Returns the Group string of variable " ""
        return self.aboutparameter(variable, "Group")

    def attributes(self, variable):
        "Returns the Attributes string of variable " ""
        return self.aboutparameter(variable, "Attributes")

    def dimension(self, variable):
        "Returns the Dimension string of variable " ""
        return self.aboutparameter(variable, "Dimension")

    def type(self, variable):
        "Returns the Type string of variable " ""
        return self.aboutparameter(variable, "Type")

    def address(self, variable):
        "Returns the Address string of variable " ""
        return self.aboutparameter(variable, "Address")

    def pyaddress(self, variable):
        "Returns the Pyaddress string of variable " ""
        return self.aboutparameter(variable, "Pyaddress")

    def unit(self, variable):
        "Returns the Unit string of variable " ""
        return self.aboutparameter(variable, "Unit")

    def comment(self, variable):
        "Returns the Comment string of variable " ""
        return self.aboutparameter(variable, "Comment")

    def searchvarname(self, string):
        """Returns all variables whose name contains string"""
        try:
            from Forthon import package, packageobject
        except:
            pass

        ret = []
        for pack in package():
            for var in packageobject(pack).varlist():
                if string.lower() in var.lower():
                    ret.append(var)
        if len(ret) == 0:
            return None
        else:
            return ret

    def search(self, string):
        """Returns list of variable with string in about

        Looks for the supplied string under Group, Attributes, and
        Comment in the about output for all variables.
        """
        try:
            from Forthon import package, packageobject
        except:
            pass

        ret = []
        for pack in package():
            for var in packageobject(pack).varlist():
                attrs = self.aboutdict(var)
                for attr in ["Group", "Attributes", "Comment"]:
                    if string.lower() in str(attrs[attr]).lower():
                        ret.append(var)
        if len(ret) == 0:
            return None
        else:
            return ret
