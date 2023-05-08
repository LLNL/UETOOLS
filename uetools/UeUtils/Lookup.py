
class Lookup:

    def getpackage(self, var, verbose=True, **kwargs):
        """ Returns the package name of variable """
        from Forthon import package, packageobject
        ret = []
        for packagestr in package():
            if var in packageobject(packagestr).varlist():
#                return packagestr
                ret.append(packagestr)
        if len(ret) > 1:
            if verbose is True:
                print('WARNING! {} found in {}!'.format(var, ret))
                print('Using {}!'.format(ret[0]))
        elif len(ret) == 0:
            return None
        return ret[0]

    def getpackobj(self, var, verbose=True, **kwargs):
        """ Returns the package object of variable """
        from Forthon import package, packageobject
        ret = []
        for packagestr in package():
            packobj = packageobject(packagestr)
            if var in packobj.varlist():
#                return packobj
                ret.append(packobj)
        if len(ret) > 1:
            if verbose is True:
                print('WARNING! {} found in {}!'.format(var, ret))
                print('Using {}!'.format(ret[0]))
        elif len(ret) == 0:
            return None
        return ret[0]


