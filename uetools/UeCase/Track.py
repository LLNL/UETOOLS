class Track:
    def __init__(self, var, package=None):
        """
        Usage: watch(<varname>)
        Create a watchpoint on named variable.
        """
        from Forthon import packageobject
        from numpy import sum

        self.vname = var
        self.pkgname = package
        self.pkgobj = packageobject(self.pkgname)
        try:
            self.varval = self.pkgobj.getpyobject(var)
        except:
            self.pkgname = self.getpackage(self.vname)
            self.varval = self.pkgobj.getpyobject(var)
        self.hash = self.hashvar(self.varval)

    def hashvar(self, var):
        from numpy import sum

        try:
            if len(self.var.shape) > 1:
                if sum(self.var.shape) > 10:
                    return hash(str(sum(var)))
        except:
            pass
        return hash(str(var))

    def changed(self):
        """
        Check current variable hash against watched hash.
        """
        newhash = self.hashvar(self.pkgobj.getpyobject(self.vname))
        if newhash != self.hash:
            return True
        else:
            return False

    def getpackage(self, var):
        """Returns the package name of variable"""
        from Forthon import package

        for packagestr in package():
            if var in packageobject(packagestr).varlist():
                return packagestr


class Tracker:
    def __init__(self):
        from Forthon import package, packageobject

        self.tracking = []
        for packstr in package():
            package = packageobject(packstr)
            for var in package.varlist():
                self.tracking.append(Track(var, packstr))

    def printchanged(self, includepkgs=False):
        for var in self.trackting:
            if var.changed():
                if includepkgs is True:
                    print(var.pname + "." + w.vname)
                else:
                    print(var.vname)

    def getchanged(self, string=False):
        changed = []
        for var in self.tracking:
            if var.changed():
                if string is False:
                    changed.append(var)
                else:
                    changed.append(var.vname)
        return changed


class watchpkg:
    def __init__(self, cpkg):
        self.watching = []
        pkg = ul.packagename2object(cpkg)
        self.pname = pkg.name()
        for v in pkg.varlist():
            self.watching.append(watch(v, pkg=pkg.name()))

    def printchanged(self, includepkgs=None):
        for w in self.watching:
            if w.changed():
                if includepkgs != None:
                    print(w.pname + "." + w.vname)
                else:
                    print(w.vname)

    def changed(self, includepkgs=None):
        changed = []
        for w in self.watching:
            if w.changed():
                if includepkgs != None:
                    changed.append(w.pname + "." + w.vname)
                else:
                    changed.append(w.vname)
        return changed


class waypoint:
    def __init__(self):
        self.var = {}
        for p in ul.packages:
            for var in p.varlist():
                pkg, pvar = findvar(var)
                if pvar is not None:
                    try:
                        self.var[p][var] = hash(str(pvar))
                    except:
                        self.var[p] = {}
                        self.var[p][var] = hash(str(pvar))

    def changed(self):
        changed = []
        for package, pkgvar in self.var.items():
            for var, hsh in pkgvar.items():
                pkg, pvar = findvar(var)
                if pvar is not None:
                    if self.var[package][var] != hash(str(pvar)):
                        changed.append(var)
        return changed
