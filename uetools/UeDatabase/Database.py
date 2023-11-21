from uetools import Case


class Database:
    # NOTE For some reason, it takes forever to create/restore Database
    # whenever polycollections are included. The pickle size is only
    # 26M, no reason for the long runtime...
    def __init__(
        self,
        database,
        savedbname=None,
        dbidentifier="_UeDB",
        sortvar="ne",
        sortlocation="midplane",
        rerun=False,
        rerun_dir=None,
        meshplot_setup=False,
        readinput=False,
    ):
        """
        
        """
        from os import getcwd 
        from os.path import isfile, isdir, exists

        if not exists(database):
            raise ValueError('Folder {} does not exist!'.format(database))
        self.cwd = getcwd()
        self.cases = {}
        self.dbidentifier = dbidentifier
        self.datbasename = database
        self.rerun = rerun
        self.rerun_dir = rerun_dir
        self.readinput = readinput
        self.create_database(database, not meshplot_setup)
        self.ixmp = self.get("ixmp")[0]
        self.iysptrx = self.get("iysptrx")[0]
        self.ixpt1 = self.get("ixpt1")[0][0]
        self.ixpt2 = self.get("ixpt2")[0][0]
        self.nx = self.get("nx")[0]
        self.ny = self.get("ny")[0]
        # TODO: Store commonly used grid locations
        # TODO: Account for different grids

        # Sort here
        try:
            self.sort(sortvar, sortlocation)
        except Exception as e:
            print(f"Error while sorting: {e}")

        # Save if requested
        if savedbname is not None:
            self.save(savedbname)

    def is_case(self, filename: str) -> bool:
        """
        Returns True if the given file is a valid HDF5 UEDGE case file

        Arguments
        ---------
        filename: string
            Path to a UEDGE HDF5 case file

        """
        from h5py import File

        try:
            with File(filename, "r") as f:
                # Check that necessary groups are present in file
                for entry in ["centered", "staggered", "restore", "grid", "setup"]:
                    if entry not in f.keys():
                        return False
            return True
        except:
            return False

    def save(self, savename: str):
        """ Save the database to a pickle file"""
        from pickle import dump

        with open("{}.db".format(savename.split(".")[0]), "wb") as f:
            dump(self, f)

    def sort(self, variable, location, increasing=True, varname=''):
        """ Sorts cases according to variable at location

        Note: This works because python dict preserves the order
        of insertion (since python 3.7).
        """
        from numpy import argsort, where, ndarray

        if isinstance(variable, str):
            self.sortvar = variable
            self.sortlocation = location

            # Make sort location more advance
            if self.sortlocation == "midplane":
                self.sortlocation = (self.ixmp, self.iysptrx + 1)
            elif isinstance(self.sortlocation, int):
                self.sortlocation = [self.sortlocation]
            elif isinstance(self.sortlocation, str):
                print(
                    'Sort location option "{}" not recognized. Aborting'.format(
                        self.sortlocation
                    )
                )
            order = self.get(self.sortvar)
            # Gets the correct indices
            if len(order.shape)>1:
                for ind in self.sortlocation:
                    order = order[:, ind]
            order = argsort(order)
        elif isinstance(variable, (list, ndarray)):
            self.sortvar = varname
            order = argsort(variable)

        if increasing is False:
            order = order[::-1]
        neworder = {}
        for i in order:
            key = list(self.cases.keys())[i]
            neworder[key] = self.cases[key]
        self.cases = neworder

        if isinstance(variable, str):
            self.scanvar = self.get(self.sortvar)
            if self.sortlocation is not None:
                for ind in self.sortlocation:
                    self.scanvar = self.scanvar[:, ind]
        elif isinstance(variable, (list, ndarray)):
            self.scanvar = variable[order] 


    def create_database(self, path, database):
        from  os.path import join, exists
        from os import walk
        from pathlib import Path
        from tqdm import tqdm

        createdb = []
        for parent, dirs, files in walk(path):
            if "ignore" in parent:
                continue
            subdir = parent.split("/")[-1]
            databases = [db for db in files if self.is_case("{}/{}".format(parent, db))]
            if self.rerun is False:
                # HDF5s identified
                if len(databases) == 1:
                    # Verify all necessary data is present
                    if len(subdir) == 0:
                        casekey = databases[0].replace(".hdf5", "")
                    else:
                        casekey = "{}/{}".format(subdir, databases[0].replace(\
                            ".hdf5", ""))
                    self.cases[casekey] = Case(
                        "{}/{}".format(parent, databases[0]),
                        inplace=True,
                        verbose=False,
                        database=database,
                    )
                elif len(databases) > 1: 
                    for db in databases:
                        if len(subdir) == 0:
                            casekey = db.replace(".hdf5", "")
                        else:
                            casekey = "{}/{}".format(subdir, db.replace(\
                                ".hdf5", ""))
                        self.cases[casekey] = Case(
                            "{}/{}".format(parent, db), 
                            inplace=True, 
                            database=database,
                            verbose=False
                        )
                # No database found, store location where input is
                # Changing dirs while executing walk breaks the call
                elif ("input.yaml" in files) and self.readinput:
                    createdb.append(parent)
            else:
                # Tries to restore cases from file
                createdb.append((parent, databases))
#                print(parent, databases) 
#                if "input.yaml" in files:
#                    createdb.append(parent)
        # Now, create and read the files
        if len(createdb) > 0:
            print("===== CREATING NEW CASE DUMPS =====")
            for dbfolder in createdb:
                if self.rerun_dir is not None:
                    subdir = dbfolder[0].split('/')
                    _ = subdir.pop(0)
                    subdir = '/'.join(subdir)                    
                    newdbfolder = join(self.rerun_dir, subdir)
                else:
                    newdbfolder = dbfolder[0]
                    subdir = newdbfolder.split("/")[-1]

                for case in dbfolder[1]:
                    scase = case.replace('.hdf5','')
                    identifier = f"{subdir}/{scase}{self.dbidentifier}"
                    hdf5_file = join(
                        newdbfolder, 
                        f"{scase}.hdf5"
                    )
                    # Re-open HDF5 file
                    savefolder = '/'.join(hdf5_file.split('/')[:-1])
                    if not exists(savefolder):
                        Path(savefolder).mkdir(\
                            parents=True, exist_ok=True
                        )
                    Case(join(dbfolder[0], case), 
                        verbose=False).save(hdf5_file)
                    self.cases[identifier] = Case(hdf5_file, inplace=True, verbose=False)

    def get(self, var, **kwargs):
        """Returns an array of var for all cases

        Note: This can fail if the variable is not the same size for
             every case e.g. different grids.

        """
        # TODO: Supress verbose output
        from numpy import array

        ret = []
        for key, case in self.cases.items():
            ret.append(case.get(var, **kwargs))
        return array(ret)

    def getcase(self, index):
        try:
            return self.cases[list(self.cases.keys())[index]]
        except:
            print("Case #{} does not exist".format(index))
            return

    def get_closest_key(self, target, var, index=None):
        from numpy import ndarray, where
        if isinstance(var, str):
            vararr = self.get(var)
        elif isinstance(var, (list, ndarray)):
            vararr = var
        if (index is None) and (len(vararr.shape) >1):
            raise KeyError('Must define index for multi-species array')
        elif isinstance(index, int):
            vararr = vararr[:, index]
        elif isinstance(index, (tuple, ndarray, list)):
            for i in index:
                vararr = vararr[:, i]
        icase = where( abs(vararr - target) == abs(vararr - target).min())[0][0]
        return list(self.cases.keys())[icase]
        

    def get_closest(self, target, var, index=None):
        return self.cases[self.get_closest_key(target, var, index)]
        


    def teITsep(self, **kwargs):
        return self.plotITsep(self.get("te") / 1.602e-19, **kwargs)

    def teOTsep(self, **kwargs):
        return self.plotOTsep(self.get("te") / 1.602e-19, **kwargs)

    def plotITsep(self, var, **kwargs):
        return self.plotscan(var, (1, self.iysptrx + 1), **kwargs)

    def plotOTsep(self, var, **kwargs):
        return self.plotscan(var, (-2, self.iysptrx + 1), **kwargs)

    def plotOTmax(self, var, inds=(None, None), **kwargs):
        from numpy import max

        return self.plotvar(
            self.scanvar, max(var[:, -2, slice(*inds)], axis=1), **kwargs
        )

    def plotOMP(self, var, **kwargs):
        return self.plotscan(var, (self.ixmp, self.iysptrx + 1), **kwargs)

    def plotscan(self, var, location=(), **kwargs):
        for index in location:
            var = var[:, index]
        return self.plotvar(self.scanvar, var, **kwargs)

    def plotvar(self, xvar, yvar, ax=None, **kwargs):
        """Plots yvar as a function of xvar for all cases"""
        from matplotlib.pyplot import subplots, Figure

        if ax is None:
            f, ax = subplots()
        elif ax is Figure:
            ax = f.get_axes()[0]

        try:
            kwargs["marker"]
        except:
            kwargs["marker"] = "."
        try:
            kwargs["linestyle"]
        except:
            kwargs["linestyle"] = ""
        try:
            kwargs["color"]
        except:
            kwargs["color"] = "k"

        ax.plot(xvar, yvar, **kwargs)
        ax.set_xlabel(self.sortvar)

        return ax.get_figure()

    def ne_2Dseries(self, **kwargs):
        return self.plot_2Dseries(self.get("ne"), **kwargs)

    def te_2Dseries(self, **kwargs):
        self.plot_2Dseries(self.get("te") / 1.602e-19, **kwargs)

    def plot_2Dseries(self, vararray, **kwargs):
        """Returns a series of figures to scroll through"""
        from matplotlib.pyplot import subplots, ion, ioff
        from matplotlib.widgets import Slider, RangeSlider

        ioff()
        f, ax = subplots(figsize=(7, 8))

        try:
            kwargs["zrange"]
            origrange = kwargs["zrange"]
        except:
            kwargs["zrange"] = (
                vararray[1:-1, 1:-1, :].min(),
                vararray[1:-1, 1:-1, :].max(),
            )
            origrange = kwargs["zrange"]

        c = self.getcase(0)
        cbar, verts = c.plotmesh(
            vararray[0], ax=ax, watermark=False, interactive=True, **kwargs
        )
        f.axes[0].set_position([0.125, 0.13, 0.55, 0.85])
        f.axes[1].set_position([0.7, 0.13, 0.82, 0.85])
        slice_position = f.add_axes([0.1, 0.02, 0.65, 0.04])
        slice_slider = Slider(
            slice_position,
            self.sortvar,
            self.scanvar.min(),
            self.scanvar.max(),
            valstep=self.scanvar,
        )
        zrange_position = f.add_axes([0.85, 0.13, 0.04, 0.85])
        zrange_slider = RangeSlider(
            zrange_position,
            "",
            vararray[1:-1, 1:-1, :].min(),
            vararray[1:-1, 1:-1, :].max(),
            valinit=(origrange),
            orientation="vertical",
        )

        def update(val):
            from numpy import where

            slce = slice_slider.val
            index = where(self.scanvar == slce)[0][0]
            verts.set_array(
                vararray[index, 1:-1, 1:-1].reshape(
                    self.getcase(index).get("nx") * self.getcase(index).get("ny")
                )
            )
            verts.set_clim(zrange_slider.val)

        slice_slider.on_changed(update)
        zrange_slider.on_changed(update)

        f.show()
        ion()
        return f, slice_slider, zrange_slider

    def animation(self):
        """Creates an animation from a series of figures"""
        print("TBD")


def restoredb(restorename):
    """Restore a database"""
    from pickle import load

    with open(restorename, "rb") as f:
        db = load(f)
    return db
