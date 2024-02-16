from uetools import Case
from uetools.UeDatabase.DB_1Dplots import DB_1DPlots
from uetools.UeDatabase.DB_2Dplots import DB_2DPlots

def rewrite_case(restorefile, savefile):
    """ Process target function """
    from uetools import Case
    c=Case(restorefile, verbose=False)
    c.save(savefile)
    del c


class Database(DB_1DPlots, DB_2DPlots):
    def __init__(
        self,
        database,
        dbidentifier="_UeDB",
        sortvar="ne",
        sortspecies=None,
        sortlocation="OMPsep",
        rerun=False,
        rerun_dir="UeDB_rerun",
    ):
        """
        
        """
        from os import getcwd 
        from os.path import isfile, isdir, exists

        self.cwd = getcwd()
        self.cases = {}
        self.dbidentifier = dbidentifier
        self.datbasename = database
        self.rerun = rerun
        self.rerun_dir = rerun_dir
        self.create_database(database)
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
            self.sort(sortvar, sortlocation, species=sortspecies)
        except Exception as e:
            print(f"Error while sorting: {e}")

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

    def concatenate(self, database):
        """ Absorbs the cases of another database into this one """
        if isinstance(database, list):
            for db in database:
                for key, case in db.cases.items():
                    self.cases[key] = case
        else:
            for key, case in database.cases.items():
                self.cases[key] = case
        # Re-sort cases based on original sorting parameters
        self.sort(self.sortvar, self.sortlocation, self.increasesort, 
            self.sortspecies, self.origvarname)


    def sort(self, variable, location="OMPsep", increasing=True, species=None, varname=None):
        """ Sorts cases according to variable at location

        Note: This works because python dict preserves the order
        of insertion (since python 3.7).
        """
        from numpy import argsort, where, ndarray, array
        from copy import deepcopy

        self.increasesort = deepcopy(increasing)
        self.origvarname = deepcopy(varname)
        self.sortspecies = species
        # Request sort based on case values
        if isinstance(variable, str):
            self.sortvar = variable
            self.sortlocation = location
            # Check for valid keyword locations
            if self.sortlocation in ["OMPsep", "OTsep", "ITsep", 
                    "ITmax", "OTmax", "ITmin", "OTmin"]:
                pass
            # Assert sort location to be list for later use
            elif isinstance(self.sortlocation, int):
                self.sortlocation = [self.sortlocation]
            # Raise error if not OK
            elif isinstance(self.sortlocation, str):
                raise ValueError(
                    'Sort location option "{}" not recognized. Aborting'.format(
                        self.sortlocation
                    )
                )
            # Get the array to sort
            order = self.get(self.sortvar)
            # 2D arrays: sort by location
            if len(order.shape)>1:
                # Keyword sorting location
                if isinstance(self.sortlocation, str):
                    order_index = []
                    # Loop through each case
                    for i in range(len(order)):
                        # Get helpers
                        case = self.getcase(i)
                        ixmp = case.get('ixmp')
                        iysptrx = case.get('iysptrx')
                        # For each case, access keyword location value
                        if self.sortlocation == "OMPsep":
                            order_index.append(order[i, ixmp, iysptrx+1])
                        elif self.sortlocation == "OTsep":
                            order_index.append(order[i, -2, iysptrx+1])
                        elif self.sortlocation == "OTmax":
                            order_index.append(max(order[i, -2]))
                        elif self.sortlocation == "OTmin":
                            order_index.append(min(order[i, -2]))
                        elif self.sortlocation == "ITsep":
                            order_index.append(order[i, 1, iysptrx+1])
                        elif self.sortlocation == "ITmax":
                            order_index.append(max(order[i, 1]))
                        elif self.sortlocation == "ITmin":
                            order_index.append(min(order[i, 1]))
                    order_index = array(order_index)
                    # Get the correct species index
                    if len(order_index.shape) == 2:
                        order_index = order_index[:,species]
                    # Get ordered index list
                    order = argsort(order_index)
                    # Get ordered values
                    order_index = [order_index[x] for x in order]
                # Sorting by given location
                else:
                    # Get non-keyword location values
                    for ind in self.sortlocation:
                        order = order[:, ind]
                    # Get the correct species index
                    if len(order.shape) == 2:
                        order = order[:,species]
                    order_index = deepcopy(order)
                    order = argsort(order)
                    order_index = [order_index[x] for x in order]
            # 1D-array to sort
            else:
                order_index = deepcopy(order)
                order = argsort(order)
                order_index = [order_index[x] for x in order]
                
        # Sort by supplied list of values
        elif isinstance(variable, (list, ndarray)):
            if varname is None:
                self.sortlabel = ""
            else:
                self.sortlabel = varname
            order = argsort(variable)
            order_index = [variable[x] for x in order]
        # Sort in decreasing order if increasing is False
        if increasing is False:
            order = order[::-1]
            order_index = order_index[::-1]
        # Re-order the case dictionary
        neworder = {}
        for i in order:
            key = list(self.cases.keys())[i]
            neworder[key] = self.cases[key]
        self.cases = neworder
        self.sortvalues = array(order_index)
        self.sortlabel = "{} {}".format(self.sortlocation, self.sortvar)


    def closedb(self):
        for key, case in self.cases.items():
            case.hdf5case.close()

    def create_database(self, path):#, database):
        from  os.path import join, exists, isdir, isfile
        from os import walk
        from pathlib import Path
        from multiprocessing import Process


        createdb = []
        databases = {}
        filelist = []
        # Gather all files to be stored in database based on path
        if isinstance(path, list):
            # List of files to be included in database
            filelist = list(path)
        elif isdir(path):
            # Path is directory, recursively get all files
            for parent, dirs, files in walk(path):
                if "ignore" in parent:
                    continue
                casefiles = [db for db in files if self.is_case( \
                                "{}/{}".format(parent, db))]
                if len(files) > 0:
                    databases[parent] = casefiles
            for parent, files in databases.items():
                for file in files:
                    if parent == '.':
                        filelist.append(file)
                    else:
                        filelist.append('/'.join([parent, file]))
        elif isfile(path):
            with open(path) as file:
                for line in file:
                    plainline = line.split("#")[0].strip()
                    if len(plainline) > 0:
                        filelist.append(plainline)
        else:
            raise ValueError('Folder/file "{}" does not exist!'.format(path))
        if self.rerun is False:
            for file in filelist:
                # NOTE: This should never happen, as any yamls are removed
                # by is_case function
                try:
                    self.cases[file.replace(".hdf5", "")] = Case(
                        file,
                        inplace=True,
                        verbose=False,
                    )
                except Exception as e:
                    print(f"Failed reading case {file}: {e}")
        else:
            # Tries to restore cases from file
            for file in filelist:
                createdb.append(file)
        # Now, create and read any files not created
        if len(createdb) > 0:
            cases = []
            dbcases = []
            print("===== CREATING NEW CASE DUMPS =====")
            for file in createdb:
                if self.rerun_dir is not None:
                    subdir = '/'.join((file.split('/'))[1:-1])
                    newdbfolder = join(self.rerun_dir, subdir)
                else:
                    newdbfolder = '/'.join(file.split('/')[:-1])
                    subdir = newdbfolder.split("/")[-1]
                if len(subdir) == 0:
                    subdir = '.'
                if len(newdbfolder) == 0:
                    newdbfolder = '.'
                scase = file.split('/')[-1].replace(".hdf5", "")
                identifier = f"{subdir}/{scase}{self.dbidentifier}"
                hdf5_file = join(
                    newdbfolder, 
                    f"{scase}{self.dbidentifier}.hdf5"
                )
                savefolder = '/'.join(hdf5_file.split('/')[:-1])
                Path(savefolder).mkdir(\
                    parents=True, exist_ok=True
                )
                cases.append(Process(
                        target=rewrite_case, 
                        args=(file, hdf5_file), 
                        kwargs=()
                ))
                cases[-1].start()
                dbcases.append((identifier, hdf5_file))
            for case in cases:
                case.join()
            for dbcase in dbcases:
                self.cases[dbcase[0]] = Case(
                        dbcase[1], 
                        inplace=True, 
                        verbose=False
                )

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

    def get_closest_key(self, target, var=None, index=None, species=None):
        from numpy import ndarray, where
        if var is None:
            vararr = self.sortvalues
        else:
            if isinstance(var, str):
                vararr = self.get(var)
            elif isinstance(var, (list, ndarray)):
                vararr = var
            if len(vararr.shape)==4:
                if species is None:
                    raise ValueError('Must define species to use for multi-species array')
                else:
                    vararr = vararr[:,:,:,species]
            if (index is None) and (len(vararr.shape) >1):
                raise KeyError('Must define index for multi-dimensional array')
            elif isinstance(index, int):
                vararr = vararr[:, index]
            elif isinstance(index, (tuple, ndarray, list)):
                for i in index:
                    vararr = vararr[:, i]
        icase = where( abs(vararr - target) == abs(vararr - target).min())[0][0]
        return list(self.cases.keys())[icase]
        

    def get_closest(self, target, var=None, index=None, species=None):
        return self.cases[self.get_closest_key(target, var, index, species)]
        

