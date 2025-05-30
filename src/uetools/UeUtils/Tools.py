class Tools:
    """Class containing useful standalone Tools for UETOOLS

    Methods
    -------
    get_bottomkeys(vardict, keys=None)
        recursively get the last keys in nested dict
    hdf5tree(fname, group=None, depth=None, pre='')
        display HDF5 file contents as a tree structure
    hdf5search(file, var)
        searches HDF5 file for var and returns values if found
    is_case(filename)
        tests whether a file is a valid UETOOLS HDF5 save file
    count_cases(path, omit=['ignore', 'archive'], duplicates=True)
        counts the number of UETOOLS HDF5 save files in path
    smooth_curve(x, y, s=1, **kwargs)
        smooths arbitrary (x,y) curve by s-order spline
    readyaml(fname, **kwargs)
        returns dict of yaml file
    """

    def get_bottomkeys(self, vardict, keys=None):
        """Returns a list of keys at the bottom of nested dict"""
        # Iterate through all dictionary entries
        for key, var in vardict.items():
            # If there is a nested dict, traverse down recursively
            if isinstance(var, dict):
                keys = self.get_bottomkeys(var, keys)
            # Otherwise, ensure that the key is a variable, not index
            elif isinstance(key, str):
                # Assert there is a list to store to, and store variable
                try:
                    keys.append(key)
                except:
                    keys = []
                    keys.append(key)
            # The YAML can define which indices to set: this catches and
            # passes on any such instances
            else:
                continue
        return keys

    def hdf5tree(self, fname, group=None, depth=None, pre=""):
        """Displays HDF5 file contents as a tree structure

        Arguments
        ---------
        fname - path to file to display

        Keyword arguments
        -----------------
        group : str (default = None)
            specific (sub)group in file to print
        dept : int (default = None)
            max recursion depth. If None, traverses to the bottom
        pre : str (default = '')
            prepended string, used for recursive evaluation, do not edit

        Returns
        -------
        None
        """
        from h5py import File, _hl
        from copy import deepcopy

        def typestr(var):
            for vt in ["int", "float", "bytes"]:
                if vt in str(type(var)):
                    return "({})".format(vt)
            return type(var)

        if depth is None:
            depth = 1000

        def h5_tree(val, pre="", maxdepth=1000, depth=0):
            # TODO implement variable depth
            items = len(val)
            depth += 1
            if depth > maxdepth:
                return
            for key, val in val.items():
                items -= 1
                if items == 0:
                    # the last item
                    if type(val) == _hl.group.Group:
                        print("{}└── {}".format(pre, key))
                        h5_tree(val, pre + "    ", maxdepth, deepcopy(depth))
                    else:
                        try:
                            len(val)
                            if len(val.shape) == 1:
                                descr = "({})".format(len(val))
                            else:
                                descr = val.shape
                            print("{}└── {} {}".format(pre, key, descr))
                        except:
                            print("{}└── {} {}".format(pre, key, typestr(val[()])))
                else:
                    if type(val) == _hl.group.Group:
                        print("{}├── {}".format(pre, key))
                        h5_tree(val, pre + "│   ", maxdepth, deepcopy(depth))
                    else:
                        try:
                            len(val)
                            if len(val.shape) == 1:
                                descr = "({})".format(len(val))
                            else:
                                descr = val.shape
                            print("{}├── {} {}".format(pre, key, descr))
                        except:
                            print("{}├── {} {}".format(pre, key, typestr(val[()])))

        if type(fname) == _hl.group.Group:
            h5_tree(fname, maxdepth=depth)
        else:
            with File(fname, "r") as h5file:
                print(h5file)
                if group is not None:
                    print(group)
                    h5_tree(h5file[group], maxdepth=depth)
                else:
                    h5_tree(h5file, maxdepth=depth)

    def hdf5search(self, file, var):
        """Searches file and returns value of var if found, None else"""
        from h5py import File, Group

        ret = None
        # Iterate through all dictionary entries
        if isinstance(file, str):
            with File(file) as f:
                for key, varval in f.items():
                    # If there is a nested dict, traverse down recursively
                    if key == var:
                        return varval[()]
                    elif isinstance(varval, Group):
                        ret = self.hdf5search(varval, var)
                        if ret is not None:
                            return ret
                    # The YAML can define which indices to set: this catches and
                    # passes on any such instances
                    else:
                        continue
        else:
            for key, varval in file.items():
                # If there is a nested dict, traverse down recursively
                if key == var:
                    return varval[()]
                elif isinstance(varval, Group):
                    ret = self.hdf5search(varval, var)
                    if ret is not None:
                        return ret
                # The YAML can define which indices to set: this catches and
                # passes on any such instances
                else:
                    continue
        return ret

    def is_case(self, filename: str) -> bool:
        """
        Returns True if the given file is a valid HDF5 UEDGE case file

        Arguments
        ---------
        filename: string of path to a UEDGE HDF5 case file

        Returns
        -------
        True if valid UETOOLS HDF5 input file, False otherwise
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

    def count_cases(self, path, omit=["ignore", "archive"], duplicates=True):
        """Returns number of UETOOLS HDF5 save files detected in path

        Arguments
        ---------
        path - string with path to top-level dir to start recursive
            search from

        Keyword arguments
        ------------------
        omit : list of strings (default = ['ignore', 'archive'])
            list of directory names to omit in recursive search
        duplicates : bool (default = True)
            whether to count duplicate save names (True) or not (False)

        Returns
        -------
        int of number of UETOOLS HDF5 save files found
        """
        from os import walk
        from os.path import join

        cases = []
        do_omit = False
        for root, dirs, files in walk(path):
            for omitdir in omit:
                if omitdir in root:
                    do_omit = True
                else:
                    do_omit = False
            if do_omit is False:
                for name in files:
                    buff = join(root, name)
                    if self.is_case(buff):
                        cases.append(name)
        if duplicates is False:
            res = []
            [res.append(x) for x in cases if x not in res]
        else:
            res = cases

        return len(res)

    def smooth_curve(self, x, y, s=1, **kwargs):
        """Returns s-order spline smoothed x and y data"""
        from scipy.interpolate import splrep, BSpline

        return x, BSpline(*splrep(x, y, s=s))(x)

    def readyaml(self, fname, **kwargs):
        """Reads a YAML file and returns a nested dict

        Arguments
        ------------
        fname : str
            path to YAML file to be read

        Returns
        ------------
        nested dict containing YAML data
        """
        from yaml import safe_load
        from pathlib import Path

        return safe_load(Path(fname).read_text())
