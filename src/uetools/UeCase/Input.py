import os
from uetools.UeUtils import Tools

try:
    from uedge import bbb, com, aph, api, svr

    uedge_is_installed = True
except:
    uedge_is_installed = False


class Input:
    """Object providing functions for reading YAML and HDF5 input files

    Methods
    -------
    radialdiff(difffname, **kwargs)
        reads and restores radial diffusivity profiles from difffname
    read(setupfile=None, restore=True, savefile=None, readinput=True,
            restoresave=False)
        Recursively reads the YAML input file setupfile
    readhdf5(fname)
        Recursively reads the setup variables in HDF5 file fname
    userdiff(fname, **kwargs)
        reads and restpre spatial diffusivity profiles from fname
    """

    def __init__(self, case):
        """Initializes Input object and links it to Case object"""
        self.setue = case.setue
        self.getue = case.getue
        self.get = case.get
        self.info = case.info
        self.mutex = case.mutex
        self.variables = case.variables
        self.tools = Tools()
        self.reload = case.reload
        self.tracker = case.tracker
        self.savefuncs = case.savefuncs
        # Makes self.populate available in input files
        self.populate = case.populate

    def readhdf5(self, fname):
        """Reads the UEDGE input deck from setup group of HDF5

        Assumes the HDF5 file was written by UETOOLS and contains
        the necessay groups and variables.

        Arguments
        -----------------
        fname : str
            path to HDF5 file to read the setup from

        Returns
        -------
        nested dict containing setup data (e.g. input file in UETOOLS type)
        """
        from h5py import File, Group
        from os.path import exists

        def recursive_readhdf5(ret, setup, group=[], root=None):
            if len(group) > 0:
                lastgroup=group[0]
                # Conditional to avoid unnecessary nesting
                if ((len(group) > 1) and \
                    (group[-1] == group[-2]) \
                ):
                    group.pop(-1)
                else:
                    for subgroup in group:
                        # Conditional to avoid unnecessary nesting
                        if subgroup not in root.keys():
                            try:
                                ret = ret[subgroup]
                            except:
                                ret[subgroup] = {}
                                ret = ret[subgroup]
            for name, content in setup.items():
                if isinstance(content, Group):
                    recursive_readhdf5(ret, content, group + [name], root)
                else:
                    ret[name] = content[()]

        ret = dict()
        if not exists(fname):
            raise OSError('File "{}" not found!'.format(fname))
        else:
            with File(fname, "r") as f:
                recursive_readhdf5(ret, f["setup"], root=ret)
        self.info["savefile"] = fname
        return ret

    def read(
        self,
        setupfile=None,
        restore=True,
        savefile=None,
        readinput=True,
        restoresave=False,
        **kwargs,
    ):
        """Reads input file and restores data to UEDGE memory

        The working heart of the UETOOLS Case object. Contains the
        YAML input file parser, which defines the YAML input file
        behavior.

        Arguments
        ---------

        Keyword arguments
        -----------------
        setupfile : str (default = None)
            Path to input file to be read. If None, looks for
            {Case.casename}.hdf5 and tries to restore it.
        restore : bool (default = True)
            Switch whether to set UEDGE parameters to the read data
        savefile : str (default = None)
            Path to save file to be read. If None, uses savefile as
            specified in the input file.
        readinput : bool (default = True)
            Switch whether to re-read input data from file. Necessary
            for assign to work, as it re-assigns the input data to UEDGE
            memory without re-reading the input data in order to track
            changes.
        restoresave : bool (default = True)
            Switch telling Case to call Case.load_state and populate
            the UEDGE memory.

        Modifies
        --------
        Case memory : reads input data to Case.varinput["setup"] and
            writes any data specified in input using special variables
        UEDGE memory : sets UEDGE memory to correspond to the input file
            and restores the variables based on the UEDGE state.

        Returns
        -------
        None
        """
        from copy import deepcopy
        from numpy import array
        from h5py import is_hdf5
        from os.path import abspath
        from Forthon import packageobject
        from os.path import exists

        # Extract user-supplied casename and diff_file
        casename = deepcopy(self.info["casename"])
        diff_file = deepcopy(self.info["diffusivity_file"])
        self.casename_set = False
        if self.mutex() is False:
            raise Exception("Case doesn't own UEDGE memory")

        if readinput is True:
            if setupfile is None:
                print("No setup file specified:")
                if self.info["casename"] is None:
                    raise ValueError("    No casename defined: aborting!")
                else:
                    print("    Using casename '{}'".format(self.info["casename"]))
                setupfile = "{}.yaml".format(self.info["casename"])
            if is_hdf5(setupfile):
                self.variables["input"]["setup"] = self.readhdf5(setupfile)
                self.info["restored_from_hdf5"] = True
                self.setue("GridFileName", abspath(setupfile))
                self.setue("isgriduehdf5", 1)
            else:
                try:
                    self.variables["input"]["setup"] = self.tools.readyaml(setupfile)
                except Exception as e:
                    raise ValueError(f"Input file could not be parsed: {e}")
        setup = deepcopy(self.variables["input"]["setup"])
        # Pop out groups that cannot be parsed by default
        if "commands" in setup:
            commands = setup.pop("commands")
        if "detected" in setup:
            detected = setup.pop("detected")

        # TODO: Add mist.dat as an optional parameter/etc to allow changing
        #       the name/path to the data file
        def setinputrecursive(dictobj, group=[]):
            if not isinstance(dictobj, dict):
                # Skip UeCase-unique parameters
                if group[-1] not in self.variables["omit"] + ["chgstate_format"]:
                    # NOTE: Not sure what to do with chgstate_format, fails for some strange reason...
                    # NOTE: Should not be an input, just skip for the time being
                    # Avoid overwriting grid path when restoring from HDF5
                    if group[-1] == "GridFileName":
                        if not self.info["restored_from_hdf5"]:
                            self.setue(
                                "GridFileName",
                                os.path.join(self.info["location"], dictobj),
                            )
                        return
                    # Circumvent the padding with nulls for strings
                    try:
                        dictobj = dictobj.ljust(len(self.getue(group[-2])[group[-1]]))
                    except:
                        pass
                    # Set all other parameters
                    if isinstance(dictobj, list):
                        # Set whole array, given a list
                        if isinstance(group[-1], int):
                            var = group[-2]
                            ind0 = group[-1]
                        else:
                            var = group[-1]
                            ind0 = 0
                        if len(self.getue(var).shape) > 1:
                            #                            print('2D arrays set')
                            ueshape = self.getue(var).shape
                            # TODO: Check compatibility here
                            # Check whether full 2D array is being set
                            if ueshape == array(dictobj).shape:
                                self.setue(var, array(dictobj))
                            elif ueshape == (array(dictobj).T).shape:
                                self.setue(var, array(dictobj).T)
                            # Only part of the array is being set
                            elif len(ueshape) != len(array(dictobj).shape):
                                for j in range(ueshape[-1]):
                                    self.getue(var, cp=False)[:, j].put(
                                        range(ind0, len(dictobj) + ind0), dictobj
                                    )
                            else:
                                raise ValueError(
                                    "!!! ERROR !!! Unable to determine "
                                    "shape of {} from input".format(var)
                                )
                        # A 1D array is being set with a list
                        else:
                            #                            print('1D array set')
                            # Edit array without copying == setue
                            # TODO: Check compatibility here
                            self.getue(var, cp=False).put(
                                range(ind0, len(dictobj) + ind0), dictobj
                            )
                    # A single entry in a 1D or 2D array is set
                    elif isinstance(group[-1], int):
                        #                        print('Setting single index')
                        if len(group) > 2 and isinstance(group[-2], int):
                            # 2D array
                            datalist = self.getue(group[-3])
                            datalist[group[-2], group[-1]] = dictobj
                            # TODO: Check compatibility here
                            self.setue(group[-3], datalist)
                        else:
                            # 1D array
                            datalist = self.getue(group[-2])
                            datalist[group[-1]] = dictobj
                            # TODO: Check compatibility here
                            self.setue(group[-2], datalist)
                    elif dictobj is None:
                        print("WARNING Unset specifier in input:", group[-1])
                    else:
                        try:
                            self.setue(group[-1], dictobj)
                        except KeyError as e:
                            print(f"WARNING Could not set '{group[-1]}' to '{dictobj}'. Reason: {e}")

                else:  # Set calls to restore diffusivities
                    if (group[-1] == "savefile") and (
                        self.info["restored_from_hdf5"] == True
                    ):
                        pass
                    elif group[-1] in ["casename", "commands", "chgstate_format"]:
                        if group[-1]=="casename":
                            self.casename_set = True
                        self.info[group[-1]] = dictobj
                    elif group[-1] in [
                            "userdifffname", 
                            "radialdifffname",
                            "diff_file",
                            "savefile",
                    ]:
                        if isinstance(dictobj, (bytes, bytearray)):
                            dictobj = dictobj.decode("UTF-8")
                        if dictobj is not False:
                            self.info[group[-1]] = "/".join(
                                [self.info["location"], dictobj]
                        )
                    elif group[-1] in self.variables["omit"]:# ["lynix", "lyphix", "lytex", "lytix"]:
                        pass
                    else:
                        if isinstance(dictobj, (bytes, bytearray)):
                            dictobj = dictobj.decode("UTF-8")
                        try:
                            self.info[group[-1]] = os.path.join(self.info["location"], dictobj)
                        except TypeError:
                            # dictobj may be e.g. bool that can't be joined with path
                            self.info[group[-1]] = self.info["location"]
            else:
                for key, value in dictobj.items():
                    dictobj = setinputrecursive(value, group + [key])
                return dictobj

        # Set group order to assure proper allocation and avoid warnings
        if "grid" in setup.keys():
            setinputrecursive(setup.pop("grid"))
        else:
            print(
                "Setup group 'grid' not detected: trying to set terms"
                + " individually from input file"
            )
            for var in [
                "mhdgeo",
                "gengrid",
                "isgriduehdf5",
                "GridFileName",
                "geometry",
                "isnonog",
            ]:
                self.setue(var, self.tools.hdf5search(setupfile, var))
        if "species" in setup.keys():
            setinputrecursive(setup.pop("species"))
        else:
            print(
                "Setup group 'species' not detected: trying to set terms"
                + " individually from input file"
            )
            for var in [
                "ngsp",
                "nhsp",
                "nhgsp",
                "isimpon",
                "nzsp",
                "isupgon",
                "ziin",
                "ismctab",
                "mcfilename",
                "znuclin",
                "nusp_imp",
            ]:
                self.setue(var, self.tools.hdf5search(setupfile, var))
        try:
            self.setue("restart", self.tools.hdf5search(setupfile, "restart"))
        except:
            pass
        packageobject("bbb").getpyobject("allocate")()
        if restore is True:
            setinputrecursive(setup)
            packageobject("bbb").getpyobject("allocate")()
            if "detected" in locals():
                setinputrecursive(detected)
            if isinstance(self.info["casename"], bytes):
                self.info["casename"] = self.info["casename"].decode("UTF-8")
            if isinstance(self.info["savefile"], bytes):
                self.info["savefile"] = self.info["savefile"].decode("UTF-8")
            if self.info["restored_from_hdf5"] is True:

                prfile = setupfile
                if len(prfile.split('/'))>3:
                    prfile = ".../{}".format('/'.join(prfile.split('/')[-3:]))
                print("=================================================")
                print("Restoring case from HDF5 file:")
                print("  Rate dirs read from .uedgerc")
                print("  Grid read from {}".format(prfile))
                self.info["diffusivity_file"] = os.path.join(self.info["location"], setupfile)

            # Override with diff_file maually defined diff_file upon
            if diff_file is not None:
                self.info["diffusivity_file"] = os.path.join(self.info["location"], diff_file)
            # Otherwise, try setting accoridng to input
            else:
                # diff_file takes precedence
                if self.info["diffusivity_file"] is None:
                    # if diff_file not set, override using old diffusion files
                    # userdifffname takes precedence over radialdifffile
                    # if both are present
                    if "radialdifffname" in self.info:
                        if (self.info["radialdifffname"] is not None) and (
                            self.info["radialdifffname"] is not False
                        ):
                            self.info["diffusivity_file"] = \
                                    self.info['radialdifffname']
                        del self.info["radialdifffname"]
                    if "userdifffname" in self.info:
                        if (self.info["userdifffname"] is not None) and (
                            self.info["userdifffname"] is not False
                        ):
                            self.info["diffusivity_file"] = \
                                    self.info["userdifffname"]
                        del self.info["userdifffname"]
            if (self.info["diffusivity_file"] is None) and (
                self.getue("isbohmcalc") in [0, 2]
            ):
                self.info["diffusivity_file"] = self.info["savefile"]

                prfile = self.info["diffusivity_file"]
                if len(prfile.split('/'))>3:
                    prfile = ".../{}".format('/'.join(prfile.split('/')[-3:]))
                print(
                    "No diffusivity-file supplied: reading from "
                    + 'save-file "{}"'.format(prfile)
                )
            # Set diffusivities based on file if model requires profiles
            if self.getue("isbohmcalc") == 0:
                prfile = self.info["diffusivity_file"]
                if len(prfile.split('/'))>3:
                    prfile = ".../{}".format('/'.join(prfile.split('/')[-3:]))
                print(
                    "  User-specified diffusivities read from HDF5 "
                    + 'file "{}"'.format(prfile)
                )
                self.userdiff(self.info["diffusivity_file"])
                try:
                    self.userdiff(self.info["diffusivity_file"])
                except Exception as e:
                    print(
                        f"WARNING: failed to read diffusivities "
                        + "from {}: {}".format(self.info["diffusivity_file"], e)
                    )
            elif self.getue("isbohmcalc") == 2:
                print(
                    "  Radial diffusivities read from HDF5 file "
                    + '"{}"'.format(self.info["diffusivity_file"])
                )
                try:
                    self.radialdiff(self.info["diffusivity_file"])
                except Exception as e:
                    print(
                        f"WARNING: failed to read radial diffusivities "
                        + "from {}: {}".format(self.info["diffusivity_file"], e)
                    )
        if self.casename_set is False:
            self.info["casename"] = casename
        if savefile is not None:
            self.info["savefile"] = savefile
        if self.info["aphdir"] is not None:
            self.setue("aphdir", self.info["aphdir"])
        if self.info["apidir"] is not None:
            self.setue("apidir", self.info["apidir"])
        if restoresave is True:
            if (self.info["savefile"] is None) and (self.get("restart") == 1):
                if self.info["casename"] is not None:
                    for suffix in ["h5","hdf5"]:
                        savefile = ".".join([self.info["casename"], suffix])
                        if exists(savefile):
                            self.info["savefile"] = savefile
                if self.info["savefile"] is None:
                    raise ValueError("No save-file supplied!")
            if exists(self.info["savefile"]):
                self.savefuncs.load_state(self.info["savefile"], **kwargs)
            else:
                raise ValueError("Could not open save-file '{}'".format(\
                    self.info["savefile"]))
        if uedge_is_installed and not self.info["inplace"]:
            self.tracker.get_uevars()
        # NOTE: Get the hashes before running any commands. This way,
        # any changes done in external scripts etc will be registered.
        # This is useful (and necessary) to capture changes to arrays
        # being modified.
        self.reload()
        if not self.info["restored_from_hdf5"]:
            if "commands" in locals():
                for command in commands:
                    try:
                        exec(command)
                    except Exception as e:
                        print(f"Command {command} failed: {e}")

    def userdiff(self, fname, **kwargs):
        """Sets user-defined diffusivities from HDF5 file

        Arguments
        ------------
        diffname : str
            HDF5 file from where to read 'diffusivities/bbb'-values

        Modifies
        --------
        UEDGE memory : modifies user-defined 2D diffusivity arrays
            dif_use, kye_use, kyi_use, and tray_use. Also restores
            coefficient fcdif if in file.

        Returns
        -------
        None
        """
        from h5py import File, is_hdf5
        from os.path import exists

        # TODO: replace with save-group function call?
        # NOTE: not sure why h5pickle throws error here?
        # No matter, we are only reading: use h5py

        if self.mutex() is False:
            raise Exception("Case doesn't own UEDGE memory")

        if not exists(fname):
            if is_hdf5(self.info["filename"]):
                fname = self.info["filename"]
            else:
                raise Exception(f"Diffusivity file '{fname}' not found!")
        with File(fname) as file:
            for variable in ["dif_use", "kye_use", "kyi_use", "tray_use"]:
                self.setue(variable, file[f"diffusivities/bbb/{variable}"][()])
            for variable in ["difni", "kye", "kye", "travis", "fcdif"]:
                if variable in file["diffusivities/bbb"]:
                    self.setue(variable, file[f"diffusivities/bbb/{variable}"][()])

    def radialdiff(self, difffname, **kwargs):
        """Sets radially varying diffusivities from HDF5 file

        Arguments
        ------------
        diffname : str
            HDF5 file from where to read 'diffusivities/bbb'-values

        Modifies
        --------
        UEDGE memory : radial diffusivity coefficients difniv, kyev,
            kyiv, and travisv

        Returns
        -------
        None
        """
        from h5py import File
        from os.path import exists

        if self.mutex() is False:
            raise Exception("Case doesn't own UEDGE memory")

        if not exists(difffname):
            difffname = self.info["filename"]
            if not exists(self.info["filename"]):
                raise Exception("Diffusivity file not found!")

        with File(difffname) as file:
            for variable in ["difniv", "kyev", "kyiv", "travisv"]:
                self.setue(variable, file["diffusivities"]["bbb"][variable][()])

    # def validate_settings(self, var, data):
    # Checks: check variable if:
    #   - Dimension is dependent on ngsp or nisp
    #   - Dimension is equal to ngspmx or nispmx
    # If either satisfied:
    #   - See if index beyond nisp/ngs is modified
    #   - If yes, raise a warning
    # pass
    # if var in self.variables['dims']:
    #    if  ('ngsp' in str(self.variables['dims'][var])) or\
    #        ('nisp' in str(self.variables['dims'][var])):
    #        raise NotImplementedError("")
    #        print('Shaped')
