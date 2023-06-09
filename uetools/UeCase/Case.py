from Forthon import packageobject
from .CasePlot import Caseplot
from .Solver import Solver
from .Track import Tracker
from .Save import Save
from uetools.UeDashboard import Dashboard
from uetools.UeUtils.Lookup import Lookup
from uetools.UeUtils.ConvergeStep import ConvergeStep
from uetools.UePostproc.Postproc import PostProcessors
from uetools.UePostproc.ADAS import ADAS
from uetools.UeConfig.Config import Config
from uedge import bbb, com, aph, api, svr, __version__

# TODO: Where/how to define impurity files and paths?
# TODO: Update data after reading save file
# TODO: make yaml read/write/get/set case-insensitive
# TODO: Consider compression of data
# TODO: implement divergence plotting/calculation
# TODO: Unify all data to be stored in the same dictionary?


class Case(
    Caseplot, Solver, Lookup, PostProcessors, ConvergeStep, Save, ADAS, Dashboard
):
    """UEDGE Case container object.

    Subclasses
    ------------
    Caseplot -- contains plotting routines and widgets
    Solver -- contains time-stepping and concergence routines

    Attributes
    ------------
    allocate : function
        allocates the UEDGE arrays based on the input
    casefname : string
        path to input file where data is read from
    inplace : boolean
        switch defining whether to read data into memory or read from
        HDF5
    vars : nested dict
        variables available based on YAML variable file:
        structure vars[package][varname]
    varinput : nested dict
        variables from YAML variable file:
        structure varinput[group][package][varname]
    packagelist : dict
        lookup dictionary listing packages associated with
        individual variables
    setup : nested dict
        variables defined in YAML input file
    userdifffname : string
        path to user-defined spatial diffusivity coefficients
    radialdifffname : string
        path to user-defined radial diffusivity coefficients
    hdf5case : h5py File object
        HDF5 file reader to read data from

    Methods
    ------------
    set_userdiff():
        sets user-defined diffusion coefficients
    save_userdiff(savefname, readvars=False):
        saves user-defined diffusivity coefficients
    restoreinput():
        sets UEDGE input parameters according to setup attribute
    readsetup(setupfile, restore=True):
        reads a UEDGE input file
    openhdf5(fname, operation):
        returns an open h5Py File object
    readyaml(fname):
        returns a nested dict read from standard YAML file
    closehdf5():
        closes hdf5case File object
    reinitializeplot():
        updates the plotting-related grid parameters
    reload(group=None):
         reads the UEDGE parameters of variables into the UeCase object
    getue(package, variable):
        returns a copy of the UEDGE variable from package
    setue(package, variable, value):
        sets the UEDGE variable in package to data
    get(variable, s=None):
        returns the data stored for variable in UeCase's var
    set_radialdiff():
        sets the radial diffusivites
    save(savefname, group):
        saves UeCase variables to HDF5 file
    savevar(savefile, groups, variable, data):
        writes data and metadata to HDF5 file
    load(casefname=None):
        reads data from HDF5 file to UeCase vars
    readgridhdf5(savefname):
        reads gridue data from HDF5 file
    readgridue(gridue='gridue'):
        reads gridue file into UEDGE and sets the case up to run on
        a gridue grid
    restore(savefname, **kwargs)
        restores a UEDGE solution from HDF5
    """

    def __init__(
        self,
        casefname=None,
        inplace=False,
        variableyamlfile=None,
        casename=None,
        assign=True,
        verbose=True,
        database=False,
        **kwargs
    ):
        """Initializes the UeCase object.

        Keyword arguments
        ------------
        casefname : str (default =  None)
            HDF5 file where to read data from. If None, data is read
            from UEDGE
        inplace : bool (default = False)
            Switch whether to read data from file into UeCase memory
            (False) or get data using file I/O at every call
        variableyamlfile : str (default = None)
            Path to YAML file containing definitions of data and
            variables to be read. If None, accesses the module defaults
            # TODO: Use .uedgerc to define the file in question?
        assign : boolean (default = True)
            Switch whether to assign the current run to the caseobject
        """
        import uetools
        from os.path import exists, abspath
        from os import getlogin, getcwd
        from socket import gethostname
        from matplotlib.pyplot import ioff, ion

        conf = Config(verbose=verbose)
        if conf.configured is False:
            return
        # TODO: add label attribute
        # Read from uedge.label and strip
        # Use to restore and/or save files
        # Initialize class attributes
        # Checksum whether to update data

        # Initialize parameters
        self.casename = casename
        self.verbose = verbose
        self.restored_from_hdf5 = False
        self.uetoolsversion = "1.0"  # UEtools version
        self.allocate = packageobject("bbb").getpyobject("allocate")
        self.casefname = casefname
        self.inplace = inplace
        self.userdifffname = None
        self.radialdifffname = None
        self.hdf5case = None
        self.database = database
        self.pyver = __version__
        self.uedge_ver = (
            packageobject("bbb").getpyobject("uedge_ver")[0].strip().decode("UTF-8")
        )
        self.user = getlogin()
        self.hostname = gethostname()
        self.unset_variables = []
        self.omitvars = [
            "userdifffname",
            "radialdifffname",
            "casename",
            "commands",
            "savefile",
        ]

        self.vars = dict()
        self.varinput = dict()
        self.packagelist = dict()
        # TODO: add hostname, aphdir, apidir, mcfilename, aphfname

        """ Set up structure for reading/writing data """
        # Load all data to object in memory
        if self.inplace is False:
            self.get = self.get_memory
            self.getue = self.getue_memory
            self.set = self.set_memory
            self.setue = self.setue_memory
            try:
                self.location = abspath(self.casefname)
            except:
                self.location = getcwd()
            self.session_id = self.getue("max_session_id") + 1
            setattr(
                packageobject("bbb"), "max_session_id", self.getue("max_session_id") + 1
            )
            self.exmain_evals = self.getue("exmain_evals")
            if assign is True:
                self.assign()
            # Read YAML to get variables to be read/stored/used
            if variableyamlfile is None:  # No YAML variable file requested
                # Use default: find package location and package YAMLs
                variableyamlpath = "{}/{}".format(
                    uetools.__path__[0], "yamls/standardvariables.yml"
                )
            else:  # YAML specified, use user input
                variableyamlpath = "{}/{}".format(path, variableyamlfile)
            self.varinput = self.readyaml(variableyamlpath)  # Read to memory
            if self.casefname is not None:
                self.restore(self.casefname)
            else:
                #                self.populate(verbose=False)
                self.reload()
        # Read all data directly from HDF5 file
        else:
            self.get = self.get_inplace
            self.set = self.getsetue_inplace
            self.getue = self.getsetue_inplace
            self.setue = self.getsetue_inplace
            if self.casefname is None:
                print("Must specify data file when inplace=True! Aborting.")
                return
            self.casefname = abspath(self.casefname)
            if exists(self.casefname):
                try:
                    self.hdf5case = self.openhdf5(self.casefname, "r")
                except:
                    print("Unable to open {}. Aborting.".format(self.casefname))
                    return
            self.load_inplace()

        # Initialize parent classes
        super(Case, self).__init__()
        return

    # NOTE: Update class data, or try reading from forthon first??
    def update(self, **kwargs):
        """Checks if UEDGE state has changed and updates if needed"""

        if self.exmain_evals != self.getue("exmain_evals"):
            self.exmain_evals = self.getue("exmain_evals")
            if self.mutex() is False:
                return
            self.reload()
            self.vertices = self.createpolycollection(self.get("rm"), self.get("zm"))

    def get_inplace(self, variable, s=None, **kwargs):
        """Returns variable from HDF5"""
        from numpy import ndarray

        try:
            retvar = self.hdf5case[self.vars[variable]][()]
        except:
            print("{} not found in {}".format(variable, self.casefname))
            return

        if isinstance(retvar, (ndarray, list)):
            if len(retvar.shape) == 3:
                if s is not None:
                    retvar = retvar[:, :, s]
        return retvar

    def set_memory(self, variable, **kwargs):
        """Returns pointer to variable that can be modified"""
        return self.getue_memory(variable, cp=False, **kwargs)

    def get_memory(self, variable, s=None, **kwargs):
        """Returns variable

        Method assumes unique variable names across all packages and
        groups.

        First, checks if UEDGE solution has changed from previous step.
        If it has, updates UeCase data. Then, looks for data in UeCase
        and returns it if found. If not found, returns data from Forthon
        memory.
        NOTE: return None, or access value from UEDGE?

        Arguments
        ------------
        variable : str
            name of variable to be returned

        Keyword arguments
        ------------
        s : int (default = None)
            species index to be returned for 3-dimensional arrays. If
            None, get returns the full 3-dimensional array. Otherwise,
            get returns the s:th index of the third dimension

        Returns
        ------------
        UeCase value of variable (array/int/str/float)
        """
        # TODO: serach input too? Store input to Vars?
        from numpy import ndarray

        self.update()  # Update results from UEDGE if they have changed
        # Switch to asses where to access data from
        try:
            retvar = self.vars[variable]
        except:
            retvar = self.getue(variable)
        # Check the size of the array, and return index if multi-species array
        if isinstance(retvar, (ndarray, list)):
            if len(retvar.shape) == 3:
                if s is not None:
                    retvar = retvar[:, :, s]
        return retvar

    def assign(self, **kwargs):
        """Assigns the UEDGE session to this object"""
        setattr(packageobject("bbb"), "session_id", self.session_id)
        try:
            # Restore input to UEDGE
            # NOTE: variables not set maintain previous values. Reset
            # all values before setting input?
            self.setinput(readinput=False)
        except:
            pass
        try:
            if self.restored_from_hdf5 is True:
                packageobject("grd").getpyobject("readgrid")(
                    self.getue("GridFileName"), self.vars["runid"].strip()
                )
        except:
            pass

    def getsetue_inplace(self, *args, **kwargs):
        """Placeholder to avoid getting/setting when reading inplace"""
        print("Cannot set/get UEDGE values when reading from HDF5 file")
        return

    def setue_memory(self, variable, value, idx=None, **kwargs):
        """Sets the Forthon variable in package to data

        Arguments
        ------------
        variable : str
            variable name
        value : array/list/float/int
            data to be written to the UEDGE variable
        """
        try:
            package = self.packagelist[variable]
        except:
            package = self.getpackage(variable, verbose=False)
        if self.mutex():
            try:
                setattr(packageobject(package), variable, value)
            except:
                raise KeyError("{} could not be set".format(variable))

    def getue_memory(self, variable, s=None, cp=True, **kwargs):
        """Retireves data from UEDGE variable in package

        Arguments
        ------------
        variable : str
            variable name

        Returns
        ------------
        value of UEDGE variable (array/str/int/float)
        """
        from copy import deepcopy
        from numpy import ndarray

        try:
            package = self.packagelist[variable]
        except:
            package = self.getpackage(variable, verbose=False)

        try:
            if cp is True:
                retvar = deepcopy(packageobject(package).getpyobject(variable))
            else:
                retvar = packageobject(package).getpyobject(variable)
        except:
            raise KeyError("{} not found".format(variable))

        if isinstance(retvar, (ndarray, list)):
            if len(retvar.shape) == 3:
                if s is not None:
                    retvar = retvar[:, :, s]
        return retvar

    def reload(self, group=None, **kwargs):
        """Reloads variables from UEDGE to UeCase

        Omits the setup file to avoid overwriting original
        setup.

        TODO: Is this desired behavior?

        Keyword arguments
        ------------
        group : str (default = None)
            group specifier to reload. If None, reloads all
            variables in vars
        """
        from copy import deepcopy
        from numpy import ndarray, int64, float64

        if self.mutex is False:
            return
        try:
            commands = self.varinput["setup"].pop("commands")
        except:
            pass

        def recursivereload(dictobj=None, group=[]):
            if dictobj is None:
                dictobj = self.varinput
            if not isinstance(dictobj, dict):
                # Reached bottom of nested dictionaries: determine format
                if isinstance(dictobj, (list, ndarray)):
                    # We have a list: either list of variables to store or
                    # list defning the variable array
                    if self.getpackage(group[-1], verbose=False) != None:
                        # Request to set array starting from index 0:
                        # just read the variable into memory
                        self.vars[group[-1]] = self.getue(group[-1])
                    # TODO: check/verify this
                    elif isinstance(group[-1], int):
                        # Setting subarray, store variable
                        self.vars[group[-2]] = self.getue(group[-2])
                    else:
                        # List of variables, store each
                        for variable in dictobj:
                            self.vars[variable] = self.getue(variable)
                elif isinstance(group[-1], int):
                    self.vars[group[-2]] = self.getue(group[-2])
                elif isinstance(dictobj, bool):
                    # TODO: Now assumed only Falses set, which do nothing
                    # In the future, we might include Trues on keywords.
                    # Such behavior goes here
                    pass
                elif isinstance(dictobj, (int, float, int64, float64)):
                    self.vars[group[-1]] = self.getue(group[-1])
                elif isinstance(dictobj, (bytes, str)):
                    try:
                        self.vars[group[-1]] = self.getue(group[-1])
                    except:
                        pass
                else:
                    self.unset_variables.append([group, dictobj])
            else:
                for key, value in dictobj.items():
                    dictobj = recursivereload(value, group + [key])
                return dictobj

        # Check whether data is read into memory
        if self.inplace:
            print('Cannot reload directly to HDF5 file with option "inplace".')
            print("Aborting")
            return
        if group is None:
            recursivereload()
        else:
            recursivereload(self.varinput[group], [group])
        if "commands" in locals():
            try:
                for command in commands:
                    exec(command)
            except:  # Exception as e:
                raise  # e
        #        try:
        #            self.varinput['setup']['commands'] = commands
        #            print('B4', bbb.albrb[16:19])
        #            for command in commands:
        #                exec(command)
        #            print('AFTER', bbb.albrb[16:19])
        #        except:
        #            pass
        # TODO: assess time lost here?
        for variable in self.vars.keys():
            try:
                self.packagelist[variable]
            except:
                self.packagelist[variable] = self.getpackage(variable, verbose=False)

    def load_inplace(self, fileobj=None, group=[]):
        """Creates dictionaries necessary for accessing HDF5 data"""
        from h5pickle import Group, File

        if fileobj is None:
            fileobj = self.hdf5case
        if isinstance(fileobj, File):
            for subgroup, data in fileobj.items():
                self.load_inplace(data, group + [subgroup])
        if isinstance(fileobj, Group):
            for subgroup, data in fileobj.items():
                self.load_inplace(data, group + [subgroup])
        else:
            self.vars[fileobj.name.split("/")[-1]] = fileobj.name

    def reinitializeplot(self, **kwargs):
        """Reinitializes the data of Subclasses for plotting"""
        # TODO: Return from UEDGE or variables?
        self.vertices = self.createpolycollection(self.getue("rm"), self.getue("zm"))

    def readyaml(self, fname, **kwargs):
        """Reads a YAML file and returns a nested dict

        Arguments
        ------------
        fname : str
            path and filename of YAML file to be read

        Returns
        ------------
        nested dict
        """
        from yaml import safe_load
        from pathlib import Path

        return safe_load(Path(fname).read_text())

    def openhdf5(self, fname, operation, **kwargs):
        """Opens HDF5 file and returns File object

        Arguments
        ------------
        fname : str
            path to/name of file to be opened
        operation : str
            operation to open the file for ('r'/'w'/'r+')

        Returns
        ------------
        h5py File object
        """
        from h5pickle import File
        from os.path import exists

        if exists(fname):
            return File(fname, operation)
        else:
            raise OSError('File "{}" not found!'.format(fname))

    def closehdf5(self, **kwargs):
        """Closes UeCase file hdf5case  that is being read"""
        try:
            self.hdf5case.close()
        except:
            raise OSError("No HDF5 file open")

    def read_hdf5_setup(self, fname):
        from h5pickle import File, Group
        from os.path import exists

        if exists(fname):
            savefile = File(fname, "r")
        else:
            raise OSError('File "{}" not found!'.format(fname))
        setup = savefile["setup"]
        ret = dict()

        def recursive_read_hdf5_setup(ret, setup, group=[]):
            if len(group) > 0:
                for subgroup in group:
                    try:
                        ret = ret[subgroup]
                    except:
                        ret[subgroup] = {}
                        ret = ret[subgroup]
            for name, content in setup.items():
                if isinstance(content, Group):
                    recursive_read_hdf5_setup(ret, content, group + [name])
                else:
                    ret[name] = content[()]

        recursive_read_hdf5_setup(ret, setup)
        return ret

    def setinput(
        self,
        setupfile=None,
        restore=True,
        savefname=None,
        readinput=True,
        restoresave=False,
        **kwargs
    ):
        """Sets all UEDGE variables from Case input"""
        """ Reads YAML input file

        Reads data from file to attribute setup.

        Arguments
        ------------
        setupfile : str
            path to/name of input file to be read

        Keyword arguments
        ------------
        restore : bool (default = True)
            switch whether to set UEDGE parameters to the read data 

        Returns
        ------------
        None
        """
        from collections import OrderedDict
        from copy import deepcopy
        from numpy import array

        if self.mutex() is False:
            return
        if readinput is True:
            if setupfile is None:
                setupfile = "{}.yaml".format(self.casename)
            try:
                self.varinput["setup"] = self.readyaml(setupfile)
                if "savefile" in self.varinput["setup"].keys():
                    savefname = self.varinput["setup"]["savefile"]
            except:
                self.varinput["setup"] = self.read_hdf5_setup(setupfile)
                self.restored_from_hdf5 = True
                savefname = setupfile
                self.setue("GridFileName", setupfile)
                self.setue("isgriduehdf5", 1)
        setup = deepcopy(self.varinput["setup"])

        # Pop out groups that cannot be parsed by default
        try:
            commands = setup.pop("commands")
        except:
            pass
        # TODO: tidy up casename definition
        try:
            self.casename = setup.pop("casename")
        except:
            pass
        try:
            self.savefname = setup.pop("savefile")
        except:
            pass
        if self.casename is None:
            self.casename = casename
        if isinstance(self.casename, bytes):
            self.casename = self.casename.decode("UTF-8")

        # TODO: Find a way to catch user-specified and radially varying
        #       diffusive coefficients when reading from file: userdifffname
        #       and radialdifffname attributes no available!
        def setinputrecursive(dictobj, group=[]):
            # NOTE: Something in this function is SLOOOW
            if not isinstance(dictobj, dict):
                # Skip UeCase-unique parameters
                if group[-1] not in ["userdifffname", "radialdifffname"]:
                    # Circumvent the padding with nulls for strings
                    try:
                        dictobj = dictobj.ljust(len(self.getue(group[-2])[group[-1]]))
                    except:
                        pass
                    # Set all other parameters
                    if isinstance(group[-1], int) and not isinstance(
                        dictobj, list
                    ):  # Set index-by-index
                        datalist = self.getue(group[-2])
                        datalist[group[-1]] = dictobj
                        self.setue(group[-2], datalist)
                    elif not isinstance(dictobj, list):
                        self.setue(group[-1], dictobj)
                    else:  # Set whole array
                        if isinstance(group[-1], int):
                            var = group[-2]
                            ind0 = group[-1]
                        else:
                            var = group[-1]
                            ind0 = 0
                        if len(self.getue(var).shape) > 1:
                            ueshape = self.getue(var).shape
                            if sum(ueshape) == sum(array(dictobj).shape):
                                try:
                                    self.setue(var, array(dictobj))
                                except:
                                    self.setue(var, array(dictobj).T)
                            elif len(ueshape) != len(array(dictobj).shape):
                                for j in range(ueshape[-1]):
                                    self.getue(var, cp=False)[:, j].put(
                                        range(ind0, len(dictobj) + ind0), dictobj
                                    )
                            else:
                                print(
                                    "!!! ERROR !!! Unable to determine "
                                    "shape of {} from input".format(var)
                                )
                                print(" *** ABORTING ***")
                        else:
                            # Edit array without copying == setue
                            self.getue(var, cp=False).put(
                                range(ind0, len(dictobj) + ind0), dictobj
                            )
                else:  # Set calls to restore diffusivities
                    setattr(self, group[-1], dictobj)
            else:
                for key, value in dictobj.items():
                    dictobj = setinputrecursive(value, group + [key])
                return dictobj

        # Set group order to assure proper allocation and avoid warnings
        for allokey in ["grid", "species"]:
            allolist = setup.pop(allokey)
            setinputrecursive(allolist)
            self.allocate()
        if restore is True:
            setinputrecursive(setup)
            self.allocate()
            if self.restored_from_hdf5 is True:
                conf = Config(verbose=False)
                print("=================================================")
                print("Restoring case from HDF5 file:")
                print("  Rate dirs read from .uedgerc")
                print("  Grid read from {}".format(setupfile))
                if self.getue("isbohmcalc") in [0, 1]:
                    print("  User-specified diffusivities read from HDF5 file")
                    self.userdifffname = setupfile
                elif self.getue("isbohmcalc") == 2:
                    print("  Radial diffusivities read from HDF5 file")
                    self.radialdifffname = setupfile
            # See if diffusivities unset despite being user-defined
            # If yes, try looking for them in the case being restored
            if (self.getue("isbohmcalc") in [0, 1]) and (self.userdifffname is None):
                self.userdifffname = self.casefname
            elif (self.getue("isbohmcalc") == 2) and (self.radialdifffname is None):
                self.radialdifffname = self.casefname
            if self.userdifffname:
                self.setuserdiff(self.userdifffname)
            if self.radialdifffname:
                self.setradialdiff(self.radialdifffname)
        # TODO: Can this mess be made somehow prettier?
        if restoresave is True:
            self.restoresave(savefname, **kwargs)
        self.reload()
        # NOTE:  Commands are executed as part of reload: don't repeat here

    #        try:
    #            for command in commands:
    #                exec(command)
    #        except:
    #            pass

    def setuserdiff(self, difffname, **kwargs):
        """Sets user-defined diffusivities

        Arguments
        ------------
        diffname : str
            HDF5 file from where to read 'diffusivities'/'bbb'/values
        """
        from h5py import File

        # TODO: replace with save-group function call?
        # NOTE: not sure why h5pickle throws error here?
        # No matter, we are only reading: use h5py

        if self.mutex() is False:
            return

        try:
            difffile = File(difffname, "r")
        except:
            difffile = File(self.casefname, "r")

        for variable in ["dif_use", "kye_use", "kyi_use", "tray_use"]:
            self.setue(variable, difffile["diffusivities"]["bbb"][variable][()])
        difffile.close()
        return

    def mutex(self, silent=False, **kwargs):
        """Returns bool whether case assigned to current UEDGE session

        Keyword parameters
        ------------------
        silent : boolean (default : False)
            Switch whether to issue mutex warning or not
        """
        if self.session_id == self.getue("session_id"):
            return True
        else:
            if silent is False:
                print(
                    "Mutex error! Object run-ID is {}, UEDGE run-ID "
                    "is {}. Aborting.".format(self.session_id, self.getue("session_id"))
                )
            return False

    def setradialdiff(self, difffname, **kwargs):
        """Sets radially varying diffusivities

        Arguments
        ------------
        diffname : str
            HDF5 file from where to read 'diffusivities'/'bbb'/values
        """

        if self.mutex() is False:
            return

        try:
            difffile = self.openhdf5(difffname, "r")
        except:
            difffile = self.openhdf5(self.casefname, "r")
        for variable in ["difniv", "kyev", "kyiv", "travisv"]:
            self.setue(variable, difffile["diffusivities"]["bbb"][variable][()])
            difffile.close()
        return

    def numvararr(self, variable):
        from numpy import array, transpose

        return transpose(
            array(variable)
            .reshape((self.get("ny") + 2, self.get("nx") + 2, self.get("numvar")))
            .T,
            (1, 2, 0),
        )

    def populate(self, silent=True, verbose=None, **kwargs):
        """Populates all UEDGE arrays by evaluating static 'time-step'"""
        from copy import deepcopy

        if self.mutex() is False:
            return

        if verbose is None:
            verbose = self.verbose

        if silent is True:
            self.setue("iprint", 0)
        issfon = deepcopy(self.getue("issfon"))
        ftol = deepcopy(self.getue("ftol"))
        self.setue("issfon", 0)
        self.setue("ftol", 1e20)
        self.exmain()
        self.setue("issfon", issfon)
        self.setue("ftol", ftol)
        self.update()
        if verbose is True:
            fnrm = sum(self.getue("yldot") ** 2) ** 0.5
            prtstr = "\n*** UEDGE arrays populated: {} ***"
            if fnrm < 10:
                print(prtstr.format("Case appears converged"))
                print("fnrm without preconditioning: {:.2e}\n".format(fnrm))
            elif fnrm < 100:
                print(prtstr.format("Warning, case may noy be fully " "converged"))
                print("fnrm without preconditioning: {:.1f}\n".format(fnrm))
            else:
                print(prtstr.format("WARNING, case NOT converged"))
                print("fnrm without preconditioning: {:.2e}\n".format(fnrm))
        if silent is True:
            self.setue("iprint", 1)

    def restore(self, inputfname=None, savefname=None, populate=True, **kwargs):
        """Restores a full case into memory and object"""
        if self.mutex() is False:
            return
        self.setinput(inputfname, savefname=savefname, restoresave=True, **kwargs)
        if populate is True:
            self.populate(silent=True, **kwargs)

    '''
    def savevar(self, savefile, groups, variable, data, **kwargs):
        """ Saves variable and metadata to HDF5 group and dataset

        Arguments
        ------------
        savefile : h5py File object
            h5py object where to save data
        groups : list of str
            list of subroups to create and store variable
        variable : str
            name of dataset to create. If same as UEDGE variable name, 
            metadata (description, units) are stored as dataset 
            attributes
        data : array/int/str/float
            data to be stored to dataset
        
        Returns
        ------------
        None
        """
        # Make group, package, file into a list and iterate?
        output = []
        for group in groups:   
            savefile.require_group(group)
            savefile = savefile[group]
            output.append(group)
        if variable not in savefile.keys():
            # Check whether variable is allocated or not
            try:
                if self.getpackobj(variable, verbose=False).allocated(\
                    variable) != 0:
                    savefile.create_dataset(variable, data=data)
                else:
                    savefile.create_dataset(variable, data=0)
                    # TODO: Omit or save False?
                    # TODO: Verbose prompt or not?
            except:
                savefile.create_dataset(variable, data=data)
        try:
            savefile.attrs[variable] = \
                packageobject(package).getvarunit(variable)
        except:
            pass
        try:
            savefile.attrs[variable] = \
                packageobject(package).getvardoc(variable)
        except:
            pass

    def savegroup(self, savename, group, append=True, **kwargs):
        # TODO: Do we want to always rewrite or try to append?
        savefile = self.openhdf5(savename, 'a'*(append is True) + \
            'w'*(append is False))
        if append is False:
            self.savemetadata(savefile)
        self.recursivesave(savefile, self.varinput[group], [group])
        savefile.close()

    def savesetup(self, savename, **kwargs):
        self.savegroup(savename, 'setup', **kwargs)
    
    def savegrid(self, savename, **kwargs): 
        self.savegroup(savename, 'grid', **kwargs)

    def savediffusivities(self, savename, **kwargs):
        self.savegroup(savename, 'diffusivities')

    def recursivesave(self, savefile, saveobj, group = [], **kwargs):
        # Bottom level of structure
        if not isinstance(saveobj, dict):
            # Special setup for saving setup parameters to store
            # actual set values of any setup parameters defined in input
            if group[0]=='setup':
                variable = group.pop(-1)
                if variable in ['userdifffname', 'radialdifffname', 
                    'casename', 'commands', 'savefile']:
                    value = saveobj
                else:
                    value = self.getue(variable)
                self.savevar(savefile, group, variable, value)
            # Store requested data 
            elif isinstance(saveobj, list):
                for variable in saveobj:
                    self.savevar(savefile, group, variable, 
                        self.getue(variable))
        # Recursively go deeper in structure
        else:
            for key, value in saveobj.items():
                if isinstance(key, int):
                    # Save the full array once only
                    variable = group.pop(-1)
                    self.savevar(savefile, group, variable,
                        self.getue(variable))
                    return
                saveobj = self.recursivesave(savefile, value, group + [key])
            return saveobj 

    def savemetadata(self, savefile, **kwargs):
        from uedge import __version__
        from time import time, ctime
        try:
            savefile.attrs['casename'] = self.casename
        except:
            pass
        try:
            savefile.attrs['savefile'] = self.savefile
        except:
            pass
        savefile.attrs['UETOOLS_version'] = self.uetoolsversion
        savefile.attrs['time'] = time()
        savefile.attrs['ctime'] = ctime()
        savefile.attrs['code'] = 'UEDGE'
        savefile.attrs['ver'] = self.uedge_ver
        savefile.attrs['pyver'] = self.pyver
        savefile.attrs['user'] = self.user
        savefile.attrs['hostname'] = self.hostname
        try:
            savefile.attrs['location'] = self.location
        except:
            pass

    def save(self, savefname, group=None, append=False, **kwargs):
        """ Saves HDF5 file containing UeCase data
        
        Arguments
        ------------
        savefname : str 
            path to/name of file to write data to 

        Keyword arguments
        ------------
        group : str (default = None)
            group identifier of group to be written to file. If None,
            all data stored in UeCase is written
        """
        if self.inplace:
            print('Data read from file, no data to save. Aborting.')
            return
        savefile = self.openhdf5(savefname, 'w'*(append is False) \
            + 'a'*(append is True)) 
        self.savemetadata(savefile)
        if group is None:
            self.recursivesave(savefile, self.varinput)
        else:
            self.recursivesave(savefile, self.varinput[group], [group])
        savefile.close()
    '''
