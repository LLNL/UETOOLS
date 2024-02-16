from .CasePlot import Caseplot
from .Solver import Solver
from .Save import Save
from uetools.UeDashboard import Dashboard
from uetools.UeUtils.Misc import Misc
from uetools.UeUtils.Tracker import Tracker
from uetools.UeUtils.Convert import Convert
from uetools.UeUtils.RadTransp import RadTransp
from uetools.UeUtils.Interpolate import Interpolate
from uetools.UeUtils.ConvergeStep import ConvergeStep
from uetools.UeUtils.AboutSetup import AboutSetup
from uetools.UePostproc.Postproc import PostProcessors
from uetools.UePostproc.ADAS import ADAS
from uetools.UeConfig.Config import Config
try:
    from uedge import bbb, com, aph, api, svr
    uedge_is_installed = True
except:
    uedge_is_installed = False
    
try:
    from uedge import __version__
except:
    __version__ = 'N/A'
try:
    from Forthon import packageobject
except:
    pass

# TODO: make yaml read/write/get/set case-insensitive?
# TODO: Consider compression of data

class Case(Misc, Save, PostProcessors, ConvergeStep, ADAS, 
    RadTransp, Interpolate, Convert, Tracker, Config, Caseplot, Solver,
    AboutSetup
):
    """ UEDGE Case container object.

    Subclasses
    ----------
    Solver -- contains time-stepping and convergence routines
    PostProcessors -- contains useful post-processing routines
    ConvergeStep -- iterative advancement using time-dependent solver
    ADAS -- routines for use with ADAS atomic data files
    RadTransp -- routines for radial transport coefficient fitting
    Interpolate -- grid interpolation routines
    Convert -- writes UETOOLS data to Python/YAML input file
    Tracker -- routines for tracking changes to UEDGE input variables
    Config -- UETOOLS configuration routines (for .uetoolsrc-file)
    Caseplot -- plotting routines using Case object functionalities
    Misc -- miscellaneous UETOOLS utilities
    Save -- routines for saving and loading UEDGE states
    AboutSetup -- utility to display information about case setup

    Attributes
    ----------
    allocate : function
        allocates the UEDGE arrays based on the input
    filename : string
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
    -------
    assign(**kwargs)
        assignes the UEDGE memory to the Case object
    get_inplace(variable, s=None, **kwargs)
        returns the value of variable from HDF5 file
    get_memory(variable, s=None, **kwargs)
        returns the value of variable from Case memory
    getsetue_inplace(*args, **kwargs)
        placeholder for returning UEDGE variables when reading inplace
    get_uememory(variable, s=None, cp=True, **kwargs)
        returns the value of variable from UEDGE memory
    load_inplace(fileobj=None, group=[])
        creates the necessary dictionaries for accessing HDF5 data
    mutex(silent=False, **kwargs)
        returns True if UEDGE memory assigned to this Case object
    openhdf5(fname, operation, **kwargs)
        opens HDF5 file and returns file object
    read_hdf5_setup(fname)
        reads the UEDGE input deck from setup group of HDF5
    readyaml(fname, **kwargs)
        returns nested dict read from YAML fname
    reload
    restore_input(fname=None, savefile=None, populate=True, **kwargs)
        sets UEDGE input parameters according to setup attribute
    restore_save(savefile, **kwargs)
        restores the UEDGE state from HDF5 savefile
    setinput(   setupfile=None, restore=True, savefile=None, 
                readinput=True, restoresave=False, **kwargs
        )
        Reads a YAML input file and sets up UEDGE case
    set_radialdiff(fname, **kwargs)
        sets radial diffusion coefficient profiles as defined HDF5 file
    set_uememory(variable, value, **kwargs)
        sets the UEDGE variable in memory to value
    set_userdiff(fname, **kwargs)
        sets user-defined diffusion coefficients as defined HDF5 file
    update    
        checks UEDGE state and updates Case variables if necessary

    Variables
    ---------
    allocate: wrapper of UEDGE function allocate
    aphdir: path to hydrogenic rates used by Case objects
    apidir: path to impurity rates used by Case objects
    casename: string identifier for case
    diff_file: path to file containing diffusivity data
    exmain_evals: the number of exmain evaluations perfored by UEDGE
    filename: path to YAML input/HDF5 file read
    get: wrapper for function to get data depending on setup
    getue: wrapper for function to get UEDGE data depending on setup
    hdf5case: HDF5 File object when running in inplace-mode
    hostname: machine name for writing figure labels
    inplace: boolean for controling case setup and behavior
    location: cwd where data read to Case object is located
    omitvars: special named variables to not be parsed from input YAML
    packagelist: list of UEDGE packages available
    pyver: version of python UEDGE package/build
    restored_from_hdf5: boolean set to true if case read from HDF5 file
    savefile: path to HDF5 file containing save variables
    session_id: mutex ID of Case object
    uedge_ver: internal UEDGE version uedge_ver
    uetoolsversion: version of UETOOLS (internal)
    uevars: nested dict of UEDGE packages, variables, and hashes
    unset_variables: UEDGE variables that have not been allocated
    use_mutex: boolean instructing Case to perform mutex if True
    user: user name for writing to figure labels
    varinput: nested dict of variables to be saved to save files
    vars: data for all variables tracked by Case object
    verbose: boolean telling Case to operate silently if Dalse
     

    Side-effects
    ------------

    - Modifies UEDGE internal state, shared among all Case objects.

    - Some methods need to change working directory when
      calling UEDGE. In all cases the change should be reverted
      before returning, so that there is no net change.

    """

    def __init__(
        self,
        filename=None,
        inplace=False,
        variableyamlfile=None,
        casename=None,
        assign=True,
        verbose=True,
        savefile=None,
        diff_file=None,
        aphdir=None,
        apidir=None,
        **kwargs,
    ):
        """Initializes the UeCase object.

        Keyword arguments
        -----------------
        filename : str (default =  None)
            Path to Case object initializer. Reads YAML input file if
            path points to YAML file. Reads data from HDF5 file if path
            points to HDF5 file. Reads data from memorty if None.
        inplace : bool (default = False)
            Switch whether to read data from file into UeCase memory
            (False) or get data using file I/O at every call (True)
        variableyamlfile : str (default = None)
            Path to YAML file containing definitions of data and
            variables to be read. If None, accesses the module defaults
        assign : bool (default = True)
            Switch whether to assign the current run to the caseobject
        verbose : bool (default = True)
            Silences Case object if verbose = False
        savefile : str (default = None)
            Path to HDF5 file containing saved UEDGE state. Redundant
            if reading case form HDF5. Read from YAML input if 
            available (kwarg takes precedence).
        diff_file : str (default = None)
            Path to HDF5 file containing diffusivity coefficients. 
            Ignored if isbohmcalc != 0 or 2
        aphdir : str (default = None)
            Path to directory with hydrogenic rate files. Kwarg 
            defintion takes precedence over input file and run command
            files. Reads aphdir from run comand files or input if None.
        apidir : str (default = None)
            Path to directory with impurity rate files. Kwarg 
            defintion takes precedence over input file and run command
            files. Reads apidir from run comand files or input if None.

        """
        import uetools
        import os
        from os.path import exists, abspath
        from os import getlogin, getcwd
        from socket import gethostname
        from matplotlib.pyplot import ioff, ion

        # Assert input is correct before proceeding
        if filename is not None:
            if not exists(filename):
                raise ValueError('File {} does not exist!'.format(\
                    filename
                ))
    
        if (inplace is False) and (uedge_is_installed is False):
            print("No working UEDGE install found: only "+\
                "inplace-evaluation is possible.")
            print("Only UETOOLS HDF5 saves can be restored.")
            print("For more information, consult UETOOLS documentation.")
            inplace = True

        self.configcase(verbose=verbose)
        # TODO: add label attribute
        # Read from uedge.label and strip
        # Use to restore and/or save files
        # Initialize class attributes
        # Checksum whether to update data

        # Initialize parameters
        self.casename = casename
        self.savefile = savefile
        self.inplace = inplace
        self.verbose = verbose
        self.restored_from_hdf5 = False
        self.uetoolsversion = "1.1.2"  # UEtools version
        try:
            self.allocate = packageobject("bbb").getpyobject("allocate")
        except:
            pass
        self.filename = filename
        self.diff_file = diff_file
        try: 
            self.apidir
            if apidir is not None:
                self.apidir = aphdir
        except:
            self.apidir = apidir
        try: 
            self.aphdir
            if aphdir is not None:
                self.aphdir = aphdir
        except:
            self.aphdir = aphdir
        self.hdf5case = None
        self.pyver = __version__
        try:
            self.uedge_ver = (
                packageobject("bbb").getpyobject("uedge_ver")[0].strip().decode("UTF-8")
            )
        except:
            self.uedge_ver = 'unknown'

        try:
            self.user = getlogin()
        except OSError:
            # Can fail on e.g cluster nodes
            self.user = "unknown"
        try:
            self.hostname = gethostname()
        except OSError:
            self.hostname = "unknown"

        self.unset_variables = []
        self.omitvars = [
            "userdifffname",
            "radialdifffname",
            "diff_file",
            "casename",
            "commands",
            "savefile",
        ]

        self.vars = dict()
        self.varinput = dict()
        self.packagelist = dict()
        # TODO: add hostname, mcfilename, aphfname

        # Set up structure for reading/writing data
        # Load all data to object in memory
        if self.inplace is False:
            self.get = self.get_memory
            self.getue = self.get_uememory
            self.setue = self.set_uememory
            try:
                # Get the directory containing the input file
                self.location = os.path.dirname(abspath(self.filename))
            except:
                self.location = getcwd()
            try:
                self.exmain_evals = self.getue("exmain_evals")
                self.use_mutex = True
            except:
                print('Variable "exmain_evals" not found!')
                print('Using UEDGE version <7, deactivate mutex')
                self.use_mutex = False
            if self.use_mutex is True:
                self.session_id = self.getue("max_session_id") + 1
                setattr(
                    packageobject("bbb"), "max_session_id", 
                    self.getue("max_session_id") + 1
                )

            if assign is True:
                self.assign()


            
            self.varinput = self.readyaml("{}/{}".format(
                uetools.__path__[0], "yamls/requiredvariables.yaml"
            ))
            # Read YAML to get variables to be read/stored/used
            if variableyamlfile is None:  # No YAML variable file requested
                if hasattr(self, "variableyamlfile"):
                    self.varinput.update(self.readyaml(self.variableyamlfile))
                else:
                # Use default: find package location and package YAMLs
                    self.varinput.update(self.readyaml("{}/{}".format(
                        uetools.__path__[0], "yamls/standardvariables.yaml"
                    )))
            else:  # YAML specified, use user input
                self.varinput.update(self.readyaml(variableyamlfile))  # Read to memory

            if self.filename is not None:
                self.restore_input(self.filename, self.savefile)
            else:
                self.reload()
                # TODO:
                if uedge_is_installed and not self.inplace:
                    self.get_uevars()
        # Read all data directly from HDF5 file
        else:
            self.get = self.get_inplace
            self.set = self.getsetue_inplace
            self.getue = self.getsetue_inplace
            self.setue = self.getsetue_inplace
            self.filename = abspath(self.filename)
            self.location = os.path.dirname(self.filename)
            if exists(self.filename):
                try:
                    self.hdf5case = self.openhdf5(self.filename, "r")
                except Exception as err:
                    print("Unable to open {}. Aborting.".format(self.filename))
                    raise err
            self.load_inplace()
            if self.filename is None:
                print("Must specify data file when inplace=True! Aborting.")
                return
        self.snull = (self.get('geometry')[0].decode('UTF-8').strip() \
            in ['uppersn', 'snull'])
#        if self.snull:
#            self.ixpt1 = self.get('ixpt1')[0]
#            self.ixpt2 = self.get('ixpt2')[0]
#            self.iysptrx = self.get('iysptrx')
#        self.nx = self.get('nx')
#        self.ny = self.get('ny')
#        self.ixmp = self.get('ixmp')
        # Initialize parent classes
        super().__init__(**kwargs)


    # NOTE: Update class data, or try reading from forthon first??
    def update(self, **kwargs):
        """Checks if UEDGE state has changed and updates as needed.

        Modifies
        --------
        Calls Case.reload, updates Case.vertices

        Returns
        -------
        """

        if self.exmain_evals != self.getue("exmain_evals"):
            self.exmain_evals = self.getue("exmain_evals")
            if self.mutex() is False:
                raise Exception("Case doesn't own UEDGE memory")
            self.reload()
            self.vertices = self.createpolycollection(
                                self.get("rm"), 
                                self.get("zm")
                            )

    def get_inplace(self, variable, s=None, **kwargs):
        """Returns variable from HDF5 file

        Arguments
        -----------------
        variable : str
            Variable to return from HDF5 file

        Keyword arguments
        -----------------
        s : int (None)
            Species index, if variable is multi-species array

        Returns
        -------
        ndarray containing data or None if not found
        """
        from numpy import ndarray
        from h5py import File
        
        # TODO: figure out why some cases spontaneously close?
        if 'Closed' in str(self.hdf5case):
            self.hdf5case = File(self.filename, 'r')

        try:
            retvar = self.hdf5case[self.vars[variable]][()]
        except:
            print("{} not found in {}".format(variable, self.filename))
            return

        if isinstance(retvar, (ndarray, list)):
            if len(retvar.shape) == 3:
                if s is not None:
                    retvar = retvar[:, :, s]
        return retvar

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
        ---------
        variable : str
            name of variable to be returned

        Keyword arguments
        -----------------
        s : int (default = None)
            species index to be returned for 3-dimensional arrays. If
            None, get returns the full 3-dimensional array. Otherwise,
            get returns the s:th index of the third dimension

        Returns
        -------
        UeCase value of variable
        """
        # TODO: search input too? Store input to Vars?
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
        """Assigns the UEDGE session to this object

        Modifies
        --------
        UEDGE memory : bbb.session_id set to Case.session_id

        Returns
        -------
        None
        """
        if self.use_mutex is True:
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
        """Placeholder to avoid getting/setting when reading inplace."""
        raise Exception("Cannot set/get UEDGE values when reading from "+\
                            "HDF5 file")

    def set_uememory(self, variable, value, **kwargs):
        """Sets the Forthon variable in package to data

        Arguments
        ---------
        variable : str
            variable name
        value : array/list/float/int
            data to be written to the UEDGE variable

        Modifies
        --------
        UEDGE memory : sets variable value to value

        Returns
        -------
        None
        """
        try:
            package = self.packagelist[variable]
        except:
            package = self.getpackage(variable, verbose=False)
        if self.mutex():
            try:
                setattr(packageobject(package), variable, value)
            except Exception as e:
                raise KeyError("{} could not be set: {}".format(variable, e))

    def get_uememory(self, variable, s=None, cp=True, **kwargs):
        """Retrieves data from UEDGE variable in package.

        Arguments
        ---------
        variable : str
            variable name

        Keyword arguments
        -----------------
        s : int (default = None)
            species index to return if requested array is 3D species dependent
        cp : boolean (default = True)
            returns a copy of the values if True, a pointer to variable
            if False

        Returns
        -------
        value of/pointer to UEDGE variable
        """
        # TODO: fix get on functions!
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

        Keyword arguments
        -----------------
        group : str (default = None)
            group specifier to reload. If None, reloads all
            variables in Case.vars

        Modifies
        --------
        self.vars : dict
            Dictionary of values, with keys specified in self.varinput

        Returns
        -------
        None

        """
        from numpy import ndarray, int64, float64

        # Check whether data is read into memory
        if self.inplace:
            raise Exception(
                'Cannot reload directly to HDF5 file with option "inplace".'
            )
        # Check that the case is assigned UEDGE memory
        if self.mutex is False:
            raise Exception("Case doesn't own UEDGE memory")

        def recursivereload(dictobj, group=[]):
            """ Recursively traverses dictionary and stores UEDGE data to self """
            if not isinstance(dictobj, dict):
                # Reached bottom of nested dictionaries: determine format
                if isinstance(dictobj, (list, ndarray)):
                    # We have a list: either list of variables to store or
                    # list defning the variable array
                    if self.getpackage(group[-1], verbose=False) != None:
                        # Request to set array starting from index 0:
                        # just read the variable into memory
                        self.vars[group[-1]] = self.getue(group[-1])
                    elif isinstance(group[-1], int):
                        # Setting subarray, store variable
                        self.vars[group[-2]] = self.getue(group[-2])
                    else:
                        # List of variables, store each
                        for variable in dictobj:
                            self.vars[variable] = self.getue(variable)
                elif isinstance(group[-1], int):
                    if len(group) > 2 and isinstance(group[-2], int):
                        self.vars[group[-3]] = self.getue(group[-3])
                    else:
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
                    recursivereload(value, group + [key])
        # Pop out any custom commands, as these cannot be reloaded (not vars)
        try:
            commands = self.varinput["setup"].pop("commands")
        except:
            pass
        # Reload the variables recurively
        if group is None:
            recursivereload(self.varinput)
        else:
            recursivereload(self.varinput[group], [group])
        # If there were any custom commands, put them back where they belong
        try:
            self.varinput["setup"]["commands"] = commands
        except:
            pass
        # Update the dict containing the package containing each variable
        for variable in self.vars.keys():
            if variable not in self.packagelist:
                self.packagelist[variable] = self.getpackage(
                                                    variable, 
                                                    verbose=False
                                            )

    def load_inplace(self, fileobj=None, group=[]):
        """Creates dictionaries necessary for accessing HDF5 data

        Recursively reads the supplied HDF5 file and maps each variable
        to its location in the file.

        Keyword arguments
        -----------------
        fileobj : HDF5 File object (default = None)
            
        Modifies
        --------
        Case.vars dictionary : adds variables as keys with paths in
            file as items.

        Returns
        -------
        None
        """
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

    def openhdf5(self, fname, operation, **kwargs):
        """Opens HDF5 file and returns File object

        Arguments
        ---------
        fname : str
            path to/name of file to be opened
        operation : str
            operation to open the file for ('r'/'w'/'r+')

        Returns
        -------
        h5py File object
        """
        from h5pickle import File
        from os.path import exists

        if exists(fname):
            return File(fname, operation)
        else:
            raise OSError('File "{}" not found!'.format(fname))

    def read_hdf5_setup(self, fname):
        """ Reads the UEDGE input deck from setup group of HDF5 

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
        from h5pickle import File, Group
        from os.path import exists

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
        ret = dict()
        if not exists(fname):
            raise OSError('File "{}" not found!'.format(fname))
        else:
            with File(fname, "r") as savefile:
                recursive_read_hdf5_setup(ret, savefile["setup"])
            del savefile
        self.savefile = fname
        return ret

    def setinput(
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
        # Extract user-supplied casename and diff_file
        casename = deepcopy(self.casename)
        diff_file = deepcopy(self.diff_file)
        if self.mutex() is False:
            raise Exception("Case doesn't own UEDGE memory")

        # TODO: check whether file is HDF5 instead of Try
        if readinput is True:
            if setupfile is None:
                print("No setup file specified:")
                if self.casename is None:
                    raise ValueError("    No casename defined: aborting!")
                else:
                    print("    Using casename '{}'".format(self.casename))
                setupfile = "{}.yaml".format(self.casename)
            if is_hdf5(setupfile):
                self.varinput["setup"] = self.read_hdf5_setup(setupfile)
                self.restored_from_hdf5 = True
                self.setue("GridFileName", setupfile)
                self.setue("isgriduehdf5", 1)
            else:
                try:
                    self.varinput["setup"] = self.readyaml(setupfile)
                except Exception as e:
                    raise ValueError(f"Input file could not be parsed: {e}")
                
#            try:
#                self.varinput["setup"] = self.readyaml(setupfile)
#            except:
#                self.varinput["setup"] = self.read_hdf5_setup(setupfile)
#                self.restored_from_hdf5 = True
#                self.setue("GridFileName", setupfile)
#                self.setue("isgriduehdf5", 1)
        setup = deepcopy(self.varinput["setup"])

        # Pop out groups that cannot be parsed by default
        try:
            commands = setup.pop("commands")
        except:
            pass
        try:
            detected = setup.pop('detected')
        except:
            pass
        # TODO: Add mist.dat as an optional parameter/etc to allow changing
        #       the name/path to the data file
        def setinputrecursive(dictobj, group=[]):
            if not isinstance(dictobj, dict):
                # Skip UeCase-unique parameters
                if group[-1] not in self.omitvars + ['chgstate_format']:
                    # NOTE: Not sure what to do with chgstate_format, fauls for some strange reason...
                    # NOTE: Should not be an input, just skip for the time being
                    # Avoid overwriting grid path when restoring from HDF5
                    if (group[-1]=='GridFileName') and\
                        (self.restored_from_hdf5 is True):
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
                                raise ValueError(
                                    "!!! ERROR !!! Unable to determine "
                                    "shape of {} from input".format(var)
                                )
                        else:
                            # Edit array without copying == setue
                            self.getue(var, cp=False).put(
                                range(ind0, len(dictobj) + ind0), dictobj
                            )
                    elif isinstance(group[-1], int):
                        # Set a single entry in a 1D or 2D array
                        if len(group) > 2 and isinstance(group[-2], int):
                            # 2D array
                            datalist = self.getue(group[-3])
                            datalist[group[-2], group[-1]] = dictobj
                            self.setue(group[-3], datalist)
                        else:
                            # 1D array
                            datalist = self.getue(group[-2])
                            datalist[group[-1]] = dictobj
                            self.setue(group[-2], datalist)
                    elif dictobj is None:
                        print("WARNING Unset specifier in input:", group[-1])
                    else:
                        self.setue(group[-1], dictobj)

                else:  # Set calls to restore diffusivities
                    if (group[-1]=='savefile') and \
                        (self.restored_from_hdf5 == True):
                        pass
                    else:
                        setattr(self, group[-1], dictobj)
            else:
                for key, value in dictobj.items():
                    dictobj = setinputrecursive(value, group + [key])
                return dictobj

        # Set group order to assure proper allocation and avoid warnings
        for allokey in ["grid", "species"]:
            try:
                allolist = setup.pop(allokey)
            except KeyError:
                print(f"WARNING: Expecting setup group '{allokey}'")
                continue
            setinputrecursive(allolist)
            self.allocate()

        if restore is True:
            setinputrecursive(setup)
            self.allocate()
            if "detected" in locals():
                self.detected = detected
                setinputrecursive(detected)
            if isinstance(self.casename, bytes):
                self.casename = self.casename.decode("UTF-8")
            if isinstance(self.savefile, bytes):
                self.savefile = self.savefile.decode("UTF-8")
            if self.restored_from_hdf5 is True:
                print("=================================================")
                print("Restoring case from HDF5 file:")
                print("  Rate dirs read from .uedgerc")
                print("  Grid read from {}".format(setupfile))
                self.diff_file = setupfile
            # Override with diff_file maually defined diff_file upon
            if diff_file is not None:
                self.diff_file = diff_file
            # Otherwise, try setting accoridng to input
            else:
                # diff_file takes precedence
                if self.diff_file is None:
                    # if diff_file not set, override using old diffusion files
                    # userdifffname takes precedence over radialdifffile
                    # if both are present
                    if hasattr(self, "radialdifffname"):
                        if (self.radialdifffname is not None) \
                            and (self.radialdifffname is not False
                        ):
                            self.diff_file = self.radialdifffname
                        del self.radialdifffname
                    if hasattr(self, "userdifffname"):
                        if (self.userdifffname is not None) \
                            and (self.userdifffname is not False
                        ):
                            self.diff_file = self.userdifffname
                        del self.userdifffname
            if (self.diff_file is None) and (self.getue("isbohmcalc") in [0,2]):
                self.diff_file = self.savefile
                print('No diffusivity-file supplied: reading from '+\
                    'save-file "{}"'.format(self.diff_file))
            # Set diffusivities based on file if model requires profiles
            if self.getue("isbohmcalc") == 0:
                print("  User-specified diffusivities read from HDF5 "+\
                    'file "{}"'.format(self.diff_file))
                try: 
                    self.setuserdiff(self.diff_file)
                except Exception as e:
                    print(
                        f"WARNING: failed to read diffusivities "+\
                            f"from {self.diff_file}: {e}"
                    )
            elif self.getue("isbohmcalc") == 2:
                print("  Radial diffusivities read from HDF5 file "+\
                    '"{}"'.format(self.diff_file))
                try:
                    self.set_radialdiff(self.diff_file)
                except Exception as e:
                    print(
                        f"WARNING: failed to read radial diffusivities "+\
                            f"from {self.diff_file}: {e}"
                    )
        if casename is not None:
            self.casename = casename
        if savefile is not None:
            self.savefile = savefile
        if self.aphdir is not None:
            self.setue("aphdir", self.aphdir)
        if self.apidir is not None:
            self.setue("apidir", self.apidir)
        if restoresave is True:
            if (self.savefile is None) and (self.get('restart') == 1):
                raise ValueError("No save-file supplied!")
            elif self.get('restart') == 1:
                self.load_state(self.savefile, **kwargs)
        if uedge_is_installed and not self.inplace:
            self.get_uevars()
        # NOTE: Get the hashes before running any commands. This way,
        # any changes done in external scripts etc will be registered.
        # This is useful (and necessary) to capture changes to arrays
        # being modified.
        self.reload()
        if not self.restored_from_hdf5:
            if "commands" in locals():
                for command in commands:
                    try:
                        exec(command)
                    except Exception as e:
                        print(f"Command {command} failed: {e}")
        

    def setuserdiff(self, difffname, **kwargs):
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
        from h5py import File
        from os.path import exists

        # TODO: replace with save-group function call?
        # NOTE: not sure why h5pickle throws error here?
        # No matter, we are only reading: use h5py

        if self.mutex() is False:
            raise Exception("Case doesn't own UEDGE memory")
        

        if not exists(difffname):
            difffname = self.filename
            if not exists(self.filename):
                raise Exception("Diffusivity file not found!")
    
        with File(difffname) as file:
            for variable in ["dif_use", "kye_use", "kyi_use", "tray_use"]:
                self.setue(variable, file[f"diffusivities/bbb/{variable}"][()])
            for variable in ["difni", "kye", "kye", "travis", "fcdif"]:
                if variable in file["diffusivities/bbb"]:
                    self.setue(variable, 
                        file[f"diffusivities/bbb/{variable}"][()]
                    )

    def mutex(self, silent=False, **kwargs):
        """Returns True if case assigned to current UEDGE session.

        Keyword parameters
        ------------------
        silent : boolean (default : False)
            Switch whether to issue mutex warning or not

        Returns
        -------
        True if Case object own UEDGE memore, False otherwise
        """
        if self.use_mutex == False:
            return True
        if self.session_id == self.getue("session_id"):
            return True
        else:
            if silent is False:
                print(
                    "Mutex error! Object run-ID is {}, UEDGE run-ID "
                    "is {}. Aborting.".format(self.session_id, self.getue("session_id"))
                )
            return False

    def set_radialdiff(self, difffname, **kwargs):
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
            difffname = self.filename
            if not exists(self.filename):
                raise Exception("Diffusivity file not found!")
        
        with File(difffname) as file:
            for variable in ["difniv", "kyev", "kyiv", "travisv"]:
                self.setue(variable, file["diffusivities"]["bbb"][variable][()])


    def populate(self, silent=True, verbose=None, **kwargs):
        """ Populates all UEDGE arrays by evaluating static 'time-step'

        Outputs prompt assessing whether case is converged or not.

        Keyword arguments
        -----------------
        silent : bool (default = True)
            Tells Case object to silence UEDGE exmain writes.
        verbose : bool (default = None)
            Tells Case to output UETOOLS prompts if True. If 
            verbose = None, uses Case.verbose default.

        Modifies
        --------
        UEDGE memory : updates all UEDGE array to correspond to 
            restored state

        Returns
        -------
        None
        """
        from copy import deepcopy
        import os

        if self.mutex() is False:
            raise Exception("Case doesn't own UEDGE memory")

        if verbose is None:
            verbose = self.verbose

        if silent is True:
            old_iprint = self.getue("iprint")
            self.setue("iprint", 0)

        issfon = deepcopy(self.getue("issfon"))
        ftol = deepcopy(self.getue("ftol"))
        try:
            self.setue("issfon", 0)
            self.setue("ftol", 1e20)
            self.exmain()
            self.setue("gengrid", 0)  # Ensure that grids aren't generated again
        finally:
            # Ensure that original settings and working directory are restored
            self.setue("issfon", issfon)
            self.setue("ftol", ftol)

        self.update()  # Reloads variables from UEDGE

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
            self.setue("iprint", old_iprint)  # Restore

    def restore_input(self, inputfname=None, savefile=None, 
        populate=True, **kwargs):
        """ Restores a full case into memory and object.

        Keyword arguments
        -----------------
        inputfname : str (default = None)
            Path to the input file to be read.
        savefile : str (default = None)
            Path to HDF5 file containing save data. If None,
            savefile is read from input file.
        populate : bool (default = True)
            Tells UETOOLS to populate all UEDGE arrays after reading
            input and save file.

        Modifies
        --------
        Case object : local arrays are updated to correspond to input
            and local UEDGE variables
        UEDGE memory : UEDGE memory is updated to correspond to input
            and saved state

        Returns
        -------
        None
        """
        if self.mutex() is False:
            raise Exception("Case doesn't own UEDGE memory")

        self.setinput(inputfname, savefile=savefile, restoresave=True, 
                **kwargs
        )
        if populate is True:
            self.populate(silent=True, **kwargs)

    def restore_save(self, savefile, **kwargs):
        """ Procedure to read saved state and restore UEDGE variables.

        Arguments
        -----------------
        savefile : str 
            Path to HDF5 file containing UEDGE state data

        Modifies
        --------
        UEDGE state

        Returns
        -------
        None
        """
        self.load_state(savefile, **kwargs)
        self.populate(**kwargs)
