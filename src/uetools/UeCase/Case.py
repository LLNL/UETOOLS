from uetools.UePlot import Caseplot
from uetools.UeSolver import Solver
from .Save import Save
from .Config import Config
from .Input import Input
from uetools.UeUtils import *
from uetools.UePostproc.Postproc import PostProcessors
from uetools.UeDiagnostics.ADAS import ADAS
import uetools
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
# TODO: Store all *var* members to var-dict 
# TODO: Re-introduce the verices

class Case:
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
    read_hdf5_setup(fname)
        reads the UEDGE input deck from setup group of HDF5
    reload
    restore_input(fname=None, savefile=None, populate=True, **kwargs)
        sets UEDGE input parameters according to setup attribute
    restore_save(savefile, **kwargs)
        restores the UEDGE state from HDF5 savefile
    setinput(   setupfile=None, restore=True, savefile=None, 
                readinput=True, restoresave=False, **kwargs
        )
        Reads a YAML input file and sets up UEDGE case
    setradialdiff(fname, **kwargs)
        sets radial diffusion coefficient profiles as defined HDF5 file
    set_uememory(variable, value, **kwargs)
        sets the UEDGE variable in memory to value
    setuserdiff(fname, **kwargs)
        sets user-defined diffusion coefficients as defined HDF5 file
    update    
        checks UEDGE state and updates Case variables if necessary

    Variables
    ---------
    aphdir: path to hydrogenic rates used by Case objects
    apidir: path to impurity rates used by Case objects
    casename: string identifier for case
    diff_file: path to file containing diffusivity data
    exmain_evals: the number of exmain evaluations perfored by UEDGE
    filename: path to YAML input/HDF5 file read
    get: wrapper for function to get data depending on setup
    getue: wrapper for function to get UEDGE data depending on setup
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
        diffusivity_file=None,
        aphdir=None,
        apidir=None,
        restoresave=True,
        store_defaults=False,
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
        restoresave : bool (default = True)
            Switch whether to restore save file or not during reading
        """
        import uetools
        import os
        from copy import deepcopy
        from os.path import exists, abspath
        from os import getlogin, getcwd
        from socket import gethostname
        from matplotlib.pyplot import ioff, ion

        # Set up nested dict for variable operations
        self.variables = {
            'stored': {},
            'input': {},
            'package': {},
            'hashes': {},
            'defaults': {},
            'unset': [],
            'omit': [
                "userdifffname",
                "radialdifffname",
                "diff_file",
                "casename",
                "commands",
                "savefile",
            ],
        }

        # Assert input file exists before proceeding
        if filename is not None:
            if not exists(filename):
                raise ValueError('File {} does not exist!'.format(\
                    filename
                ))
        # Check whether UEDGE is installed and enter inplace mode w/
        # msg if not
        if (inplace is False) and (uedge_is_installed is False):
            print("No working UEDGE install found: only "+\
                "inplace-evaluation is possible.")
            print("Only UETOOLS HDF5 saves can be restored.")
            print("For more information, consult UETOOLS documentation.")
            inplace = True
        # Parse the casename
        if casename is None:
            try:
                casename = "/".join(".".join(filename.split(".")[:-1]).split("/")[-2:])
            except:
                casename = "/".join((".".join(getcwd().split(".")[:-1]).split("/")[-1],"Case()"))
        # Get the current UEDGE version
        try:
            uedge_ver = (
                packageobject("bbb").getpyobject("uedge_ver")[0].strip().decode("UTF-8")
            )
        except:
            uedge_ver = 'unknown'
        # Get the user-name
        try:
            user = getlogin()
        except OSError:
            # Can fail on e.g cluster nodes
            user = "unknown"
        # Get the hostname
        try:
            hostname = gethostname()
        except OSError:
            hostname = "unknown"
        # Get the name of the file being read
        if filename is None:
            location = getcwd()
        else:
            location = abspath(os.path.dirname(filename))
            filename = abspath(filename)
        # Check whether a separate diffusivity file is requested
        try:
            diffusivity_file = abspath(diffusivity_file)
        except:
            pass
        # Get the location of the case
        try:
            if exists('/'.join([location, savefile])):
                savefile = '/'.join([location, savefile])  
        except:
            # NOTE: Shoul probably raise an error here?
            pass
        # Store case data and information in dictionary
        self.info = {
            'casename': casename,
            'uetoolsversion': uetools.__version__,
            'uedge_ver':    uedge_ver,
            'pyver':    __version__,
            'user':     user,
            'hostname': hostname,
            'location': location,
            'inplace':  inplace,
            'restored_from_hdf5': False,
            'verbose':  verbose,
            'savefile': savefile,
            'filename': filename,
            'diffusivity_file': diffusivity_file,
            'aphdir': None,
            'apidir': None,
            'session_id': None,
            
        }
        # Link top-level classes to Case
        self.tools = Tools()
        self.search = Lookup()
        # Set up the functions to get/set data from Case and UEDGE
        # Reading from an HDF5 file
        if inplace:
            # Get the absolute path to the HDF5 file 
            self.info['filename'] = abspath(self.info['filename'])
            # Ensure the file exists
            if not exists(self.info['filename']):
                raise Exception("File {} not found. Aborting!".format(self.info['filename']))
            # Get the variables available in the HDF5 and link their locations
            self.load_inplace()
            self.getset = GetSetInplace(self)
        else:
            self.getset = GetSetMemory(self)

        # Link commands
        self.get = self.getset.get
        self.getue = self.getset.getue
        self.setue = self.getset.setue

        if not inplace:
            self.exmain_evals = self.getue("exmain_evals")
            # Assign mutex checks, unless at early UEDGE version
            try:
                self.exmain_evals = self.getue("exmain_evals")
                self.use_mutex = True
            except:
                print('Variable "exmain_evals" not found!')
                print('Using UEDGE version <7, deactivate mutex')
                self.use_mutex = False


        # Link all other Classes to Case
        self.tracker = Tracker(self)
        self.plot = Caseplot(self)
        self.postproc = PostProcessors(self)
        self.savefuncs = Save(self)
        self.save = self.savefuncs.save
        self.solver = Solver(self)
        self.populate = self.solver.populate
        self.utils = Misc(self)
        self.config = Config()
        self.convert = Convert(self)
        self.exmain = self.solver.exmain
        self.adas = ADAS(self)
#        self.radtransp = RadTransp(self)
        self.interpolate = Interpolate(self)
        self.about = AboutSetup(self)
        self.input = Input(self)
        # Set up paths from config file
        self.config.case(verbose=False)
        for key, value in self.config.configs.items():
            self.info[key] = value
        # Set up rate paths: Case-level paths take precedence over
        # config paths
        if aphdir is not None:
            self.info['aphdir'] = aphdir
        if apidir is not None:
            self.info['apidir'] = apidir

        # Perform additional operations, requiring linked packages
        if inplace is False:
            # Parse
            self.variables['input'] = self.tools.readyaml("{}/{}".format(
                uetools.__path__[0], "yamls/requiredvariables.yaml"
            ))
            # Read YAML to get variables to be read/stored/used
            if variableyamlfile is None:  # No YAML variable file requested
                if hasattr(self, "variableyamlfile"):
                    self.variables['input'].update(self.tools.readyaml(self.variableyamlfile))
                else:
                # Use default: find package location and package YAMLs
                    self.variables['input'].update(self.tools.readyaml("{}/{}".format(
                        uetools.__path__[0], "yamls/standardvariables.yaml"
                    )))
            else:  # YAML specified, use user input
                self.variables['input'].update(self.tools.readyaml(variableyamlfile))  # Read to memory
            if self.use_mutex is True:
                self.info['session_id'] = self.getue("max_session_id") + 1
                setattr(
                    packageobject("bbb"), "max_session_id", 
                    self.getue("max_session_id") + 1
                )
            if assign is True:
                self.assign()
            if self.info['filename'] is not None:
                self.restore_input(self.info['filename'], self.info['savefile'], 
                    restoresave=restoresave)
            else:
                self.reload()
                # TODO:
                if uedge_is_installed and not self.info['inplace']:
                    self.tracker.get_uevars()
                    if store_defaults:
                        # Track potential inputs too, just to be safe
                        for pkg in ['input', 'maybeinput']:
                            for key, _ in self.variables['hashes'][pkg].items():
                                self.variables['defaults'][key] = deepcopy(self.getue(key))
        self.plot = Caseplot(self)

    # NOTE: Update class data, or try reading from forthon first??
    def update(self, **kwargs):
        """Checks if UEDGE state has changed and updates as needed.

        Modifies
        --------
        Calls Case.reload, updates Case.vertices

        Returns
        -------
        """

        if self.use_mutex:
            if self.exmain_evals != self.getue("exmain_evals"):
                self.exmain_evals = self.getue("exmain_evals")
                if self.mutex() is False:
                    raise Exception("Case doesn't own UEDGE memory")
                self.reload()
                self.plot = Caseplot(self, rm=self.get("rm"), zm=self.get("zm"))
        else:
            self.reload()
            try:
                self.plot = Caseplot(self, rm=self.get("rm"), zm=self.get("zm"))
            except:
                pass



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
            setattr(packageobject("bbb"), "session_id", self.info['session_id'])
#        try:
#            # Restore input to UEDGE
#            # NOTE: variables not set maintain previous values. Reset
#            # all values before setting input?
#            print('assing')
#            self.setinput(readinput=False)
#        except:
#            pass
        try:
            if self.info['restored_from_hdf5'] is True:
                packageobject("grd").getpyobject("readgrid")(
                    self.getue("GridFileName"), self.variables['stored']["runid"].strip()
                )
        except:
            pass



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
            Dictionary of values, with keys specified in self.variables['input']

        Returns
        -------
        None

        """
        from numpy import ndarray, int64, float64

        # Check whether data is read into memory
        if self.info['inplace']:
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
                    if self.search.getpackage(group[-1], verbose=False) != None:
                        # Request to set array starting from index 0:
                        # just read the variable into memory
                        self.variables['stored'][group[-1]] = self.getue(group[-1])
                    elif isinstance(group[-1], int):
                        # Setting subarray, store variable
                        self.variables['stored'][group[-2]] = self.getue(group[-2])
                    else:
                        # List of variables, store each
                        for variable in dictobj:
                            self.variables['stored'][variable] = self.getue(variable)
                elif isinstance(group[-1], int):
                    if len(group) > 2 and isinstance(group[-2], int):
                        self.variables['stored'][group[-3]] = self.getue(group[-3])
                    else:
                        self.variables['stored'][group[-2]] = self.getue(group[-2])
                elif isinstance(dictobj, bool):
                    # TODO: Now assumed only Falses set, which do nothing
                    # In the future, we might include Trues on keywords.
                    # Such behavior goes here
                    pass
                elif isinstance(dictobj, (int, float, int64, float64)):
                    self.variables['stored'][group[-1]] = self.getue(group[-1])
                elif isinstance(dictobj, (bytes, str)):
                    try:
                        self.variables['stored'][group[-1]] = self.getue(group[-1])
                    except:
                        pass
                else:
                    self.variables['unset'].append([group, dictobj])
            else:
                for key, value in dictobj.items():
                    recursivereload(value, group + [key])
        # Pop out any custom commands, as these cannot be reloaded (not vars)
        try:
            commands = self.variables['input']["setup"].pop("commands")
        except:
            pass
        # Reload the variables recurively
        if group is None:
            recursivereload(self.variables['input'])
        else:
            recursivereload(self.variables['input'][group], [group])
        # If there were any custom commands, put them back where they belong
        try:
            self.variables['input']["setup"]["commands"] = commands
        except:
            pass
        # Update the dict containing the package containing each variable
        for variable in self.variables['stored'].keys():
            if variable not in self.variables['package']:
                self.variables['package'][variable] = self.search.getpackage(
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
        from h5py import Group, File
        if fileobj is None:
            with File(self.info['filename'], 'r') as f:
                for subgroup, data in f.items():
                    self.load_inplace(data, group + [subgroup])
        elif isinstance(fileobj, File):
            for subgroup, data in fileobj.items():
                self.load_inplace(data, group + [subgroup])
        elif isinstance(fileobj, Group):
            for subgroup, data in fileobj.items():
                self.load_inplace(data, group + [subgroup])
        else:
            self.variables['stored'][fileobj.name.split("/")[-1]] = fileobj.name
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
        if self.info['session_id'] == self.getue("session_id"):
            return True
        else:
            if silent is False:
                print(
                    "Mutex error! Object run-ID is {}, UEDGE run-ID "
                    "is {}. Aborting.".format(self.info['session_id'], self.getue("session_id"))
                )
            return False

    def restore_input(self, inputfname=None, savefile=None, 
        populate=True, restoresave=True, **kwargs):
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

        self.input.read(inputfname, savefile=savefile, restoresave=restoresave, 
                **kwargs
        )
        if (restoresave is True) and (populate is True):
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
        self.savefuncs.load_state(savefile, **kwargs)
        self.populate(**kwargs)

    def add_spectrometer(self, specname=None, **kwargs):
        from uetools.UeDiagnostics import Spectrometer
        try:
            self.diagnostics
        except:
            self.diagnostics = {}
        if specname is None:
            specname = 'spec{}'.format(len(self.diagnostics)+1)
        self.diagnostics[specname] = Spectrometer(self, **kwargs)
        return self.diagnostics[specname]


    def dashboard(self):
        """ Opens a Dashboard for Self """
        from uetools import StandaloneDashboard 
        from PyQt5.QtWidgets import QApplication
        import sys
        app = QApplication([])
        win = StandaloneDashboard(self)    
        win.show()
        app.exec_()

        
class GetSetMemory:
    def __init__(self, case):
        self.update = case.update
        self.variables = case.variables
        self.search = Lookup()
        self.mutex = case.mutex


    def get(self, variable, s=None, **kwargs):
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
            retvar = self.variables['stored'][variable]
        except:
            retvar = self.getue(variable)
        # Check the size of the array, and return index if multi-species array
        if isinstance(retvar, (ndarray, list)):
            if len(retvar.shape) == 3:
                if s is not None:
                    retvar = retvar[:, :, s]
        return retvar

    def setue(self, variable, value, **kwargs):
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
            package = self.variables['package'][variable]
        except:
            package = self.search.getpackage(variable, verbose=False)
        if self.mutex():
            try:
                setattr(packageobject(package), variable, value)
            except Exception as e:
                raise KeyError("{} could not be set: {}".format(variable, e))

    def getue(self, variable, s=None, cp=True, **kwargs):
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
            package = self.variables['package'][variable]
        except:
            package = self.search.getpackage(variable, verbose=False)

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


class GetSetInplace:
    def __init__(self, case):
        self.info = case.info
        self.variables = case.variables

    def get(self, variable, s=None, verbose=True, **kwargs):
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
        
        try:
            with File(self.info['filename'], 'r') as f:
                retvar = f[self.variables['stored'][variable]][()]
        except:
            if verbose:
                print("{} not found in {}".format(variable, self.info['filename']))
            return

        if isinstance(retvar, (ndarray, list)):
            if len(retvar.shape) == 3:
                if s is not None:
                    retvar = retvar[:, :, s]
        return retvar

    def getue(self, *args, **kwargs):
        """Placeholder to avoid getting/setting when reading inplace."""
        raise Exception("Cannot get UEDGE values when reading from "+\
                            "HDF5 file")

    def setue(self, *args, **kwargs):
        """Placeholder to avoid getting/setting when reading inplace."""
        raise Exception("Cannot set UEDGE values when reading from "+\
                            "HDF5 file")


