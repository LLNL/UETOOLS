from Forthon import packageobject
from .CasePlot import Caseplot
from .Solver import Solver
from .Track import Tracker
from uetools.UeLookup.Lookup import Lookup
from uetools.UePostproc.Postproc import PostProcessors
from uetools.UeConfig.Config import Config
from uedge import bbb, com, aph, api, svr

# TODO: Where/how to define impurity files and paths?
# TODO: Update data after reading save file
# TODO: make yaml read/write/get/set case-insensitive
# TODO: Consider compression of data
# TODO: implement divergence plotting/calculation
# TODO: Unify all data to be stored in the same dictionary?

class Case(Caseplot, Solver, Lookup, PostProcessors):
    """ UEDGE Case container object.

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
        switch defining wheter to read data into memory or read from 
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
    grouplist : dict 
        lookup dictionary listing groups associated with 
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
    createhelperdicts():
        initializes packagelist and grouplist
    closehdf5():
        closes hdf5case File object
    createvarsdict():
        initializes vars from varinput
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

    def __init__(self, casefname=None, inplace=False, variableyamlfile = None,
        casename=None, assign=True, **kwargs):
        """ Initializes the UeCase object.

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
        
        conf = Config()
        if conf.configured is False:
            return
        # TODO: add label attribute
        # Read from uedge.label and strip
        # Use to restore and/or save files
        # Initialize class attributes
        # Checksum whether to update data
        self.session_id = self.getue('max_session_id') + 1
        setattr(packageobject('bbb'), 'max_session_id', 
            self.getue('max_session_id')+1)
        self.exmain_evals = self.getue('exmain_evals')
        if assign is True:
            self.assign()

        self.casename = casename
        self.uetoolsversion = '0.1' # UEtools version
        self.allocate = packageobject('bbb').getpyobject('allocate')
        self.casefname = casefname
        self.inplace = inplace
        for selfdict in ['vars', 'varinput', 'packagelist', 'grouplist',
            'setup', 'commands']:
            setattr(self, selfdict, dict())
            #setattr(self, selfdict, self.createdict())
        self.userdifffname = None
        self.radialdifffname = None
        self.hdf5case = None
        # Read YAML to get variables to be read/stored/used
        if self.casefname is None: # No saved case restored, read from memory
            if variableyamlfile is None: # No YAML variable file requested
                # Use default: find package location and package YAMLs
                variableyamlpath = '{}/{}'.format(uetools.__path__[0], 
                    'yamls/standardvariables.yml') 
            else: # YAML specified, use user input
                variableyamlpath = '{}/{}'.format(path, variableyamlfile)
            self.varinput = self.readyaml(variableyamlpath) # Read to memory
            # Create dictionaries based on YAML input
            self.createhelperdicts()
            self.createvarsdict()
        # Read data from stored file
        else: # Restore data from previous save file
            self.hdf5case = self.openhdf5(self.casefname, 'r')
            if self.inplace: # Access HDF5 file for data
                self.createhelperdicts()
            else: # Read HDF5 data into memory
                self.createvarsdictfromhdf5()
                self.closehdf5()
        # Initialize parent classes
        super(Case, self).__init__()
        return

    # NOTE: Update class data, or try reading from forthon first??
    def update(self, **kwargs):
        ''' Checks if UEDGE state has changed and updates if needed '''

        if self.exmain_evals != self.getue('exmain_evals'):
            if self.mutex() is False:
                return
            self.exmain_evals = self.getue('exmain_evals')
            self.reload()
            self.vertices = self.createpolycollection(self.get('rm'), 
                self.get('zm'))
        
    def get(self, variable, s=None, **kwargs):
        """ Returns variable 
        
        Method assumes unique variable names across all packages and 
        groups. Returns data from vars is inplace=False, from HDF5 
        specified at creation if inplace=True.

        First, checks if UEDGE solution has changed from previous step.
        If it has, updates UeCase data. Then, looks for data in UeCase 
        and returns it if found. If not found, the data is accessed from
        Forhton memory and returned.

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
        self.update() # Update results from UEDGE if they have changed
        # Switch to asses where to access data from
        if self.inplace is True:
            retvar = self.hdf5case[self.grouplist[variable]]\
                [self.packagelist[variable]][variable][()]
        else:
            try:
                retvar = self.vars[self.packagelist[variable]][variable]
            except:
                retvar = self.getue(variable)
        # Check the size of the array, and return index if multi-species array
        if isinstance(retvar, (ndarray, list)):
            if len(retvar.shape) == 3:
                if s is not None:
                    retvar = retvar[:, :, s]
        return retvar

    def assign(self, **kwargs):
        ''' Assigns the UEDGE session to this object '''
        setattr(packageobject('bbb'), 'session_id', self.session_id)
        try:
            self.setinput(readyaml=False) 
        except:
            pass

    def setue(self, variable, value, idx=None, **kwargs):
        """ Sets the Forthon variable in package to data 

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
            package = self.getpackage(variable)
        if self.mutex():
            setattr(packageobject(package), variable, value)

    def getue(self, variable, s=None, **kwargs):
        """ Retireves data from UEDGE variable in package 

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
            package = self.getpackage(variable)

        retvar = deepcopy(packageobject(package).getpyobject(variable))

        if isinstance(retvar, (ndarray, list)):
            if len(retvar.shape) == 3:
                if s is not None:
                    retvar = retvar[:, :, s]
        return retvar

    def reload(self, group=None, **kwargs):
        """ Reloads variables from UEDGE to UeCase

        Keyword arguments
        ------------
        group : str (default = None)
            group specifier to reload. If None, reloads all
            variables in vars
        """
        from copy import deepcopy

        if self.mutex is False:
            return
        try:
            casename = self.varinput['setup'].pop('casename')
        except:
            pass
        try:
            commands = self.varinput['setup'].pop('commands')
        except:
            pass

        def recursivereload(dictobj = None, group=[]):
            if dictobj is None:
                dictobj = self.varinput
            if not isinstance(dictobj, dict):
                if isinstance(dictobj, list):
                    for variable in dictobj:
                        try:
                            self.vars[group[-1]]
                        except:
                            self.vars[group[-1]] = dict()
                        self.vars[group[-1]][variable] = \
                            self.getue(variable)
                elif isinstance(group[-1], int):
                    # Store the whole variable
                    try:
                        self.vars[group[-3]]
                    except:
                        self.vars[group[-3]] = dict()
                    self.vars[group[-3]][group[-2]] = \
                        self.getue(group[-2])
            else:
                for key, value in dictobj.items():
                    dictobj = recursivereload( value, group + [key])
                return dictobj 

        # Check whether data is read into memory
        if self.inplace:
            print('Cannot reload directly to HDF5 file with option "inplace".')
            print('Aborting')
            return
        if group is None:
            recursivereload()
        else:
            recursivereload(self.varinput[group], [group])
        try:
            self.varinput['setup']['casename'] = casename
        except:
            pass
        try:
            self.varinput['setup']['commands'] = self.commands
        except:
            pass

    def reinitializeplot(self, **kwargs):
        """ Reinitializes the data of Subclasses for plotting """
        # TODO: Return from UEDGE or variables?
        self.vertices = self.createpolycollection(self.getue('rm'), 
            self.getue('zm'))

    def createvarsdict(self, **kwargs):
        """ Rearranges self.vars into structure [package][variable] """
        # NOTE: Check before loop - fewer conditionals, more repeated code
        # NOTE: Check in loop - less code, what kind of slowdown?
        
        if self.casefname is None:
            for group, packages in self.varinput.items():
                for package, variables in packages.items():
                    for variable in variables:
                        try:
                            self.vars[package]
                        except:
                            self.vars[package] = dict()
                        self.vars[package][variable] = \
                            self.getue(variable)
        else:
            for group, packages in self.varinput.items():
                for package, variables in packages.items():
                    for variable in variables:
                        try:
                            self.vars[package]
                        except:
                            self.vars[package] = dict()
                        self.vars[package][variable] = \
                            self.varinput[group][package][variable]

    def createvarsdictfromhdf5(self, **kwargs):
        ''' Creates a dictionary of variables in YAMLs '''
        self.varinput = self.gethdf5data(self.hdf5case)
        self.vars = {}
        def  recursivecreatevars(dictobj = None, group=[], keylist = None):
            from numpy import ndarray
            if not isinstance(dictobj, dict):
                try:
                    self.vars[group[-2]]
                except:
                    self.vars[group[-2]] = {}
                self.vars[group[-2]][group[-1]] = dictobj
            else:
                keylist = list(dictobj.keys())
                for key, value in dictobj.items():
                    if not isinstance(value, dict):
                        if group[0] != 'setup':
                            vardict = self.varinput
                            for g in group[:-1]:
                                if g not in list(vardict.keys()):
                                    vardict[g] = {}
                                vardict = vardict[g]
                            vardict[group[-1]] = keylist
                    dictobj = recursivecreatevars( value, group + [key], keylist)
                return dictobj 
        recursivecreatevars(self.varinput)
        # Try setting input
        self.setinput(readyaml=False)
        self.createhelperdicts()

    def createhelperdicts(self, **kwargs):
        """ Initializes grouplist and packagelist """
        if self.inplace is True: # If no save to be read, use local data
            datasource = self.hdf5case
        else: # Get data from HDF5
            datasource = self.varinput
        # NOTE: assumption is all variables have unique names: no 
        # duplicates names in different packages
        for group, packages in datasource.items():
            for package, variables in packages.items():
                for variable in variables:
                    # TODO: Add error/warning for duplicates?
                    self.packagelist[variable] = package
                    self.grouplist[variable] = group

    def readyaml(self, fname, **kwargs):
        """ Reads a YAML file and returns a nested dict 

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
        """ Opens HDF5 file and returns File object

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
        from h5py import File
        try:
            return File(fname, operation)
        except:
            print('File "{}" not found!'.format(fname))

    def closehdf5(self, **kwargs):
        """ Closes UeCase file hdf5case  that is being read """
        try:
            self.hdf5case.close()
        except:
            print('No HDF5 file open')
        

            
    def setinput(self, setupfile=None, restore=True, allocate=True, 
        savefname=None, readyaml=True, restoresave=False, **kwargs):
        ''' Sets all UEDGE variables from Case input '''
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
        if self.mutex() is False:
            return
        if readyaml is True:
            if setupfile is None:
                setupfile = '{}.yaml'.format(self.casename)
            self.varinput['setup'] = self.readyaml(setupfile)


        # Pop out groups that cannot be parsed by default
        try:
            self.commands = self.varinput['setup'].pop('commands')
        except:
            pass
        # TODO: tidy up casename definition
        try:
            self.casename = self.varinput['setup'].pop('casename')
        except:
            pass
        if self.casename is None:
            self.casename = casename

        # TODO: Find a way to catch user-specified and radially varying
        #       diffusive coefficients when reading from file: userdifffname
        #       and radialdifffname attributes no available!
        def setinputrecursive(dictobj, group=[], allocate=True):
            # NOTE: Something in this function is SLOOOW
            if not isinstance(dictobj, dict):
                # Skip UeCase-unique parameters
                if group[-1] not in ['userdifffname', 'radialdifffname']:
                    # Circumvent the padding with nulls for strings
                    try:
                        dictobj = dictobj.ljust(len(self.getue(
                            group[-2])[group[-1]]))
                    except:
                        pass
                    # Set all other parameters
                    if isinstance(group[-1], int): # Set index-by-index
                        datalist = self.getue(group[-2])
                        datalist[group[-1]] = dictobj
                        try:
                            self.setue(group[-2], datalist)
                        except:
                            if allocate is True:
                                self.allocate()
                            self.setue(group[-2], datalist)
                    else: # Set whole array
                        try:
                            self.setue(group[-1], dictobj)
                        except:
                            if allocate is True:
                                self.allocate()
                            self.setue(group[-1], dictobj)
                else: # Set calls to restore diffusivities
                    setattr(self, group[-1], dictobj)
            else:
                for key, value in dictobj.items():
                    dictobj = setinputrecursive( value, group + [key], 
                        allocate)
                return dictobj 
        # Set group order to assure proper allocation and avoid warnings
        neworder = OrderedDict()
        setupkeys = list(self.varinput['setup'].keys())
        for key in ['grid', 'boundaryconditions', 'allocate', 'diffusivities']:
            setupkeys.remove(key)
        for key in ['grid', 'boundaryconditions', 'allocate']:
            setupkeys.insert(0, key)
        setupkeys.append('diffusivities')
        for key in setupkeys:
            neworder[key] = self.varinput['setup'][key]
        self.varinput['setup'] = neworder 
        # Set UEDGE parameters to input if restore is True
        if restore is True:
            setinputrecursive(self.varinput['setup'], allocate=allocate)
            # See if diffusivities unset despite being user-defined
            # If yes, try looking for them in the case being restored
            if (self.varinput['setup']['diffusivities']['isbohmcalc'] in [0, 1]) \
                and (self.userdifffname is None):
                    self.userdifffname = self.casefname
            elif (self.varinput['setup']['diffusivities']['isbohmcalc'] == 2) \
                and (self.radialdifffname is None):
                    self.radialdifffname = self.casefname

            if allocate is True:
                self.allocate()
            if self.userdifffname: 
                self.setuserdiff(self.userdifffname)
            if self.radialdifffname: 
                self.setradialdiff(self.radialdifffname)
        # TODO: Can this mess be made somehow prettier?
        if restoresave is True:
            self.restoresave(savefname, **kwargs)
        self.reload()
        for command in self.commands:
            exec(command)
        # TODO: don't put this back in, do something else with this?
        try:
            self.varinput['setup']['casename'] = self.casename
        except:
            pass
        try:
            self.varinput['setup']['commands'] = self.commands
        except:
            pass

    def setuserdiff(self, difffname, **kwargs):
        """ Sets user-defined diffusivities

        Arguments
        ------------
        diffname : str
            HDF5 file from where to read 'diffusivities'/'bbb'/values
         """
        # TODO: replace with save-group function call?

        if self.mutex() is False:
            return

        try:
            difffile = self.openhdf5(difffname, 'r')
        except:
            difffile = self.openhdf5(self.casefname, 'r')
            
        for variable in ['dif_use', 'kye_use', 'kyi_use', 'tray_use']:
            self.setue(variable, 
                difffile['diffusivities']['bbb'][variable][()])
        difffile.close()
        return

    def mutex(self, silent=False, **kwargs):
        ''' Returns bool whether case assigned to current UEDGE session

        Keyword parameters
        ------------------
        silent : boolean (default : False)
            Switch whether to issue mutex warning or not
        '''
        if self.session_id == self.getue('session_id'):
            return True
        else:
            if silent is False:
                print('Mutex error! Object run-ID is {}, UEDGE run-ID is {}. Aborting.'.format\
                    (self.session_id, self.getue('session_id')))
            return False

    def setradialdiff(self, difffname, **kwargs):
        """ Sets radially varying diffusivities 

        Arguments
        ------------
        diffname : str
            HDF5 file from where to read 'diffusivities'/'bbb'/values
        """

        if self.mutex() is False:
            return

        try:
            difffile = self.openhdf5(difffname, 'r')
        except:
            difffile = self.openhdf5(self.casefname, 'r')
        for variable in ['difniv', 'kyev', 'kyiv', 'travisv']:
            self.setue(variable, 
                difffile['diffusivities']['bbb'][variable][()])
            difffile.close()
        return

    def setgroup(self, group, savefname = None, **kwargs):
        """ Sets UEDGE variables listed in group to UeCase values

        Arguments
        ------------
        group : str
            group to be set in UEDGE from UeCase values
        
        Keyword arguments
        ------------
        savefname : str (default = None)
            name of HDF5 file from where to read data. If None,
            data is read from UeCase object
        """

        if self.mutex() is False:
            return

        # NOTE: Is this function obsolete? Can it be superseded
        #       by some more general function
        # TODO: Add option to read from other source file
        if isinstance(savefname, str):
            datafile = self.openhdf5(savefname, 'r')
            datasource = datafile['restore']
            for package, variables in datasource.items():
                for variable, data in variables.items():
                    self.setue(variable, data[()])
            datafile.close()
        elif self.inplace: # Read from file
            datasource = self.hdf5case['restore']
            for package, variables in datasource.items():
                for variable, data in variables.items():
                    self.setue(variable, data[()])
        else:
            datasource = self.vars['restore']
            for package, variables in datasource.items():
                for variable, data in variables:
                    self.setue(variable, data)
        # Restore every variable classed as restore
        return



    def gethdf5data(self, fileobj, **kwargs):
        """ Returns all data from fileobj in nested dict """
 
        # Create new home for variable
        savedict = dict()       
        # Return stored data if we are at the bottom of loop
        try:
            return fileobj[()]
        # If not, recurse deeper
        except:
            # First, make sure we have an open object
            try:
                save = self.openhdf5(fileobj, 'r')
            except:
                save = fileobj
            # Next, loop through the values
            for key, value in save.items():
                # Test whether element is set
                try:
                    key = int(key)
                except:
                    pass
                # Recurse deeper into function
                savedict[key] = self.gethdf5data(value)
            # Return constructed dictionary
            return savedict


    def restoresave(self, savefname=None, **kwargs):
        """ Restores a saved solution 
        
        Keyword arguments
        ------------
        savefname : str (default = None)
            HDF5 file to read stored solution from. If None, solution
            is read from UeCase object
        **kwargs
            passed to setgroup
        """
        if self.mutex() is False:
            return

        if savefname is None:
            savefname = '{}.hdf5'.format(self.casename)
        savefile = self.openhdf5(savefname, 'r')
        # Try reading new, subdivided save file
        try:
            # Don't override user-specified name for case by reading from file
            if casefname is None:
                self.casename = savefile.attrs['casename']
        except:
            pass
        try:
            for group, variables in savefile['restore'].items():
                for variable, value in variables.items():
                    self.setue(variable, value[()])
                    self.vars[group][variable] = value[()]
        # If not, try reading old-style save file
        except:
            for group, variables in self.varinput['restore'].items():
                for variable in variables:
                    self.setue(variable, savefile[group][variable][()])
                    self.vars[group][variable] = savefile[group][variable][()]
        return

    def populate(self, silent=False, **kwargs):
        """ Populates all UEDGE arrays by evaluating static 'time-step' """
        from copy import deepcopy

        if self.mutex() is False:
            return

        if silent is True:
            self.setue('iprint', 0)
        issfon = deepcopy(self.getue('issfon'))
        ftol = deepcopy(self.getue('ftol'))
        self.setue('issfon', 0)
        self.setue('ftol', 1e20)
        self.exmain()
        self.setue('issfon', issfon)
        self.setue('ftol', ftol)
        self.update()
        if silent is True:
            fnrm = sum(self.getue('yldot')**2)**0.5
            prtstr = '\n*** UEDGE arrays populated: {} ***'
            if fnrm < 10:
                print(prtstr.format('Case appears converged'))
                print('fnrm without preconditioning: {:.2e}\n'.format(fnrm))
            elif fnrm < 100:
                print(prtstr.format('Warning, case may noy be fully converged')) 
                print('fnrm without preconditioning: {:.1f}\n'.format(fnrm))
            else:
                print(prtstr.format('WARNING, case NOT converged'))
                print('fnrm without preconditioning: {:.2e}\n'.format(fnrm))
            self.setue('iprint', 1)

    def restore(self, inputfname=None, savefname=None, **kwargs):
        """ Restores a full case into memory and object"""
        if self.mutex() is False:
            return
        self.setinput(inputfname, savefname=savefname, restoresave=True,
            **kwargs)
        self.populate(silent=True, **kwargs)

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
        try:
            savefile.create_dataset(variable, data=data)
        except:
            # Catch exception of unset data and store False
            # TODO: Omit or save False?
            # TODO: Verbose prompt or not?
            savefile.create_dataset(variable, data=False)
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
        if not isinstance(saveobj, dict):
            if isinstance(saveobj, list):
                for variable in saveobj:
                    self.savevar(savefile, group, variable, 
                        self.getue(variable))
            elif group[-1] not in ['userdifffname', 'radialdifffname']:
                variable = group.pop(-1)
                self.savevar(savefile, group, variable, 
                    self.getue(variable))
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
            import __version__ as pyv
            pyver = pyv.__version__
        except:
            pyver = __version__
        savefile.attrs['casename'] = self.casename
        savefile.attrs['time'] = time()
        savefile.attrs['ctime'] = ctime()
        savefile.attrs['code'] = 'UEDGE'
        savefile.attrs['ver'] = packageobject('bbb').getpyobject('uedge_ver')
        savefile.attrs['pyver'] = pyver
        # TODO: add user, hostname, aphdir, apidir mcfilename, aphfname

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
        # TODO: always overwrite or have a switch?
        # Open file to save to
        # TODO: add commands to save
        # TODO: add savefname to save
        savefile = self.openhdf5(savefname, 'w'*(append is False) \
            + 'a'*(append is True)) 
        self.savemetadata(savefile)
        if group is None:
            self.recursivesave(savefile, self.varinput)
        else:
            self.recursivesave(savefile, self.varinput[group], [group])
        savefile.close()
 
    '''
    def readgridue(self, gridue='gridue'):
        """ Reads a gridue file and primes UEDGE for reading grid

        Keyword arguments
        ------------
        gridue : str (default = 'gridue')
            path to/name of gridue file
        """
        # NOTE: interim solution for testing!
        self.setue('bbb', 'mhdgeo', 1)
        self.setue('com', 'geometry', 'snull')
        self.setue('bbb', 'gengrid', 0)
        self.getue('grd', 'readgrid')(gridue, '')
        self.reload('grid')
        # TODO: put read gridue here

    def readgridhdf5(self, savefname):
        """ Reads a grid HDF5 file and stores to memory

        Arguments
        ------------
        savefname : str
            path to/name of HDF5 file to read grid from
        """
        # Read the grid dimension to allocate
        # TODO: check source of grid data
        gridfile = self.openhdf5(savefname, 'r') 
        self.setue('com', 'nx', gridfile['grid']['com']['nx'][()])
        self.setue('com', 'nx', gridfile['grid']['com']['ny'][()])
        self.setue('com', 'nx', gridfile['grid']['com']['iysptrx'][()])
        self.setue('com', 'nx', gridfile['grid']['com']['ixpt1'][()])
        self.setue('com', 'nx', gridfile['grid']['com']['ixpt2'][()])
        # Drive the grid restoration scripts from here
        self.readpackage('grid')   
        return 

    def createdict(self):
        """ Creates dict that assigns unspecified keys on call. """
        from collections import defaultdict
        return defaultdict(lambda: defaultdict(dict)) 

    def load(self, casefname=None):
        """ Loads data from HDF5 to UeCase
        
        Keyword arguments
        ------------
        casefname : str (default = None)
            HDF5 file to restore 
        
        Returns
        ------------
        None
        """
        # TODO: verify this function is working as intended!

        if casefname is None:
            casefname = self.casefname
        casefile = self.openhdf5(casefname, 'r')
        for group, packages in  casefile.items():
            for package, variables in packages.items():
                for variable, data in variables.items():
                    self.varinput[group][package][variable] = data[()]
        # TODO: print relevant file metadata
        self.createhelperdicts()
        self.createvarsdict()
        self.reinitializeplot()
        casefile.close()
        return

    def setdimensions(self):
        """ Allocates UEDGE arrays to right dimensions """
        self.setue('com', 'ngsp', self.vars('com', 'ngsp'))
        self.setue('com', 'nhsp', self.vars('com', 'ngsp'))
        self.setue('com', 'ngsp', self.vars('com', 'nzsp'))
        self.setue('com', 'nx', self.vars('com', 'nx'))
        self.setue('com', 'ny', self.vars('com', 'ny'))
'''
