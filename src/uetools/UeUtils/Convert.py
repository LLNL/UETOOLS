from .Lookup import Lookup

class Convert:
    def __init__(self, case):
        self.getue = case.getue
        self.setue = case.setue
        self.reload = case.reload
        self.varinput = case.varinput
        self.get = case.get
        self.save = case.save
        self.info = case.info
        self.defaults = case.defaults
        self.populate = case.populate
        self.record_changes = case.tracker.record_changes
        self.get_bottomkeys = case.utils.get_bottomkeys
        self.getpackage = Lookup().getpackage

    def write_py(self, fname):
        """ Writes a standalone Python input file for the Case
        """
        try:
            from Forthon import package, packageobject
        except:
            pass

        def recursive_lineread(dictobj, lines=None, fails=None):
            """ Recusively parses setup variables to input lines """
            from Forthon import getpackage
            # TODO: Implement checks on string lengths & shapes
            #       2D arrays are being set by 1D arrays work for YAMLs
            #       but not in python...
            if lines is None:
                lines = []
            if fails is None:
                fails = {}
            # Loop through all the entries at the current level
            for key, value in dictobj.items():
                # Detect and omit keyword inputs: treated separately
                if key in [
                            'casename', 
                            'savefile',
                            'commands', 
                            'userdifffname',
                            'radialdifffname',
                            'diff_file'
                ]:
                    # Pass using fails
                    fails[key] = value
                    continue
                # Catch any variables being set
                elif self.getpackage(key) is not None:
                    # Check if indices are being set separately
                    if not isinstance(value, dict):
                        # Check if whole array is being set 
                        if isinstance(value, list):
                            # Transpose ordering for setting 2D arrays
                            # (primarily albedolb and albedorb)
                            if len(self.getue(key).shape)>1:
                                valbuff = []
                                # Wrap in nested list to transpose
                                for v in value:
                                    valbuff.append([v])
                                value=valbuff
                            # Set the first N elements corresponding to
                            # the list size
                            lines.append('{}.{}[:{}]={}'.format(
                                            self.getpackage(key), key, 
                                            len(value), value
                                        ))
                        # Omitted variable - raise/gather them?
                        elif value is False:
                            continue
                        # Variable a string, include "
                        elif isinstance(value, str):
                            lines.append('{}.{}="{}"'.format(
                                            self.getpackage(key), key, value
                                        ))
                        # Any other variable is a direct int/float, just set it
                        else:
                            line = '{}.{}={}'.format(
                                            self.getpackage(key), key, value)
                            # NOTE: Lazy protection againstduplicating input
                            # definitions for nested input decks. I have no
                            # clue why nested decks end up in this loop 
                            # twice... This is a makeshift solution
                            if line not in lines:
                                lines.append(line)
                    # Setting a non-starting index defined on separate line
                    # Create placeholder and fill later on
                    elif len(value) > 1:
                        lines, fails = recursive_lineread(value, lines, fails)
                    else:
                        lines.append('{}.{}[{}]={}'.format(
                                            self.getpackage(key), key, 
                                            '{}', '{}'
                                        ))
                # If the variable is a string but not a UEDGE var,
                # assume it is a 'header' and set as comment in input file
                elif isinstance(key, str):
                    lines.append('\n# ==== {} ===='.format(key.upper()))
                # If the variable is an integer, it's a non-zero starting
                # index for a variable defined on the previous line
                elif isinstance(key, int):
                    # If several indices are set starting at a non-zero
                    # index, use the span as index
                    if isinstance(value, list):
                        index = '{}:{}'.format(key, key+len(value)) 
                    # Otherwise, it's a single index being set
                    else:
                        index = key
                    # If the value is a string, include "
                    if isinstance(value, str):
                        value = '"{}"'.format(value)
                    # Backfill the previous line with the correct values
                    if isinstance(value, dict):
                        subindices = [str(index)]
                        while isinstance(value, dict):
                            for index, value in value.items():
                                subindices.append(str(index))
                        lines[-1] = lines[-1].format(','.join(subindices), value)
                    else:
                        lines[-1] = lines[-1].format(index, value)
                # Catch anything falling through 
                else:
                    fails[key] = value
                # The first line has been parsed. If it was a header, the
                # subdictionaries should be set. Call recursively
                if isinstance(value, dict):
                    # TODO: move to be called after header?
                    # Should be logically equivalent
                    lines, fails = recursive_lineread(value, lines, fails)
            # Pass the lines and fails up to the recursive parent
            return lines, fails
        # Get the lines parsed and those that could not be parsed
        lines, fails = recursive_lineread(self.varinput['setup'])
        # Pop out values defined in Case object
        try:
            fails.pop('casename')
        except:
            pass
        try:
            fails.pop('savefile')
        except:
            pass
        # Open file for writing
        with open(fname, 'w') as f:
            # Import required packages
            f.write('from uedge import *\n')
            f.write('from h5py import File\n')
            f.write('from uedge.hdf5 import hdf5_restore\n\n')
            # Set casename variable
            f.write('casename = "{}"\n'.format(self.info['casename']))
            # Set paths to rate files as defined in configuration file
            f.write('aph.aphdir = "{}"\n'.format(self.get('aphdir')[0].decode('UTF-8').strip()))
            f.write('api.apidir = "{}"\n'.format(self.get('apidir')[0].decode('UTF-8').strip()))
            # WRITE LINES
            for line in lines:
                f.write(line)
                f.write('\n')
            # Allocate to make space for restore variables
            f.write('\nbbb.allocate()\n'.format(self.info['savefile']))
            # Restore solution: use inplace version to ensure compatibility
            # with both native UEDGE saves and UETOOLS saves
            f.write('\n# ==== RESTORE SOLUTION ====\n')
            f.write('with File("{}") as f:\n'.format(self.info['savefile']))
            f.write('    try:\n')
            f.write('        savegroup = f["restore/bbb"]\n')
            f.write('    except:\n')
            f.write('        savegroup = f["bbb"]\n')
            f.write('    for var in ["ngs", "nis", "phis", "tes", '
                        '"tgs", "tis", "ups"]:\n')
            f.write('        setattr(bbb, var, savegroup[var][()])\n')
            # Now, add any manual commands there might be
            if 'commands' in fails:
                f.write('\n# ==== USER-SPECIFIED COMMANDS ====\n')
                # Pop to remove from fails as they are being set
                for command in fails.pop('commands'):
                    # Populate is a Case-call: replicate it using UEDGE
                    if 'self.populate' in command:
                        f.write('# Populate UEDGE arrays\n')
                        f.write(
                            'bbb.issfon=0\n'
                            'bbb.ftol=1e20\n'
                            'bbb.exmain()\n'
                            'bbb.issfon=1\n'
                            'bbb.ftol=1e-8\n'
                        )
                    # Write commands to file
                    else:
                        f.write(command)
                        f.write('\n')
            # If a file is used to set the radial transport coefficients,
            # manually for the whole domain, read and restore them here
            if self.info['diffusivity_file'] is not None:
                file = self.info['diffusivity_file']
                f.write('\n# ==== SET USER-SPECIFIED DIFFUSIVITIES ====\n')
                f.write('with File("{}", "r") as f:\n'.format(file))
                f.write('    try:\n')
                f.write('        bbb.difniv=f["diffusivities/bbb/difniv"][()]\n')
                f.write('    except:\n')
                f.write('        pass\n')
                f.write('    try:\n')
                f.write('        bbb.kyev=f["diffusivities/bbb/kyev"][()]\n')
                f.write('    except:\n')
                f.write('        pass\n')
                f.write('    try:\n')
                f.write('        bbb.kyiv=f["diffusivities/bbb/kyiv"][()]\n')
                f.write('    except:\n')
                f.write('        pass\n')
                f.write('    try:\n')
                f.write('        bbb.travisv=f["diffusivities/bbb/travisv"][()]\n')
                f.write('    except:\n')
                f.write('        pass\n')
                f.write('    try:\n')
                f.write('        bbb.dif_use=f["diffusivities/bbb/dif_use"][()]\n')
                f.write('    except:\n')
                f.write('        pass\n')
                f.write('    try:\n')
                f.write('        bbb.kye_use=f["diffusivities/bbb/kye_use"][()]\n')
                f.write('    except:\n')
                f.write('        pass\n')
                f.write('    try:\n')
                f.write('        bbb.kyi_use=f["diffusivities/bbb/kyi_use"][()]\n')
                f.write('    except:\n')
                f.write('        pass\n')
                f.write('    try:\n')
                f.write('        bbb.tray_use=f["diffusivities/bbb/tray_use"][()]\n')
                f.write('    except:\n')
                f.write('        pass\n')
                f.write('    try:\n')
                f.write('        bbb.fcdif=f["diffusivities/bbb/fcdif"][()]\n')
                f.write('    except:\n')
                f.write('        pass\n')
        # Check whether there still are lines that could not be parsed 
        # into the input, and output them to the promt
        if len(fails)>0:
            print('WARNING! Some variables could not be parsed into input!')
            print('These variables/groups were:')
            for key, value in fails.items():
                print('    -{}'.format(key))


    def write_yaml(self, fname):
        from yaml import dump
        with open(fname, 'w') as f:
            dump(self.varinput['setup'], f, sort_keys=False,default_style=None, default_flow_style=False)
        # TODO: add capability to also write autodetected changes



    def write_py2yaml(self, fnamepy, fnameyaml, savename, casename=None, 
        additional_vars = ['isupimpap', 'nusp_imp', 'nurlx', 'numvarbwpad']):
        """ Writes a YAML input file based on arbitrary Python script

        Use: 
            - open a new Python session (!!! CRITICAL !!!)
            - import and create an uninitialized Case object storing default values
                - from uetools import Case;c = Case(store_defaults=True)
            - execute the external script to restore the solution
            - call write_py2yaml

        Workings:
            - Tracks changes to all input values
            - Cross-references any changes to the Python input files, 
              keeping variables referenced
            - 

        - Remember to include any and all files used to set up the case (even
          those called internally) to get the complete input deck
        - If variables are missing, please reach out to the developers to
          add them to future releases 

        """
        from re import sub
        from yaml import dump

        if isinstance(fnamepy, str):
            files = [fnamepy]
        elif isinstance(fnamepy, list):
            files = fnamepy
        else:
            raise Exception("Unknown Python file specifier {}".format(str(type(fnamepy))))
        if casename is None:
            casename = savename.split('.')[0]

        vetted_changes = {}
        # Ensure all variables are populated
        self.populate()
        # Get a list of all variables changed during reading input and restoring
        self.record_changes()
        for file in files:
            with open(file) as f:
                # Dump all non-comment data to string
                filedump = sub(r'(?m)^ *#.*\n?', '', f.read())
            for var, value in self.varinput['setup']['detected'].items():
                if ".".join([self.getpackage(var), var]) in filedump:
                    vetted_changes[var] = value
        
        inputdeck = {
            'casename': casename,
            'savefile': savename,
            'restart': 1,
            'species': {},
            'grid': {'gengrid': 0, 'isgriduehdf5': 1, 'GridFileName': savename},
            'diffusivities': {},
            'variables': {},
            'equations': {},
        }
        
        self.setue('gengrid', 0)
        self.setue('isgriduehdf5', 1)
        self.setue('restart', 1)
        self.setue('GridFileName', savename)

        # Set up the species array
        for var in ['nhsp', 'ngsp', 'isimpon', 'nzsp']:
            if var in vetted_changes:
                self.write_changes(var, vetted_changes[var], inputdeck['species'])
                # Remove to avoid duplication
                del( vetted_changes[var] )
        # Set up the grid deck
        for var in ['mhdgeo', 'geometry', 'isnonog']: 
            if var in vetted_changes:
                self.write_changes(var, vetted_changes[var], inputdeck['grid'])
                del( vetted_changes[var] )
            else:
                self.write_changes(var, self.getue(var), inputdeck['grid'])
        # Set up the diffusivity deck
        for var in ['isbohmcalc',  'kye',  'kyi',  'difni',  'fcdif',  'travis']:
            if var in vetted_changes:
                self.write_changes(var, vetted_changes[var], inputdeck['diffusivities'])
                del( vetted_changes[var])
        # Remove arrays restored from save-file
        for var in [ 'kye_use',  'kyi_use',  'dif_use',  'tray_use',  'kyev',  'kyiv',  'difniv',  'travisv']:
            if var in vetted_changes:
                del( vetted_changes[var] )
        # Set up the equations deck
        for var in ['isteon', 'istion', 'isnion', 'isupon', 'isphion', 
                    'isphiofft', 'isngon', 'isupgon', 'istgon']:
            if var in vetted_changes:
                self.write_changes(var, vetted_changes[var], inputdeck['equations'])
                del( vetted_changes[var])
        # Set up any remaining variables
        for var, value in vetted_changes.items():
            self.write_changes(var, value, inputdeck['variables'])
        # Patch for missing input flags for parameters: always include to be safe
        for var in additional_vars:
            self.write_changes(var, self.getue(var), inputdeck['variables'])

        self.strip_numpy(inputdeck)
        # TODO: write out YAML file
        filedump = dump(inputdeck, sort_keys=False, default_style=None, default_flow_style=False)
        # Fix styling of lists
        filedump = filedump.replace("'[","[").replace("]'", "]")
        with open(fnameyaml, 'w') as f:
            f.write(filedump)
        self.varinput['setup'] = inputdeck
        self.reload()
        self.save(savename)


    def strip_numpy(self, struct, parent=[]):
        from numpy import int64, bytes_, float64
        from copy import deepcopy
        for key, substruct in struct.items():
            keys = deepcopy(parent)
            keys.append(key)
            if isinstance(substruct, dict):
                self.strip_numpy(substruct, keys)
            elif isinstance(substruct, list):
                lst = []
                for element in substruct:
                    if isinstance(element, (int, int64)):
                        lst.append(int(element))
                    elif isinstance(element, (float)):
                        lst.append(float(element))
                struct[key] = str(lst)
            elif isinstance(substruct, (int, int64)):
                struct[key] = int(substruct)
            elif isinstance(substruct, str):
                continue
            elif isinstance(substruct, float):
                struct[key] = float(substruct)
            elif isinstance(substruct, bytes_):
                struct[key] = substruct.decode('UTF-8')
            else:
                print("Through: ", keys, key, type(substruct), isinstance(substruct, dict), substruct)

        
    def set1Dbuff(self, orig, value): 
        from numpy import equal, all, int64
       # Initialize dict
        buff = {}
        # Get boolean array of changes
        elements = equal(orig, value)
        # Create dict of indices and new values
        for i in range(len(value)):
            if not elements[i]:
                buff[i] = value[i]
        # All values have changed
        if len(buff) == len(value):
            # All set to same value
            if all(value == value[0]):
                buff = value[0]
            # All set to different values
            else:
                vallist = []
                for _, val in buff.items():
                    vallist.append(val)
                buff = vallist
        # Multiple entries changed, try to compress
        if not isinstance(buff, (int, float, int64, list)):
            if len(buff)>1:
                keys = list(buff.keys())
                if keys == list(range(keys[0], keys[-1]+1)):
                    vallist = []
                    for _, val in buff.items():
                        vallist.append(val)
                    buff = {keys[0]: vallist}
        return buff

    def set2Dbuff(self, value, orig):
        from copy import deepcopy
        from numpy import all, equal
        # Grab an entry
        comp = deepcopy(value)
        for i in range(len(value.shape)):
            comp = comp[0]
        # If all entries are the same, just set the wole array
        if all(value == comp):
            return comp
        # Shape assumes target setup
        elif (len(value.shape) == 2) and (value.shape[-1] == 2):
            # Symmetrical settings
            if sum(abs(value[:,1]-value[:,0]))/sum(sum(abs(value))) < 1e-6:
                return self.set1Dbuff(orig[:,0], value[:,0])
            else:
                return False
        # TODO: implement 2D/dynamic handling
        else:
            return False

    def write_changes(self, var, value, location):
        from numpy import ndarray, bytes_
        from copy import deepcopy
        # Set simple int/float by value: easy
        if isinstance(value, (int, float)):
            location[var] = value
            return
        # Array wit changes: complicated
        elif isinstance(value, ndarray):
            # Check whether the array has been dynamically changed/allocated 
            if (len(value) == 1) and isinstance(value[0], bytes_):
                buff = value[0].decode('UTF-8').strip()
            elif len(self.defaults[var]) != len(value):
                buff = self.set2Dbuff(value, self.defaults[var])
            # Check whether the input is multi-dimensional (2D assumed max)
            elif len(value.shape) > 1:
                buff = self.set2Dbuff(value, self.defaults[var])
            # 1D array of settings detected
            else:
                buff = self.set1Dbuff(self.defaults[var], value)
            if buff is not False:
                location[var] = buff
            else:
                print("Could not write automatic input for ", var)
            return
        else:
            print("Unknown type '{}': could not write automatic input for '{}'".format(type(value), var))



