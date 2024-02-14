
class Convert:
    def write_py(self, fname):
        """ Writes a standalone Python input file for the Case
        """
        try:
            from Forthon import package, packageobject
        except:
            pass

        def recursive_lineread(dictobj, lines=None, fails=None):
            """ Recusively parses setup variables to input lines """
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
            f.write('casename = "{}"\n'.format(self.casename))
            # Set paths to rate files as defined in configuration file
            f.write('aph.aphdir = "{}"\n'.format(self.get('aphdir')[0].decode('UTF-8').strip()))
            f.write('api.apidir = "{}"\n'.format(self.get('apidir')[0].decode('UTF-8').strip()))
            # WRITE LINES
            for line in lines:
                f.write(line)
                f.write('\n')
            # Allocate to make space for restore variables
            f.write('\nbbb.allocate()\n'.format(self.savefile))
            # Restore solution: use inplace version to ensure compatibility
            # with both native UEDGE saves and UETOOLS saves
            f.write('\n# ==== RESTORE SOLUTION ====\n')
            f.write('with File("{}") as f:\n'.format(self.savefile))
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
            if self.diff_file is not None:
                file = self.diff_file
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


    def py2yaml(self, fnamepy, fnameyaml, blockseparator='===='):
        try:
            from Forthon import package
        except:
            pass
        from yaml import dump
        yamlinput = {}
        packages = package()
        block = ''
        others = []
        index = None
        with open(fnamepy, 'r') as f:
            for line in f:
                # Omit empty lines
                if len(line.split()) == 0:
                    continue
                elif 'bbb.allocate' in line:
                    continue
                elif 'casename' in line:
                    yamlinput['casename'] = line.split('=')[-1].split('#')[0].strip()
                # Line is a variable
                elif line.split('.')[0] in packages:
                    if 'USER-SPECIFIED'.lower() in block:
                        others.append(line)
                    else:
                        var, value = line.split('=')
                        value=value.strip()
                        var = var.split('.')[1]
                        # Variable setting indices
                        if '[' in var:
                            indices = var.split('[')[-1].split(']')[0]
                            var = var.split('[')[0]
                            # Setting range
                            if ':' in indices:
                                # Set from start of array
                                if indices[0] == ':':
                                    index = None
                                # Set to start from specific index
                                else:
                                    index = int(indices.split(':')[0])
                            # Only one entry is set
                            else:
                                index = int(indices)
                            if '[[' in value:
                                value=value.replace('[[','[').\
                                    replace(']]',']').replace(' ','').\
                                    replace('],[',', ')
                            if index is None:
                                yamlinput[block][var] = value
                            else:
                                yamlinput[block][var] = {index: value}
                        else:
                            # Set undefined parameters
                            if len(block) == 0:
                                yamlinput[var] = value
                            # Write non-indexed paramters to dict
                            else:
                                yamlinput[block][var] = value
                # Line contains a block
                elif blockseparator in line:    
                    # Set block header: check from beginning and end
                    block = line.split(blockseparator)[-2].strip().lower()
                    if block == '#':
                        block = line.split(blockseparator)[1].strip().lower()
                    # Set block title
                    if ('user-specified' not in block) and ('restore solution' not in block):
                        yamlinput[block] = {}
                # Omit comment lines
                elif line.strip()[0] == '#':
                    continue
                elif 'restore solution' in block:
                    if 'File'  in line:
                        yamlinput['savefile'] = line.split('"')[1]
                else:
                    # Gather any fallen-through lines here --> commands?
                    others.append(line)

        yamlinput['commands'] = []
        for line in others:
            yamlinput['commands'].append(line.replace('\n',''))
        with open(fnameyaml, 'w') as f:
            dump(yamlinput, f, sort_keys=False,default_style=None, default_flow_style=False)
        lines = []
        commands = False
        with open(fnameyaml, 'r') as f:
            for line in f:
                lines.append(line)
        with open(fnameyaml, 'w') as f:
            for line in lines:
                if 'e' in line:
                    try:
                        linevar, lineval = line.split(':')
                    except:
                        linevar, lineval = ' ', ' '
                        pass
                    if lineval.split('e')[0][-1] in [str(x) for x in range(10)]:
                        if '.' not in lineval:
                            lineval = lineval.replace('e','.e')
                            line = '{}:{}'.format(linevar, lineval)
                if 'commands' in line:
                    commands=True
                if commands is False:
                    f.write(line.replace("'",""))
                else: 
                    if 'commands' in line:
                        f.write(line)
                    else:
                        f.write("{}'".format(line.strip().replace("'","").replace("- ", "- '")))
                        f.write('\n')
        



