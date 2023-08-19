
class Convert:
    def write_py(self,fname):
        """ Writes a standalone Python input file for the Case
        """
        from Forthon import package, packageobject

        def recursive_lineread(dictobj, lines=None, fails=None):
            """ Recusively parses setup variables to input lines """
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
                            lines.append('{}.{}={}'.format(
                                            self.getpackage(key), key, value
                                        ))
                    # Setting a non-starting index defined on separate line
                    # Create placeholder and fill later on
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
            f.write('\nbbb.allocate()\n'.format(self.savefname))
            # Restore solution: use inplace version to ensure compatibility
            # with both native UEDGE saves and UETOOLS saves
            f.write('\n# ==== RESTORE SOLUTION ====\n')
            f.write('with File("{}") as f:\n'.format(self.savefname))
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
            if 'userdifffname' in fails:
                file = fails.pop('userdifffname')
                f.write('\n# ==== SET USER-SPECIFIED DIFFUSIVITIES ====\n')
                f.write('with File("{}", "r") as f:\n'.format(file))
                f.write('    bbb.dif_use=f["diffusivities/bbb/dif_use"][()]\n')
                f.write('    bbb.kye_use=f["diffusivities/bbb/kye_use"][()]\n')
                f.write('    bbb.kyi_use=f["diffusivities/bbb/kyi_use"][()]\n')
                f.write('    bbb.tray_use=f["diffusivities/bbb/tray_use"][()]\n')
            # If a file is used to set the radial transport coefficients,
            # manually in the radial direction, uniform poloidally, read 
            # and restore them here
            if 'radialdifffname' in fails:
                file = fails.pop('radialdifffname')
                f.write('\n# ==== SET USER-SPECIFIED DIFFUSIVITIES ====\n')
                f.write('with File("{}", "r") as f:\n'.format(file))
                f.write('    bbb.difniv=f["diffusivities/bbb/difniv"][()]\n')
                f.write('    bbb.kyev=f["diffusivities/bbb/kyev"][()]\n')
                f.write('    bbb.kyiv=f["diffusivities/bbb/kyiv"][()]\n')
                f.write('    bbb.travisv=f["diffusivities/bbb/travisv"][()]\n')
        # Check whether there still are lines that could not be parsed 
        # into the input, and output them to the promt
        if len(fails)>0:
            print('WARNING! Some variables could not be parsed into input!')
            print('These variables/groups were:')
            for key, value in fails.items():
                print('    -{}'.format(key))


    # TODO: Implement the same routine for writing YAML files
    # TODO: Implement routine that converts python to yaml inpace/parses .py file to yaml


