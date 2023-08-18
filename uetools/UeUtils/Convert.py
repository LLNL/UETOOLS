
class Convert:
    def write_py(self,fname):
        from Forthon import package, packageobject

        def recursive_lineread(dictobj, lines=None, fails=None):
            if lines is None:
                lines = []
            if fails is None:
                fails = {}
            for key, value in dictobj.items():
                if key in [
                            'casename', 
                            'savefile',
                            'commands', 
                            'userdifffname',
                            'radialdifffname',
                ]:
                    fails[key] = value
                    continue
                elif self.getpackage(key) is not None:
                    if not isinstance(value, dict):
                        if isinstance(value, list):
                            if len(self.getue(key).shape)>1:
                                valbuff = []
                                for v in value:
                                    valbuff.append([v])
                                value=valbuff
                            lines.append('{}.{}[:{}]={}'.format(
                                            self.getpackage(key), key, 
                                            len(value), value
                                        ))
                        elif value is False:
                            continue
                        elif isinstance(value, str):
                            lines.append('{}.{}="{}"'.format(
                                            self.getpackage(key), key, value
                                        ))
                        else:
                            lines.append('{}.{}={}'.format(
                                            self.getpackage(key), key, value
                                        ))
                    else:
                        lines.append('{}.{}[{}]={}'.format(
                                            self.getpackage(key), key, 
                                            '{}', '{}'
                                        ))
                elif isinstance(key, str):
                    lines.append('\n# ==== {} ===='.format(key.upper()))
                elif isinstance(key, int):
                    if isinstance(value, list):
                        index = '{}:{}'.format(key, key+len(value)) 
                    else:
                        index = key
                    if isinstance(value, str):
                        value = '"{}"'.format(value)
                    lines[-1] = lines[-1].format(index, value)
                else:
                    print('ERROR! Bottomed out after line: ', lines[-1])
                if isinstance(value, dict):
                    lines, fails = recursive_lineread(value, lines, fails)
                
                    
            return lines, fails

        lines, fails = recursive_lineread(self.varinput['setup'])
        
        try:
            fails.pop('casename')
        except:
            pass
        try:
            fails.pop('savefile')
        except:
            pass

        with open(fname, 'w') as f:
            f.write('from uedge import *\n')
            f.write('from h5py import File\n')
            f.write('from uedge.hdf5 import hdf5_restore\n\n')
            f.write('casename = "{}"\n'.format(self.casename))
            f.write('aph.aphdir = "{}"\n'.format(self.get('aphdir')[0].decode('UTF-8').strip()))
            f.write('api.apidir = "{}"\n'.format(self.get('apidir')[0].decode('UTF-8').strip()))
            # WRITE LINES
            for line in lines:
                f.write(line)
                f.write('\n')
            
            # TODO: allocate after every block?
            f.write('\nbbb.allocate()\n'.format(self.savefname))
            f.write('\n# ==== RESTORE SOLUTION ====\n')
            f.write('with File("{}") as f:\n'.format(self.savefname))
            f.write('    try:\n')
            f.write('        savegroup = f["restore/bbb"]\n')
            f.write('    except:\n')
            f.write('        savegroup = f["bbb"]\n')
            f.write('    for var in ["ngs", "nis", "phis", "tes", '
                        '"tgs", "tis", "ups"]:\n')
            f.write('        setattr(bbb, var, savegroup[var][()])\n')

            if 'commands' in fails:
                f.write('\n# ==== USER-SPECIFIED COMMANDS ====\n')
                for command in fails.pop('commands'):
                    if 'self.populate' in command:
                        f.write('# Populate UEDGE arrays\n')
                        f.write('bbb.issfon=0\nbbb.ftol=1e20\nbbb.exmain()\nbbb.issfon=1\nbbb.ftol=1e-8\n')
                    else:
                        f.write(command)
                        f.write('\n')
            
            if 'userdifffname' in fails:
                file = fails.pop('userdifffname')
                f.write('\n# ==== SET USER-SPECIFIED DIFFUSIVITIES ====\n')
                f.write('with File("{}", "r") as f:\n'.format(file))
                f.write('    bbb.dif_use=f["diffusivities/bbb/dif_use"][()]\n')
                f.write('    bbb.kye_use=f["diffusivities/bbb/kye_use"][()]\n')
                f.write('    bbb.kyi_use=f["diffusivities/bbb/kyi_use"][()]\n')
                f.write('    bbb.tray_use=f["diffusivities/bbb/tray_use"][()]\n')

            if 'radialdifffname' in fails:
                file = fails.pop('radialdifffname')
                f.write('\n# ==== SET USER-SPECIFIED DIFFUSIVITIES ====\n')
                f.write('with File("{}", "r") as f:\n'.format(file))
                f.write('    bbb.difniv=f["diffusivities/bbb/difniv"][()]\n')
                f.write('    bbb.kyev=f["diffusivities/bbb/kyev"][()]\n')
                f.write('    bbb.kyiv=f["diffusivities/bbb/kyiv"][()]\n')
                f.write('    bbb.travisv=f["diffusivities/bbb/travisv"][()]\n')


        if len(fails)>0:
            print('WARNING! Some variables could not be parsed into input!')
            print('These variables/groups were:')
            for key, value in fails.items():
                print('    -{}'.format(key))




        '''

        for header, variables in self.varinput['setup'].items():
            if header in ['casename', 'savefile']:
                continue # Write these from self
            elif self.getpackage is None:
                lines.append('# ==== {} ===='.format(header.upper()))
                for var, values in variables.items():
                    if isinstance(values, dict):
                        fails.append(values)
#                    else:
#                        print(var)
            else:
                print(header, variables)
                

        return lines
            
                

        with open(oldinput) as f:
            for line in f:
                pline = line.replace(':','=').replace('\n','').split('#')[0]
                var = pline.split('=')[0]
                print(var)
                if var in list(key.keys()):
                    lines.append('{}.{}'.format(key[var], pline))
                elif var in [str(x) for x in range(100)]:
                        val = pline.split('=')[1]
                        lines[-1] = lines[-1].replace('=','[{}]={}'.format(var,val))
                        print(lines[-1])
                else:
                    lines.append(line)
                    fails.append(var)

        with open(newinput,'w') as f:
            for l in lines:
                f.write(l)
                f.write('\n')
        return fails, lines
        '''
