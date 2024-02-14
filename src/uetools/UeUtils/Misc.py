

class Misc():
    def __init__(self, *args, **kwargs):
        """ Stores connection length in Case.lcon """
        from numpy import cumsum
        try:
            # Calc and store connection lengths
            self.lcon = cumsum(1/(self.get('rr')*self.get('gx')+1e-20), axis=0)
            # No connlen for core cells
            self.lcon[self.ixpt1+1:self.ixpt2+1,:self.iysptrx+1] = 0
            self.lcon[self.ixpt2+1:, :self.iysptrx+1] = \
                cumsum(1/(self.get('rr')*self.get('gx')+1e-20)[\
                self.ixpt2+1:,:self.iysptrx+1], axis=0
            )
        except:
            pass
        super().__init__(*args, **kwargs)
        return
 

    def smooth_curve(self, x, y, s=1, **kwargs):
        """ Returns smoothed x and y data """
        from scipy.interpolate import splrep, BSpline
        return x, BSpline(*splrep(x, y, s=s))(x)
        
        
    def get_bottomkeys(self, vardict, keys=None):
        """ Returns a list of keys at the bottom of nested dict """
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

    def potent_1dsol(self):
        from numpy import zeros_like
        phi = self.get('phi')
        phis = zeros_like(phi)
        gx = self.get('gx')
        ixp1 = self.get('ixp1')
        ex = self.get('ex')


        phis[1,:] = self.get('kappal')[:,0]*self.get('te')[1]/self.get('qe')
        if self.get('isudsym') ==0:
            for iy in range(1, self.ny+1):
                for ix in range(1,self.nx+1):
                    ix1 = ixp1[ix, iy]
                    dxf = 0.5*(gx[ix, iy] + gx[ix1, iy])/(gx[ix, iy]*gx[ix1, iy])
                    phis[ix1, iy] = phis[ix, iy] - ex[ix, iy] * dxf
                
            for ix in range(self.ixpt1+1, self.ixpt2+1):
                for iy in range(self.iysptrx+1):
                    phis[ix, iy] = phis[ix, self.iysptrx+1]
                phis[ix, self.ny+1] = phis[ix, self.ny]
        else:
            nxc = self.get('nxc')
            for iy in range(1, self.ny+1):
                for ix in range(0, nxc):
                    ix1= ixp1[ix, iy]
                    dxf = 0.5*(gx[ix, iy] + gx[ix1, iy])/(gx[ix, iy]*gx[ix1, iy])
                    phis[ix1, iy] = phis[ix, iy] - ex[ix, iy] * dxf
                # NOTE: no overlap for the two subdomains in poloidal direction
                for ix in range(self.nx, nxc, -1):
                    ix1= ixm1[ix, iy]
                    ix2= ixp1[ix, iy]
                    dxf = 0.5*(gx[ix, iy] + gx[ix1, iy])/(gx[ix, iy]*gx[ix1, iy])
                    phis[ix, iy] = phis[ix2, iy] - ex[ix, iy] * dxf
                    
            

        self.setue('phis', phis)
        self.setue('phi', phis)

    def get_detected_changes(self, savefile):
        from h5py import File
        from os.path import exists
        if not exists(savefile):
            raise ValueError(f"File {savefile} does not exist.")
        with File(savefile, 'r') as f:
            try: 
                detected = f['setup/detected']
            except:
                raise KeyError(f"File {savefile} is not a UETOOLS save")
            for key, item in detected.items():
                print(key)

    def hdf5tree(self, fname, group=None, depth=None, pre=''):
        from h5py import File, _hl
        from copy import deepcopy
        def typestr(var):
            for vt in ['int', 'float', 'bytes']:
                if vt in str(type(var)):
                    return '({})'.format(vt)
            return type(var)
        if depth is None:
            depth=1000
        def h5_tree(val, pre='', maxdepth=1000, depth=0):
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
                        print('{}└── {}'.format(pre, key))
                        h5_tree(val, pre+'    ', maxdepth, deepcopy(depth))
                    else:
                        try:
                            len(val)
                            if len(val.shape) == 1:
                                descr = '({})'.format(len(val))
                            else:
                                descr = val.shape
                            print('{}└── {} {}'.format(pre, key, descr))
                        except:
                            print('{}└── {} {}'.format(
                                                        pre, 
                                                        key, 
                                                        typestr(val[()])
                            ))
                else:
                    if type(val) == _hl.group.Group:
                        print('{}├── {}'.format(pre, key))
                        h5_tree(val, pre+'│   ', maxdepth, deepcopy(depth))
                    else:
                        try:
                            len(val)
                            if len(val.shape) == 1:
                                descr = '({})'.format(len(val))
                            else:
                                descr = val.shape
                            print('{}├── {} {}'.format(
                                                        pre, 
                                                        key, 
                                                        descr
                            ))
                        except:
                            print('{}├── {} {}'.format( pre, 
                                                        key, 
                                                        typestr(val[()])
                            ))

        if type(fname) == _hl.group.Group:
            h5_tree(fname, maxdepth=depth)
        else:
            with File(fname, 'r') as h5file:
                print(h5file)
                if group is not None:
                    print(group)
                    h5_tree(h5file[group], maxdepth=depth)
                else:
                    h5_tree(h5file, maxdepth=depth)


    def hdf5search(self, file, var):
        """ Searches file and returns value of var """
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
               


