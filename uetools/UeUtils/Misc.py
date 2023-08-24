

class Misc():


    def about_setup(self):
        """ Function outputting data about setup """

        nisp = self.get('nisp')
        ngsp = self.get('ngsp')
        nhsp = self.get('nhsp')
        nhgsp = self.get('nhgsp')
    


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
