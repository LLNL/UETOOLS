

class Misc():
    def __init__(self, case):
        self.get = case.get
        self.setue = case.setue
 
    def calc_lcon(self):
        """ Stores connection length in Case.lcon """
        from numpy import cumsum
        try:
            #TODO: How to extend to dnulls????
            # Calc and store connection lengths
            lcon = cumsum(1/(self.get('rr')*self.get('gx')+1e-20), axis=0)
            # No connlen for core cells
            lcon[self.get('ixpt1')[0]+1:self.get('ixpt2')[0]+1,
                    :self.get('iysptrx')+1] = 0
            lcon[self.get('ixpt2')[0]+1:, :self.get('iysptrx')+1] = \
                cumsum(1/(self.get('rr')*self.get('gx')+1e-20)[\
                self.get('ixpt2')[0]+1:,:self.get('iysptrx')[0]+1], axis=0
            )
        except:
            pass
        return lcon


        
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
        #TODO: How to extend to dnulls????
        from numpy import zeros_like
        phi = self.get('phi')
        phis = zeros_like(phi)
        gx = self.get('gx')
        ixp1 = self.get('ixp1')
        ex = self.get('ex')


        phis[1,:] = self.get('kappal')[:,0]*self.get('te')[1]/self.get('qe')
        if self.get('isudsym') ==0:
            for iy in range(1, self.get('ny')+1):
                for ix in range(1,self.get('nx')+1):
                    ix1 = ixp1[ix, iy]
                    dxf = 0.5*(gx[ix, iy] + gx[ix1, iy])/(gx[ix, iy]*gx[ix1, iy])
                    phis[ix1, iy] = phis[ix, iy] - ex[ix, iy] * dxf
                
            for ix in range(self.get('ixpt1')[0]+1, self.get('ixpt2')[0]+1):
                for iy in range(self.get('iysptrx')+1):
                    phis[ix, iy] = phis[ix, self.get('iysptrx')+1]
                phis[ix, self.get('ny')+1] = phis[ix, self.get('ny')]
        else:
            nxc = self.get('nxc')
            for iy in range(1, self.get('ny')+1):
                for ix in range(0, nxc):
                    ix1= ixp1[ix, iy]
                    dxf = 0.5*(gx[ix, iy] + gx[ix1, iy])/(gx[ix, iy]*gx[ix1, iy])
                    phis[ix1, iy] = phis[ix, iy] - ex[ix, iy] * dxf
                # NOTE: no overlap for the two subdomains in poloidal direction
                for ix in range(self.get('nx'), nxc, -1):
                    ix1= ixm1[ix, iy]
                    ix2= ixp1[ix, iy]
                    dxf = 0.5*(gx[ix, iy] + gx[ix1, iy])/(gx[ix, iy]*gx[ix1, iy])
                    phis[ix, iy] = phis[ix2, iy] - ex[ix, iy] * dxf
                    
            

        self.setue('phis', phis)
        self.setue('phi', phis)

    def squareinterp(self, data, r=None, z=None, method='linear', 
        resolution=(500j, 800j), mask=False, fill=float('NaN')):
        from scipy.interpolate import griddata
        from numpy import floor, ceil, mgrid, array, nan, concatenate
        from shapely.geometry import Point
        from shapely.geometry.polygon import Polygon
        # Get R,Z coords
        rm = self.get("rm")
        zm = self.get("zm")
        if self.get("geometry")[0].strip().lower().decode("UTF-8") == "     uppersn":        
            zm = self.disp - zm
        if mask:
            # Mask out external points
            x = []
            y = []
            for j in range(rm.shape[1]):
                x.append(rm[1, j, 1])
                y.append(zm[1, j, 1])
            for i in range(rm.shape[0]):
                x.append(rm[i, -1, 1])
                y.append(zm[i, -1, 1])
            for j in range(rm.shape[1]-1, -1, -1):
                x.append(rm[-1, j, 1])
                y.append(zm[-1, j, 1])
            for i in range(rm.shape[0]-1, self.get('ixpt2')[0], -1):
                x.append(rm[i, 0, 4])
                y.append(zm[i, 0, 4])
            for i in range(self.get(['ixpt1'])[0], -1, -1):
                x.append(rm[i, 0, 4])
                y.append(zm[i, 0, 4])
            if self.get("geometry")[0].strip().lower().decode("UTF-8") == "     uppersn":        
                y = [self.disp-a for a in y]
            outline = Polygon(concatenate([[array(x), array(y)]], axis=0).T)
            

        # Automatically set sqare grid boundaries
        if r is None:
            r = (0.01*floor(rm.min()*100), 0.01*ceil(rm.max()*100))
        if z is None:
            z = (0.01*floor(zm.min()*100), 0.01*ceil(zm.max()*100))
        # Get cell centers and unravel data
        gx, gy = mgrid[r[0]:r[1]:resolution[0], z[0]:z[1]:resolution[1]]
                
        interp = griddata( (rm[:,:,0].ravel(), zm[:,:,0].ravel()), data.ravel(), (gx, gy), method=method, fill_value=fill)
        if mask:
            for i in range(gx.shape[0]):
                for j in range(gx.shape[1]):
                    if not outline.contains(Point(gx[i,j], gy[i,j])):
                        interp[i,j] = nan

        # Perform interpolation
        return gx, gy, interp



#     from uetools.UePostproc.ADASclass import ADASSpecies
#     self.rates = ADASSpecies(path, species, ratetype, **kwargs)


    
    
