# Object for plotting 
from matplotlib.pyplot import ion
ion()

class Plot():
    def __init__(self, rm=None, zm=None):
        ''' Constructs patches objects
        rm - UEDGE R-node object
        zm - UEDGE Z-node object
        '''
        # Intialize empty
        if (rm is None) or (zm is None):
            try:
                rm = self.get('rm')
                zm = self.get('zm')
            except:
                return
            
        # TODO: figure out why createpolycollection bogs down Datbase?
        if self.database is not True:
            self.createvertices(rm, zm)
        return
 
    def createvertices(self, rm, zm):
        self.vertices = self.createpolycollection(rm, zm)
        if self.get('geometry')[0].strip().lower().decode('UTF-8') == \
            'uppersn':
            self.disp = 0
            if self.get('rmagx') + self.get('zmagx') == 0:
                self.disp = -(-zm).min()
            else:
                self.disp = 2*self.get('zmagx')
            self.uppersnvertices = self.createpolycollection(rm, \
                -zm + self.disp, setparams=False)
 

    def createpolycollection(self, rm, zm, margins=0.05, setparams=True):
        ''' Creates a poly collection and records boundaries
        ''' 
        from matplotlib.collections import PolyCollection
        from uedge import com

        if setparams is True:
            self.nx = rm.shape[0]-2
            self.ny = rm.shape[1]-2
            # TODO: where to pass/find xpoint indices?
            self.isepr = rm[1:self.get('ixpt1')[0]+2,self.get('iysptrx')+1, 1]
            self.isepz = zm[1:self.get('ixpt1')[0]+2,self.get('iysptrx')+1, 1]
            self.osepr = rm[self.get('ixpt2')[0]:-1,self.get('iysptrx')+1, 2]
            self.osepz = zm[self.get('ixpt2')[0]:-1,self.get('iysptrx')+1, 2]

        vertices = []
        # Loop through all cells, omitting guard cells
        for i in range(1,len(rm)-1):
            for j in range(1,len(rm[i])-1):
                vert = []
                for k in [1, 2, 4, 3]:
                    vert.append([rm[i, j, k], zm[i, j, k]]) 
                vertices.append(vert)
        return PolyCollection(vertices)

    def checkusn(self, array, flip=False):
        if flip is False:
            return array
        elif (self.get('geometry')[0].strip().lower().decode('UTF-8') == \
            'uppersn'):
            return -array + self.disp

    def plotprofile(self, x, y, ax=None, xlim=(None, None), ylim=(0, None),
        figsize=(6,5), xlabel=None, ylabel=None, title=None, logx=False, 
        logy=False, color='k', watermark=True, **kwargs):
        ''' Plots y as function of x '''
        from matplotlib.pyplot import figure, Axes, Figure

        if ax is None: # Create figure if no axis supplied
            f = figure(title, figsize=figsize)
            ax = f.add_subplot()
        elif ax is Figure:
            ax = ax.get_axes()[0]
        # Switch to identify requested plot type
        if logx and logy:
            plot = ax.loglog
        elif logx and not logy:
            plot = ax.semilogx
        elif logy and not logx:
            plot = ax.semilogy
        else:
            plot = ax.plot
    
        plot(x, y, color=color, **kwargs)

        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        if watermark is True:
            self.watermark(ax.get_figure())

        return ax.get_figure()

    def plotmesh(self, z=None, rm=None, zm=None, ax=None, linewidth=0.05,
        linecolor='k', aspect='equal', figsize=(5,7), cmap='magma', units='', 
        xlim=(None, None), ylim=(None, None), zrange=(None, None), 
        log=False, vessel=True, plates=True, lcfs=True, title=None, 
        grid=False, flip=False, watermark=True, mask=None, colorbar=True,
        interactive = False):
        ''' General plotting function
        z - values, if any. If None, plots grid
        rm, zm - radial and horizontal nodes
        '''
        from matplotlib.pyplot import figure, Figure
        from matplotlib.colors import LogNorm
        from matplotlib.collections import PolyCollection
        from copy import deepcopy
        from numpy import array
        from uedge import com, bbb, grd        

        try:
            self.vertices
        except:
            self.createvertices(self.get('rm'), self.get('zm'))

        if ax is None:
            f = figure(title, figsize=figsize)
            ax = f.add_subplot()
        if ax is Figure:
            ax = ax.get_axes()[0]
        if (rm is None) or (zm is None):
            # Use stored PolyCollection
            if (self.get('geometry')[0].strip().lower().decode('UTF-8') == \
                'uppersn') and (flip is True): 
                vertices = deepcopy(self.uppersnvertices)
            else:
                vertices = deepcopy(self.vertices)
        else: # Create collection from data
            vertices = self.createpolycollection(rm, zm)
        if grid is False:
            vertices.set_linewidths(1)
            vertices.set_edgecolors('face')
        else:
            vertices.set_edgecolors(linecolor)
            vertices.set_linewidths(linewidth)
        if z is None: # Plot grid
            vertices.set_facecolor((0, 0, 0, 0))
            vertices.set_edgecolors(linecolor)
            vertices.set_linewidths(linewidth)
        else:
            vertices.set_cmap(cmap)
            vertices.set_array(z[1:-1,1:-1].reshape(self.nx*self.ny))
            vertices.set_clim(*zrange)
            if log is True:
                vertices.set_norm(LogNorm())
                vertices.set_clim(*zrange)
        if mask is not None:
            vertices.set_alpha(mask)
            
        ax.add_collection(vertices)
        # TODO: devise scheme to look for variables in memory, from 
        # Forthon, from HDF5
        if lcfs is True:
            self.plotlcfs(ax, flip)
        if vessel is True:
            self.plotvessel(ax, flip)
        if plates is True:
            self.plotplates(ax, flip)
        ax.autoscale_view()
        ax.set_title(title)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel('R [m]')
        ax.set_ylabel('Z [m]')
        ax.set_aspect(aspect)
        if (z is not None) and (colorbar is True):
            cbar = ax.get_figure().colorbar(vertices, ax=ax)
            cbar.ax.set_ylabel(units, va='bottom')
    
        if watermark is True:
            self.watermark(ax.get_figure(), bottom=0.1, left=0.02, 
                right=0.95)

        if interactive is True:
            return cbar, vertices
        else:
            return ax.get_figure()
        
           

    def plotlcfs(self, ax, flip=False, color='grey', linewidth=0.5,
        **kwargs):
        """ Plots LCFS on ax """
        from uedge import com, bbb, grd        
        
        if self.get('geometry')[0].strip().lower().decode('UTF-8') == 'dnull':
            rm = self.get('rm')
            zm = self.get('zm')
            iysptrx1 = self.get('iysptrx1')
            iysptrx2 = self.get('iysptrx2')
            ixrb = self.get('ixrb')
            ixlb = self.get('ixlb')
            ixpt1 = self.get('ixpt1')
            ixpt2 = self.get('ixpt2')

            ax.plot(rm[:ixrb[0]+1, iysptrx1[0]+1, 1],
                zm[:ixrb[0]+1, iysptrx1[0]+1, 1], color=color, 
                linewidth=linewidth)
            ax.plot(rm[ixlb[1]:, iysptrx1[0]+1, 1],
                zm[ixlb[1]:, iysptrx1[0]+1, 1], color=color, 
                linewidth=linewidth)

            ax.plot(rm[:ixpt1[0]+1, iysptrx2[0]+1, 2],
                zm[:ixpt1[0]+1, iysptrx2[0]+1, 2], 
                color=color, linewidth=linewidth)
            ax.plot(rm[ixpt2[1]+1:, iysptrx2[0]+1, 1],
                zm[ixpt2[1]+1:, iysptrx2[0]+1, 1], 
                color=color, linewidth=linewidth)
        
            ax.plot(rm[ixpt1[0]+1:ixrb[0]+1, iysptrx2[0]+1, 1],
                zm[ixpt1[0]+1:ixrb[0]+1, iysptrx2[0]+1, 1], 
                color=color, linewidth=linewidth)
            ax.plot(rm[ixlb[1]:ixpt2[1]+1, iysptrx2[0]+1, 2],
                zm[ixlb[1]:ixpt2[1]+1, iysptrx2[0]+1, 2], 
                color=color, linewidth=linewidth)

        else:
            plotted = False
            try:
                ax.plot(com.rbdry, self.checkusn(com.zbdry, flip), 
                    color=color, linewidth=linewidth)
                ax.plot(self.isepr, self.checkusn(self.isepz, flip), 
                    color=color, linewidth=linewidth)
                ax.plot(self.osepr, self.checkusn(self.osepz, flip), 
                color=color, linewidth=linewidth)
                plotted = True
            except:
                pass
            try:
                if self.get('rbdry') != False:
                    ax.plot(self.get('rbdry'), self.checkusn(self.get(\
                        'zbdry'), flip), color=color, 
                        linewidth=linewidth)
                    ax.plot(self.isepr, self.checkusn(self.isepz, flip), 
                        color=color, linewidth=linewidth)
                    ax.plot(self.osepr, self.checkusn(self.osepz, flip), 
                        color=color, linewidth=linewidth)
                    plotted = True
            except:
                pass
            if plotted is False:
                ax.plot(self.get('rm')[com.ixpt1[0]:com.ixpt2[0]+1,\
                    com.iysptrx+1,2], self.checkusn(self.get('zm')[\
                    com.ixpt1[0]:com.ixpt2[0]+1,com.iysptrx+1,2], flip),
                    color=color, linewidth=linewidth)
                ax.plot(self.isepr, self.checkusn(self.isepz, flip), 
                    color=color, linewidth=linewidth)
                ax.plot(self.osepr, self.checkusn(self.osepz, flip), 
                    color=color, linewidth=linewidth)

    def plotvessel(self, ax, flip=False):
        """ Plots vessel on ax """
        from uedge import com, bbb, grd        
        try:
            ax.plot(com.xlim, self.checkusn(com.ylim, flip), 'k-', 
                linewidth=3)
            ax.plot(com.xlim, self.checkusn(com.ylim, flip), 
                'y-', linewidth=1)
        except:
            pass
        try:
            if self.get('xlim') != False:
                ax.plot(self.get('xlim'), self.checkusn(self.get('ylim'), 
                    flip), 'k-', linewidth=3)
                ax.plot(self.get('xlim'), self.checkusn(self.get('ylim'), 
                    flip), 'y-', linewidth=1)
        except:
            pass

    def plotplates(self, ax, flip=False):
        """ Plot plates on ax """    
        from uedge import com, bbb, grd        
        try:
            ax.plot(grd.rplate1, self.checkusn(grd.zplate1, flip), 'b-', 
                linewidth=1.5)
            ax.plot(grd.rplate2, self.checkusn(grd.zplate2, flip), 'r-', 
                linewidth=1.5)
        except:
            pass
        try:
            if self.get('rplate1') != False:
                ax.plot(self.get('rplate1'), self.checkusn(self.get(\
                    'zplate1'), flip), 'b-', linewidth=1.5)
                ax.plot(self.get('rplate2'), self.checkusn(self.get(\
                    'zplate2'), flip), 'r-', linewidth=1.5)
        except:
            pass
    
    def watermark(self, figure, bottom=0.15, top=0.95, left=0.09, right=0.98):
        """ Adds metadata to figure """
        from uedge import __version__
        from time import ctime
        label = '{}, case "{}"\n'.format(ctime(), self.casename)
        label += 'UEDGE {} v{}, UETOOLS v{}, user "{}", hostname "{}"\n'\
            .format(self.uedge_ver.replace('$','\$'), self.pyver, 
            self.uetoolsversion, self.user, self.hostname)
        try:
            label += 'cwd "{}"'.format(self.location)
        except:
            label += 'cwd "{}"'.format(self.casefname)
        figure.subplots_adjust(bottom=bottom, top=top, left=left, 
            right=right)
        figure.text(0.995, 0.005, label, fontsize=4, 
            horizontalalignment='right')
        
        return


    def plotmesh_masked(self, z, zmask, maskvalues, **kwargs):
        from copy import deepcopy

        try:
            rm = kwargs['rm']
        except:
            rm = None
        try:
            zm = kwargs['zm']
        except:
            zm = None
        try:
            grid = kwargs['grid']
        except:
            grid = False
        try:
            log = kwargs['log']
        except:
            log = False
        try:
            cmap = kwargs['cmap']
        except:
            cmap = 'magma'
        try:
            zrange = kwargs['zrange']
        except:
            zrange = (None, None)
        
        if (rm is None) or (zm is None):
            # Use stored PolyCollection
            if (self.get('geometry')[0].strip().lower().decode('UTF-8') == \
                'uppersn') and (flip is True): 
                vertices = deepcopy(self.uppersnvertices)
            else:
                vertices = deepcopy(self.vertices)
        else: # Create collection from data
            vertices = self.createpolycollection(rm, zm)
        if grid is False:
            vertices.set_linewidths(1)
            vertices.set_edgecolors('face')
        else:
            vertices.set_edgecolors(linecolor)
            vertices.set_linewidths(linewidth)
        vertices.set_cmap(cmap)
        vertices.set_array(z[1:-1,1:-1].reshape(self.nx*self.ny))
        mask = z[1:-1,1:-1].reshape(self.nx*self.ny)
        vertices.set_alpha( [1*( (x<maskvalues[0]) or (x>maskvalues[1])) for x in mask])
        vertices.set_clim(*zrange)
        if log is True:
            vertices.set_norm(LogNorm())
            vertices.set_clim(*zrange)
        self.plotmesh(vertices, **kwargs)
        return mask
        
    def streamline(self, pol, rad, resolution=(500j,800j), linewidth=0.5,  
        broken_streamlines=False, color='k', maxlength=0.4, mask=True, 
        **kwargs):
        from numpy import zeros, sum, transpose, mgrid, nan, array
        from scipy.interpolate import griddata, bisplrep
        from matplotlib.patches import Polygon
        rm = self.get('rm')
        zm = self.get('zm')
        nx = self.get('nx')
        ny = self.get('ny')
        nodes = zeros((nx+2, ny+2, 5, 2))
        nodes[:, :, :, 0] = rm
        nodes[:, :, :, 1] = zm
        nodes = transpose(nodes, (2, 3, 0, 1))
        # Create polygons for masking
        outerx = []
        outerx = outerx + list(rm[::-1][-self.get('ixpt1')[0]:,0,2])
        outerx = outerx + list(rm[0,:,1])
        outerx = outerx + list(rm[:,-1,3])
        outerx = outerx + list(rm[:,::-1][-1,:,4])
        outerx = outerx + list(rm[::-1][:nx-self.get('ixpt2')[0],0,1])
        outery = []
        outery = outery + list(zm[::-1][-self.get('ixpt1')[0]:,0,2])
        outery = outery + list(zm[0,:,1])
        outery = outery + list(zm[:,-1,3])
        outery = outery + list(zm[:,::-1][-1,:,4])
        outery = outery + list(zm[::-1][:nx-self.get('ixpt2')[0],0,1])

        innerx = rm[self.get('ixpt1')[0]+1:self.get('ixpt2')[0]+1,0,1]
        innery = zm[self.get('ixpt1')[0]+1:self.get('ixpt2')[0]+1,0,1]

        outer = Polygon(array([outerx, outery]).transpose(), closed=True,
            facecolor='white', edgecolor='none')
        inner = Polygon(array([innerx, innery]).transpose(), closed=True,
            facecolor='white', edgecolor='none')

        # Find midpoints of y-faces
        symid = zeros((2, 2, nx+2, ny+2))
        symid[0] = (nodes[2] + nodes[1])/2 # Lower face center
        symid[1] = (nodes[4] + nodes[3])/2 # Upper face center
        # Find midpoints of x-faces
        sxmid = zeros((2, 2,  nx+2, ny+2))
        sxmid[0] = (nodes[3] + nodes[1])/2 # Left face center
        sxmid[1] = (nodes[4] + nodes[2])/2 # Right face center
        # Get normal vectors in y-direction of each cell
        ynormal = zeros((2, nx+2, ny+2))
        ynormal = symid[1] - symid[0]
        # Get normal vectors in x-direction of each cell
        xnormal = zeros((2, nx+2, ny+2))
        xnormal = sxmid[1] - sxmid[0]
        # Get radial unit vectors
        ynormalhat = zeros((2, nx+2, ny+2))
        for i in range(2):
            ynormalhat[i] = ynormal[i]/sum(ynormal**2, axis=0)**0.5
        # Get poloidal unit vectors
        xnormalhat = zeros((2, nx+2, ny+2))
        for i in range(2):
            xnormalhat[i] = xnormal[i]/sum(xnormal**2, axis=0)**0.5

        x = pol * xnormalhat[0] + rad * ynormalhat[0]
        y = pol * xnormalhat[1] + rad * ynormalhat[1]

        gx, gy = mgrid[rm.min():rm.max():resolution[0], 
            zm.min():zm.max():resolution[1]]
    
        # Previous implementation, which (incorrectly) assumed values to 
        # be given for cell-centers
#        xinterp = griddata( (rm[1:-1,1:-1,0].ravel(), zm[1:-1,1:-1,0].ravel()), 
#            x[1:-1,1:-1].ravel(), (gx, gy)) 
#        yinterp = griddata( (rm[1:-1,1:-1,0].ravel(), zm[1:-1,1:-1,0].ravel()), 
#            y[1:-1,1:-1].ravel(), (gx, gy)) 
        xinterp = griddata( (sxmid[1,0,1:-1,1:-1].ravel(), sxmid[1, 1, 1:-1, 1:-1].ravel()),
            x[1:-1,1:-1].ravel(), (gx, gy)) 
        yinterp = griddata( (symid[1,0,1:-1,1:-1].ravel(), symid[1, 1, 1:-1, 1:-1].ravel()),
            y[1:-1,1:-1].ravel(), (gx, gy)) 

        if mask is True:
            for i in range(gx.shape[0]):
                for j in range(gx.shape[1]):
                    p = (gx[i,j], gy[i,j])
                    if (inner.contains_point(p)) or (not outer.contains_point(p)):
                        xinterp[i,j] = nan
                        yinterp[i,j] = nan

        f=self.plotmesh()
        if linewidth == 'magnitude':
            linewidth = (xinterp**2 + yinterp**2)**0.5
            linewidth /= liewidth.max()
        
        f.get_axes()[0].streamplot(gx.transpose(), gy.transpose(), 
            xinterp.transpose(), yinterp.transpose(), linewidth=linewidth,
            broken_streamlines=broken_streamlines, color=color, 
            maxlength=maxlength, **kwargs)


        return f

