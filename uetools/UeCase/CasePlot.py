# Plotting routines for case objects
from uetools.UePlot import Plot

class Caseplot(Plot):
    def __init__(self):
        super(Caseplot, self).__init__()
        # Calculate distances along targets
        self.otdistance = self.get('yyrb')
        self.itdistance = self.get('yylb')
        return


    # TODO: implement profile plots


    # TODO: fix this function 
    def radialdistance(self, ind):
        from numpy import cumsum
        # Coordinates of edge of grid
        dr0 = self.get('rm', 2)[ind, 1]
        dz0 = self.get('zm', 2)[ind, 1]
        # Halfway point of cell centers
        drplatecenter = 0.5*(self.get('rm', 2)[ind] + self.get('rm', 4)[ind])
        dzplatecenter = 0.5*(self.get('zm', 2)[ind] + self.get('zm', 4)[ind])
        # Calculate strike point coordinates
        iysptrx = self.get('iysptrx')
        rsp = self.get('rm', 4)[ind, iysptrx] - drplatecenter[iysptrx]
        zsp = self.get('zm', 4)[ind, iysptrx] - dzplatecenter[iysptrx]
        # Normalize to grid edge
        drplatecenter -= dr0
        dzplatecenter -= dz0
        # Calculate cumulative distance along target
        dist = cumsum((drplatecenter**2 + dzplatecenter**2)**0.5)
        # Calculate distance along plate to strike point
        dsep = dist[iysptrx] + (rsp**2 + zsp**2)**0.5
        return dist-dsep

    def it(self, variable, marksep=True, staggered=False, **kwargs):
        # Exchange YYC for working radialdistance
        fig = self.plotprofile(self.itdistance[1:-1], 
            variable[0**staggered, 1:-1], **kwargs)
        # Add Sep location if requested
        if marksep is True:
            fig.get_axes()[0].axvline(0, color='g', linewidth=1)
        return fig

    def ot(self, variable, marksep=True, **kwargs):
        # Exchange YYC for working radialdistance
        fig = self.plotprofile(self.otdistance[1:-1], 
            variable[-2, 1:-1], **kwargs)
        # Add Sep location if requested
        if marksep is True:
            fig.get_axes()[0].axvline(0, color='g', linewidth=1)
        return fig

    def omp(self):
        return

    def imp(self):
        return

    def row(self):
        return

    def neOT(self, s=None, **kwargs):
        return self.ot(self.get('ne', s), **kwargs)

    def teOT(self, s=None, **kwargs):
        return self.ot(self.get('te', s)/self.get('ev'), **kwargs)

    def niOT(self, s=None, **kwargs):
        return self.ot(self.get('ni', s), **kwargs)

    def tiOT(self, s=None, **kwargs):
        return self.ot(self.get('ti', s)/self.get('ev'), **kwargs)

    # Expand the 2D plot list
    def grid(self, **kwargs):
        return self.plotmesh(**kwargs)

    def heatmap(self, var, s=None, **kwargs):
        return self.plotmesh(self.get(var, s), **kwargs)

    def ne2D(self, **kwargs):
        return self.variablemesh(self.get('ne'), **kwargs)

    def te2D(self, s=None, **kwargs):
        return self.variablemesh(self.get('te', s)/self.get('ev'), **kwargs)

    def ni2D(self, s, **kwargs):
        return self.variablemesh(self.get('ni', s), **kwargs)

    def ti2D(self, **kwargs):
        return self.variablemesh(self.get('te')/self.get('ev'), **kwargs)

    def ng2D(self, s, **kwargs):
            return self.variablemesh(self.get('ng', s), **kwargs)

    def tg2D(self, s, **kwargs):
        return self.plotmesh(self.get('tg', s)/self.get('ev'), **kwargs)

    def CIII_emission_2D(self, fname, **kwargs):
        self.emission_CIII(fname)
        return self.plotmesh(self.CIII_emission, **kwargs)

    def masked_CIII_2D(self, fname, z, maskvalues, **kwargs):
        self.emission_CIII(fname)
        mask = self.CIII_emission[1:-1,1:-1].reshape(self.nx*self.ny)
        mask = [1*( (x<maskvalues[0]) or (x>maskvalues[1])) for x in mask]
        kwargs['mask'] = mask
        return self.plotmesh(z, **kwargs)

    def CIIImasked_flow(self, fname, maskvalues, **kwargs):
        z = (self.get('upi')**2 + self.get('vy')**2)**0.5
        z = self.get('upi')
#        return self.masked_CIII_2D(fname, self.get('ne'), maskvalues, **kwargs)
        return self.masked_CIII_2D(fname, z[:,:,4], maskvalues, **kwargs)


    def variablemesh(self, z=None, **kwargs):
        from matplotlib.pyplot import ion, ioff, subplots
        from matplotlib.widgets import Slider, RangeSlider
        from matplotlib import is_interactive
        ioff()
        
        f, ax = subplots(figsize=(5,8))
        try:
            kwargs['zrange']
            origrange = kwargs['zrange']
        except:
            kwargs['zrange'] = (z.min(), z.max())
            origrange = kwargs['zrange']
        cbar, verts = self.plotmesh(z, ax=ax, watermark=False, 
            interactive=True, **kwargs)
        zrange_position = f.add_axes([0.95, 0.1, 0.04, 0.8])
        zrange_slider = RangeSlider(zrange_position, '', z.min(), z.max(),
            valinit=(origrange), orientation='vertical',
            valstep = round((z.max()-z.min())/100)) 

        def update(val):
            from numpy import floor, ceil, linspace
            zrange = zrange_slider.val
            verts.set_clim(zrange)

        zrange_slider.on_changed(update)
        f.show() 


        ion()

