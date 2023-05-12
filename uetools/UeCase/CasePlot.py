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

    def selector2D(self, var, interactive, **kwargs):
        if interactive is True:
            return self.variablemesh(var, **kwargs)
        else:
            return self.plotmesh(var, **kwargs)
        

    def ne2D(self, interactive=False, **kwargs):
        return self.selector2D(self.get('ne'), interactive, **kwargs)

    def te2D(self, interactive=False, **kwargs):
        return self.selectro2D(self.get('te')/self.get('ev'), interactive, **kwargs)

    def ni2D(self, s, interactive=False, **kwargs):
        return self.selector2D(self.get('ni', s), interactive, **kwargs)

    def ti2D(self, interactive=False, **kwargs):
        return self.selector2D(self.get('te')/self.get('ev'), interactive, **kwargs)

    def ng2D(self, s, interactive=False, **kwargs):
            return self.selector2D(self.get('ng', s), interactive, **kwargs)

    def tg2D(self, s, interactive=False, **kwargs):
        return self.selector2D(self.get('tg', s)/self.get('ev'), interactive, **kwargs)

    # TODO: implement masking selector too?
    def CIII_emission_2D(self, fname, interactive=False, **kwargs):
        self.emission_CIII(fname)
        return self.selector2D(self.CIII_emission, interactive, **kwargs)

    def masked_CIII_2D(self, fname, z, maskvalues, interactive, **kwargs):
        self.emission_CIII(fname)
        mask = self.CIII_emission[1:-1,1:-1].reshape(self.nx*self.ny)
        mask = [1*( (x<maskvalues[0]) or (x>maskvalues[1])) for x in mask]
        if interactive is False:
            return self.plotmesh(z, **kwargs)
        else:
            kwargs['mask'] = mask
            kwargs['mvs'] = maskvalues
            return self.variablemaskedmesh(z, **kwargs)

    def CIIImasked_flow(self, fname, maskvalues, interactive=False, **kwargs):
        z = (self.get('upi')**2 + self.get('vy')**2)**0.5
        z = self.get('upi')[:,:,4]
#        return self.masked_CIII_2D(fname, self.get('ne'), maskvalues, **kwargs)
        return self.masked_CIII_2D(fname, z, maskvalues, interactive, **kwargs)



    def variablemaskedmesh(self, z=None, **kwargs):
        from matplotlib.pyplot import ion, ioff, subplots
        from matplotlib.widgets import Slider, RangeSlider
        from matplotlib import is_interactive
        ioff()
        
        f, ax = subplots(figsize=(7,8))
        try:
            kwargs['zrange']
            origrange = kwargs['zrange']
        except:
            kwargs['zrange'] = (z[1:-1,1:-1].min(), z[1:-1,1:-1].max())
            origrange = kwargs['zrange']
        mvs = kwargs.pop('mvs')
        mask = self.CIII_emission[1:-1,1:-1]
        cbar, verts = self.plotmesh(z, ax=ax, watermark=False, 
            interactive=True, **kwargs)
        f.axes[0].set_position([0.125, 0.13, 0.55, 0.85])
        f.axes[1].set_position([0.7, 0.13, 0.82, 0.85])
        mask_position = f.add_axes([0.1, 0.02, 0.6, 0.04])
        mask_slider = RangeSlider(mask_position, 'Mask', mask.min(), 
            mask.max(), valinit=(mvs))
        zrange_position = f.add_axes([0.85, 0.13, 0.04, 0.85])
        zrange_slider = RangeSlider(zrange_position, '', z[1:-1,1:-1].min(), 
            z[1:-1,1:-1].max(), valinit=(origrange), orientation='vertical')

        def update(val):
            from numpy import floor, ceil, linspace
            mask = self.CIII_emission[1:-1,1:-1].reshape(self.nx*self.ny)
            mv = mask_slider.val
            mask = [1*( (x<mv[0]) or (x>mv[1])) for x in mask]
            verts.set_clim(zrange_slider.val)
            verts.set_alpha(mask)

        zrange_slider.on_changed(update)
        mask_slider.on_changed(update)
        f.show() 


        ion()
        return f, zrange_slider, mask_slider


    def variablemesh(self, z=None, **kwargs):
        from matplotlib.pyplot import ion, ioff, subplots
        from matplotlib.widgets import Slider, RangeSlider
        from matplotlib import is_interactive
        ioff()
        
        f, ax = subplots(figsize=(7,8))
        try:
            kwargs['zrange']
            origrange = kwargs['zrange']
        except:
            kwargs['zrange'] = (z[1:-1,1:-1].min(), z[1:-1,1:-1].max())
            origrange = kwargs['zrange']
        cbar, verts = self.plotmesh(z, ax=ax, watermark=False, 
            interactive=True, **kwargs)
        f.axes[0].set_position([0.125, 0.13, 0.55, 0.85])
        f.axes[1].set_position([0.7, 0.13, 0.82, 0.85])
        zrange_position = f.add_axes([0.85, 0.13, 0.04, 0.85])
        zrange_slider = RangeSlider(zrange_position, '', z[1:-1,1:-1].min(), 
            z[1:-1,1:-1].max(), valinit=(origrange), orientation='vertical')

        def update(val):
            from numpy import floor, ceil, linspace
            verts.set_clim(zrange_slider.val)

        zrange_slider.on_changed(update)
        f.show() 


        ion()
        return f, zrange_slider

