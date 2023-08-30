# Plotting routines for case objects
from uetools.UePlot import Plot


class Caseplot(Plot):
    def __init__(self, **kwargs):
#        super(Caseplot, self).__init__()
        super().__init__(**kwargs)
        self.createvertices(self.get('rm'), self.get('zm')) 
        # TODO: initialize parent properly instead
        # of explicit call - how??
        # Calculate distances along targets
        self.otdistance = self.get("yyrb")
        self.itdistance = self.get("yylb")
        return

    # TODO: implement profile plots

    # TODO: fix this function
    def radialdistance(self, ind):
        from numpy import cumsum

        # Coordinates of edge of grid
        dr0 = self.get("rm", 2)[ind, 1]
        dz0 = self.get("zm", 2)[ind, 1]
        # Halfway point of cell centers
        drplatecenter = 0.5 * (self.get("rm", 2)[ind] + self.get("rm", 4)[ind])
        dzplatecenter = 0.5 * (self.get("zm", 2)[ind] + self.get("zm", 4)[ind])
        # Calculate strike point coordinates
        iysptrx = self.get("iysptrx")
        rsp = self.get("rm", 4)[ind, iysptrx] - drplatecenter[iysptrx]
        zsp = self.get("zm", 4)[ind, iysptrx] - dzplatecenter[iysptrx]
        # Normalize to grid edge
        drplatecenter -= dr0
        dzplatecenter -= dz0
        # Calculate cumulative distance along target
        dist = cumsum((drplatecenter ** 2 + dzplatecenter ** 2) ** 0.5)
        # Calculate distance along plate to strike point
        dsep = dist[iysptrx] + (rsp ** 2 + zsp ** 2) ** 0.5
        return dist - dsep

    def it(self, variable, marksep=True, staggered=False, **kwargs):
        # Exchange YYC for working radialdistance
        fig = self.plotprofile(
            self.itdistance[1:-1], variable[0 ** staggered, 1:-1], **kwargs
        )
        # Add Sep location if requested
        if marksep is True:
            fig.get_axes()[0].axvline(0, color="g", linewidth=1)
        return fig

    def ot(self, variable, marksep=True, **kwargs):
        # Exchange YYC for working radialdistance
        fig = self.plotprofile(self.otdistance[1:-1], variable[-2, 1:-1], **kwargs)
        # Add Sep location if requested
        if marksep is True:
            fig.get_axes()[0].axvline(0, color="g", linewidth=1)
        return fig

    def omp(self):
        return

    def imp(self):
        return

    def row(self):
        return

    def neOT(self, s=None, **kwargs):
        f = self.ot(self.get('ne', s), **kwargs)
        f.get_axes()[0].set_ylabel(r'$\rm n_e^{LFS-t}~[m^{-3}]$')
        f.get_axes()[0].set_xlabel('Distance along LFS-t [m]')
        return f

    def teOT(self, s=None, **kwargs):
        f = self.ot(self.get('te', s)/self.get('ev'), **kwargs)
        f.get_axes()[0].set_ylabel(r'$\rm T_e^{LFS-t}~[eV]$')
        f.get_axes()[0].set_xlabel('Distance along LFS-t [m]')
        return f

    def niOT(self, s=None, **kwargs):
        f = self.ot(self.get('ni', s), **kwargs)
        f.get_axes()[0].set_ylabel(r'$\rm n_i^{LFS-t}~[m^{-3}]$')
        f.get_axes()[0].set_xlabel('Distance along LFS-t [m]')
        return f

    def tiOT(self, s=None, **kwargs):
        f = self.ot(self.get('ti', s)/self.get('ev'), **kwargs)
        f.get_axes()[0].set_ylabel(r'$\rm T_i^{LFS-t}~[eV]$')
        f.get_axes()[0].set_xlabel('Distance along LFS-t [m]')
        return f

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
        return self.selector2D(self.get("ne"), interactive, **kwargs)

    def te2D(self, interactive=False, **kwargs):
        return self.selector2D(self.get("te") / self.get("ev"), interactive, **kwargs)

    def ni2D(self, s, interactive=False, **kwargs):
        return self.selector2D(self.get("ni", s), interactive, **kwargs)

    def ti2D(self, interactive=False, **kwargs):
        return self.selector2D(self.get("ti") / self.get("ev"), interactive, **kwargs)

    def ng2D(self, s, interactive=False, **kwargs):
        return self.selector2D(self.get("ng", s), interactive, **kwargs)

    def tg2D(self, s, interactive=False, **kwargs):
        return self.selector2D(
            self.get("tg", s) / self.get("ev"), interactive, **kwargs
        )

    # TODO: implement masking selector too?
    def CIII_emission_2D(self, fname, interactive=False, **kwargs):
        self.emission_CIII(fname)
        return self.selector2D(self.CIII_emission, interactive, **kwargs)

    def masked_CIII_2D(self, fname, z, maskvalues, interactive, **kwargs):
        self.emission_CIII(fname)
        mask = self.CIII_emission[1:-1, 1:-1].reshape(self.nx * self.ny)
        mask = [1 * ((x < maskvalues[0]) or (x > maskvalues[1])) for x in mask]
        if interactive is False:
            return self.plotmesh(z, **kwargs)
        else:
            kwargs["mask"] = mask
            kwargs["mvs"] = maskvalues
            return self.variablemaskedmesh(z, **kwargs)

    def CIIImasked_flow(self, fname, maskvalues, interactive=False, species=4,
        **kwargs):
        z = (self.get("upi") ** 2 + self.get("vy") ** 2) ** 0.5
        z = self.get("upi")[:, :, species]
        #        return self.masked_CIII_2D(fname, self.get('ne'), maskvalues, **kwargs)
        return self.masked_CIII_2D(fname, z, maskvalues, interactive, **kwargs)
    
    def i_vel(self,pec_fname,i_idx,wavel,floor=0,plot=True,zrange=(-30,30),**kwargs):
        from numpy import nan
        self.i_emiss(pec_fname,i_idx,wavel,plot=False)
        t = self.get_tname(i_idx)
        emiss = self.__getattribute__("{}_emission".format(t))
        try:
            v = self.get("upi")[:,:,i_idx]/1000 #km/s
        except:
            v = self.get("up")[:,:,i_idx]/1000 #km/s
        v[emiss<floor*emiss.max()] = nan
        if plot:
            title = "{} velocity".format(t)
            return self.plotmesh(v,cmap='seismic',zrange=zrange,**kwargs)
        else:
            return v*1000
        

    def variablemaskedmesh(self, z=None, **kwargs):
        from matplotlib.pyplot import ion, ioff, subplots
        from matplotlib.widgets import Slider, RangeSlider
        from matplotlib import is_interactive

        ioff()

        f, ax = subplots(figsize=(7, 8))
        try:
            kwargs["zrange"]
            origrange = kwargs["zrange"]
        except:
            kwargs["zrange"] = (z[1:-1, 1:-1].min(), z[1:-1, 1:-1].max())
            origrange = kwargs["zrange"]
        mvs = kwargs.pop("mvs")
        mask = self.CIII_emission[1:-1, 1:-1]
        cbar, verts = self.plotmesh(
            z, ax=ax, watermark=False, interactive=True, **kwargs
        )
        f.axes[0].set_position([0.125, 0.13, 0.55, 0.85])
        f.axes[1].set_position([0.7, 0.13, 0.82, 0.85])
        mask_position = f.add_axes([0.1, 0.02, 0.6, 0.04])
        mask_slider = RangeSlider(
            mask_position, "Mask", mask.min(), mask.max(), valinit=(mvs)
        )
        zrange_position = f.add_axes([0.85, 0.13, 0.04, 0.85])
        zrange_slider = RangeSlider(
            zrange_position,
            "",
            z[1:-1, 1:-1].min(),
            z[1:-1, 1:-1].max(),
            valinit=(origrange),
            orientation="vertical",
        )

        def update(val):
            from numpy import floor, ceil, linspace

            mask = self.CIII_emission[1:-1, 1:-1].reshape(self.nx * self.ny)
            mv = mask_slider.val
            mask = [1 * ((x < mv[0]) or (x > mv[1])) for x in mask]
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

        f, ax = subplots(figsize=(7, 8))
        try:
            kwargs["zrange"]
            origrange = kwargs["zrange"]
        except:
            kwargs["zrange"] = (z[1:-1, 1:-1].min(), z[1:-1, 1:-1].max())
            origrange = kwargs["zrange"]
        cbar, verts = self.plotmesh(
            z, ax=ax, watermark=False, interactive=True, **kwargs
        )
        f.axes[0].set_position([0.125, 0.13, 0.55, 0.85])
        f.axes[1].set_position([0.7, 0.13, 0.82, 0.85])
        zrange_position = f.add_axes([0.85, 0.13, 0.04, 0.85])
        zrange_slider = RangeSlider(
            zrange_position,
            "",
            z[1:-1, 1:-1].min(),
            z[1:-1, 1:-1].max(),
            valinit=(origrange),
            orientation="vertical",
        )

        def update(val):
            from numpy import floor, ceil, linspace

            verts.set_clim(zrange_slider.val)

        zrange_slider.on_changed(update)
        f.show()

        ion()
        return f, zrange_slider

    def ionvel(self, s, **kwargs):
        return self.plot_streamline("upi", "vy", s, surfnorm=False, **kwargs)

    def gasvel(self, s, **kwargs):
        return self.plot_streamline("uug", "vyg", s, surfnorm=False, **kwargs)

    def ionflow(self, s, surfnorm=True, **kwargs):
        return self.plot_streamline("fnix", "fniy", s, surfnorm, **kwargs)

    def gaseneflow(self, s, surfnorm=True, **kwargs):
        return self.plot_streamline("fegx", "fegy", s, surfnorm, **kwargs)

    def eleceneflow(self, s, surfnorm=True, **kwargs):
        return self.plot_streamline("feex", "feey", s, surfnorm, **kwargs)

    def ioneneflow(self, s, surfnorm=True, **kwargs):
        return self.plot_streamline("feix", "feiy", s, surfnorm, **kwargs)

    def gasflow(self, s, surfnorm=True, **kwargs):
        return self.plot_streamline("fngx", "fngy", s, surfnorm, **kwargs)

    def plot_driftdirection(self, ax=None, width=0.02, color='k', **kwargs):
        ''' Plots the drift direction on the requested axis '''
        from numpy import sum, mean
        if ax is None:
            f = self.grid()
            ax = f.get_axes()[0]
        pol = (self.get('v2cb')[:,:,0]*self.get('rbfbt'))
        rad = (self.get('vycb'))[:,:,0]
        x0 = mean(self.get('rm')[self.get(\
            'ixpt1')[0]+1:self.get('ixpt2')[0]+1, 0, 0])
        y0 = mean(self.get('zm')[self.get(\
            'ixpt1')[0]+1:self.get('ixpt2')[0]+1, 0, 0])

        x = pol * self.eastnormaln[0] + rad * self.northnormaln[0]
        y = pol * self.eastnormaln[1] + rad * self.northnormaln[1]
        x = sum(x)
        y = sum(y)
        norm = max(abs(x), abs(y))
        x /= norm
        y /= norm
        ax.arrow(x0, y0, x/3, y/3, width=width, color=color, **kwargs)


    def plot_streamline(self, varpol, varrad, s=None, surfnorm=True, **kwargs):
        from numpy import zeros_like
        pol = self.get(varpol, s) / self.get('sx')**surfnorm
        rad = self.get(varrad, s) / self.get('sy')**surfnorm

        pol_use = zeros_like(pol)
        rad_use = zeros_like(rad)

        return self.streamline(pol, rad, **kwargs)

    def plot_2Dyldot(self, **kwargs):
        """Returns a series of figures to scroll through"""
        from matplotlib.pyplot import subplots, ion, ioff
        from matplotlib.widgets import Slider

        ioff()
        f, ax = subplots(figsize=(7, 8))

        vararray = self.numvararr((self.get("yldot") * self.get("sfscal"))[:-2])

        try:
            kwargs["zrange"]
            origrange = max([abs(x) for x in kwargs["zrange"]])
        except:
            kwargs["zrange"] = (
                vararray[1:-1, 1:-1, :].min(),
                vararray[1:-1, 1:-1, :].max(),
            )
            origrange = max([abs(x) for x in kwargs["zrange"]])
        kwargs["zrange"] = (-origrange, origrange)
        try:
            kwargs["cmap"]
        except:
            kwargs["cmap"] = "bwr"

        cbar, verts = self.plotmesh(
            vararray[:, :, 0], ax=ax, watermark=False, interactive=True, **kwargs
        )
        f.axes[0].set_position([0.125, 0.13, 0.55, 0.85])
        f.axes[1].set_position([0.7, 0.13, 0.82, 0.85])
        slice_position = f.add_axes([0.1, 0.02, 0.65, 0.04])
        slice_slider = Slider(slice_position, "", 0, vararray.shape[-1], valstep=1)
        zrange_position = f.add_axes([0.85, 0.13, 0.04, 0.85])
        zrange_slider = Slider(
            zrange_position,
            "",
            0.0,
            vararray[1:-1, 1:-1, :].max(),
            valinit=origrange,
            orientation="vertical",
        )

        def update(val):
            from numpy import where

            slce = slice_slider.val
            verts.set_array(
                vararray[1:-1, 1:-1, slce - 1].reshape(
                    (self.get("nx") * self.get("ny"))
                )
            )
            verts.set_clim((-zrange_slider.val, zrange_slider.val))

        slice_slider.on_changed(update)
        zrange_slider.on_changed(update)

        f.show()
        ion()
        return f, (slice_slider, zrange_slider)
