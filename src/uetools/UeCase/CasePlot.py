# Plotting routines for case objects
from uetools.UePlot import Plot


class Caseplot(Plot):
    # TODO: implement profile plots

    def it(self, variable, ylabel=None, marksep=True, staggered=False, 
                xlim=(None, None), ylim=(None, None), **kwargs
            ):
        """ Plots variable at inner plate as distance along the plate

        Arguments
        ---------
        variable : ndarray
            2D array with dimension (nx+2, ny+2) to be plotted

        Keyword arguments
        -----------------
        ylabel : str (default = None)
            Label for Y-axis
        xlim : tuple of floats (default = (None, None)
            X-axis limits
        ylim : tuple of floats (default = (None, None)
            Y-axis limits
        marksep : bool (default = True)
            Marks the separatrix with a vertical line if True
        staggered : bool (default = False)
            Offsets the variable cell by -1 to account for staggered grid
            for variable (such velocities, fluxes, etc.) if True
            
        Returns
        -------
        Figure 
        """
        fig = self.plotprofile(
            self.get('yylb')[1:-1], variable[0 ** staggered, 1:-1], **kwargs
        )
        # Add Sep location if requested
        if marksep is True:
            fig.get_axes()[0].axvline(0, color="grey", linewidth=1)
        fig.get_axes()[0].set_xlim(xlim)
        fig.get_axes()[0].set_ylim(ylim)
        fig.get_axes()[0].set_ylabel(ylabel)
        fig.get_axes()[0].set_xlabel("Distance from separatrix along IT [m]")
        return fig

    def ot(self, variable, ylabel=None, marksep=True, xlim=(None, None), 
            ylim=(None, None), **kwargs
            ):
        """ Plots variable at outer plate as distance along the plate

        Arguments
        ---------
        variable : ndarray
            2D array with dimension (nx+2, ny+2) to be plotted

        Keyword arguments
        -----------------
        ylabel : str (default = None)
            Label for Y-axis
        xlim : tuple of floats (default = (None, None)
            X-axis limits
        ylim : tuple of floats (default = (None, None)
            Y-axis limits
        marksep : bool (default = True)
            Marks the separatrix with a vertical line if True
            
        Returns
        -------
        Figure 
        """
        fig = self.plotprofile(self.get('yyrb')[1:-1], variable[-2, 1:-1], **kwargs)
        # Add Sep location if requested
        if marksep is True:
            fig.get_axes()[0].axvline(0, color="grey", linewidth=1)
        fig.get_axes()[0].set_ylabel(ylabel)
        fig.get_axes()[0].set_xlim(xlim)
        fig.get_axes()[0].set_ylim(ylim)
        fig.get_axes()[0].set_xlabel("Distance from separatrix along OT [m]")
        return fig

    def omp(self, variable, ylabel=None, marksep=True, xlim=(None, None), 
                ylim=(None, None), **kwargs):
        """ Plots variable at OMP as function of distance from separatrix

        Arguments
        ---------
        variable : ndarray
            2D array with dimension (nx+2, ny+2) to be plotted

        Keyword arguments
        -----------------
        ylabel : str (default = None)
            Label for Y-axis
        xlim : tuple of floats (default = (None, None)
            X-axis limits
        ylim : tuple of floats (default = (None, None)
            Y-axis limits
        marksep : bool (default = True)
            Marks the separatrix with a vertical line if True
            
        Returns
        -------
        Figure
        """
        fig = self.plotprofile(self.get('yyc')[1:-1], variable[self.get('ixmp'), 1:-1], **kwargs)
        # Add Sep location if requested
        if marksep is True:
            fig.get_axes()[0].axvline(0, color="grey", linewidth=1)
        fig.get_axes()[0].set_xlabel("Distance from separatrix along OMP [m]")
        fig.get_axes()[0].set_ylabel(ylabel)
        return fig

    def row(self, variable, ix, marksep=True, ylabel=None, 
                xlim=(None, None), ylim=(None, None), 
                **kwargs
            ):
        """ Plots variable at ix vs distance from separatrix projected to OMP

        Plots the radial profile at variable at poloidal index ix as a function
        of distance from separatrix when plotted to the outer midplane. 

        Arguments
        ---------
        variable : ndarray
            2D array with dimension (nx+2, ny+2) to be plotted
        ix : int
            Poloidal index of radial profile to be plotted

        Keyword arguments
        -----------------
        ylabel : str (default = None)
            Label for Y-axis
        xlim : tuple of floats (default = (None, None)
            X-axis limits
        ylim : tuple of floats (default = (None, None)
            Y-axis limits
        marksep : bool (default = True)
            Marks the separatrix with a vertical line if True
            
        Returns
        -------
        Figure
        """
        fig = self.plotprofile(self.get('yyc')[1:-1], variable[ix, 1:-1], **kwargs)
        # Add Sep location if requested
        if marksep is True:
            fig.get_axes()[0].axvline(0, color="grey", linewidth=1)
        fig.get_axes()[0].set_ylabel(ylabel)
        fig.get_axes()[0].set_xlim(xlim)
        fig.get_axes()[0].set_ylim(ylim)
        fig.get_axes()[0].set_xlabel("Distance from separatrix "+\
            "projected to OMP [m]")
        return fig

    def ft(self, variable, iy, markxpts=True, **kwargs):
        """ Plots variable at iy vs distance from separatrix projected to OMP

        Plots the radial profile at variable at poloidal index ix as a function
        of distance from separatrix when plotted to the outer midplane. 

        Arguments
        ---------
        variable : ndarray
            2D array with dimension (nx+2, ny+2) to be plotted
        ix : int
            Poloidal index of radial profile to be plotted

        Keyword arguments
        -----------------
        ylabel : str (default = None)
            Label for Y-axis
        xlim : tuple of floats (default = (None, None)
            X-axis limits
        ylim : tuple of floats (default = (None, None)
            Y-axis limits
        marksep : bool (default = True)
            Marks the separatrix with a vertical line if True
            
        Returns
        -------
        Figure
        """
        from numpy import concatenate
        # Check if we are in the PFR
        ixpt1 = self.get('ixpt1')[0]
        ixpt2 = self.get('ixpt2')[0]
        if (self.get('ixp1')[ixpt1, iy] - ixpt1) == 1:
            x = self.lcon[1:-1, iy]
            fig = self.plotprofile(x, variable[1:-1, iy], **kwargs)
            # Add Sep location if requested
            if markxpts is True:
                xpt1 = 0.5*sum(self.lcon[ixpt1:ixpt1+2, iy])
                xpt2 = 0.5*sum(self.lcon[ixpt2:ixpt2+2, iy])
                fig.get_axes()[0].axvline(xpt1, color="grey", linewidth=1)
                fig.get_axes()[0].axvline(xpt2, color="grey", linewidth=1)
        else:
            xpt = self.lcon[ixpt1, iy]  
            x = concatenate((self.lcon[:ixpt1+1, iy], self.lcon[ixpt2+1:, iy]+xpt))
            y = concatenate((variable[:ixpt1+1, iy], variable[ixpt2+1:, iy]))
            fig = self.plotprofile(x, y, **kwargs)
            if markxpts is True:
                xpt = 0.5*(self.lcon[ixpt1,iy] + self.lcon[ixpt2+1, iy])
                fig.get_axes()[0].axvline(xpt, color="grey", linewidth=1)

        return fig

    def sep(self, variable, markxpts=True, **kwargs):
        """

        Arguments
        ---------

        Keyword arguments
        -----------------
            
        Modifies
        --------

        Returns
        -------
        """
        return self.ft(variable, self.get('iysptrx')+1, markxpts=markxpts, **kwargs)
        


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

    def ngOT(self, s=None, **kwargs):
        f = self.ot(self.get('ng', s), **kwargs)
        f.get_axes()[0].set_ylabel(r'$\rm n_i^{LFS-t}~[m^{-3}]$')
        f.get_axes()[0].set_xlabel('Distance along LFS-t [m]')
        return f

    def neIT(self, s=None, **kwargs):
        f = self.it(self.get('ne', s), **kwargs)
        f.get_axes()[0].set_ylabel(r'$\rm n_e^{HFS-t}~[m^{-3}]$')
        f.get_axes()[0].set_xlabel('Distance along HFS-t [m]')
        return f

    def teIT(self, s=None, **kwargs):
        f = self.it(self.get('te', s)/self.get('ev'), **kwargs)
        f.get_axes()[0].set_ylabel(r'$\rm T_e^{HFS-t}~[eV]$')
        f.get_axes()[0].set_xlabel('Distance along HFS-t [m]')
        return f

    def niIT(self, s=None, **kwargs):
        f = self.it(self.get('ni', s), **kwargs)
        f.get_axes()[0].set_ylabel(r'$\rm n_i^{HFS-t}~[m^{-3}]$')
        f.get_axes()[0].set_xlabel('Distance along HFS-t [m]')
        return f

    def tiIT(self, s=None, **kwargs):
        f = self.it(self.get('ti', s)/self.get('ev'), **kwargs)
        f.get_axes()[0].set_ylabel(r'$\rm T_i^{HFS-t}~[eV]$')
        f.get_axes()[0].set_xlabel('Distance along HFS-t [m]')
        return f

    def ngIT(self, s=None, **kwargs):
        f = self.it(self.get('ng', s), **kwargs)
        f.get_axes()[0].set_ylabel(r'$\rm n_i^{LFS-t}~[m^{-3}]$')
        f.get_axes()[0].set_xlabel('Distance along LFS-t [m]')
        return f



    # Expand the 2D plot list
    def grid(self, gridue=None, **kwargs):
        return self.plotmesh(gridue, **kwargs)

    def heatmap(self, var, s=None, **kwargs):
        if isinstance(var, str):
            return self.plotmesh(self.get(var, s), **kwargs)
        else:
            return self.plotmesh(var, **kwargs)

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
        """

        Arguments
        ---------

        Keyword arguments
        -----------------
            
        Modifies
        --------

        Returns
        -------
        """
        self.emission_CIII(fname)
        return self.selector2D(self.CIII_emission, interactive, **kwargs)

    def masked_CIII_2D(self, fname, z, maskvalues, interactive, **kwargs):
        """

        Arguments
        ---------

        Keyword arguments
        -----------------
            
        Modifies
        --------

        Returns
        -------
        """
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
        """

        Arguments
        ---------

        Keyword arguments
        -----------------
            
        Modifies
        --------

        Returns
        -------
        """
        z = (self.get("upi") ** 2 + self.get("vy") ** 2) ** 0.5
        z = self.get("upi")[:, :, species]
        #        return self.masked_CIII_2D(fname, self.get('ne'), maskvalues, **kwargs)
        return self.masked_CIII_2D(fname, z, maskvalues, interactive, **kwargs)

    def variablemaskedmesh(self, z=None, **kwargs):
        """

        Arguments
        ---------

        Keyword arguments
        -----------------
            
        Modifies
        --------

        Returns
        -------
        """
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
        """

        Arguments
        ---------

        Keyword arguments
        -----------------
            
        Modifies
        --------

        Returns
        -------
        """
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

    def plot_driftdirection(self, ax=None, width=0.02, color='k', flip=False,
        **kwargs):
        ''' Plots the drift direction on the requested axis '''
        """

        Arguments
        ---------

        Keyword arguments
        -----------------
            
        Modifies
        --------

        Returns
        -------
        """
        from numpy import sum, mean
        if ax is None:
            f = self.grid()
            ax = f.get_axes()[0]
        pol = (self.get('v2cb')[:,:,0]*self.get('rbfbt'))
        rad = (self.get('vycb'))[:,:,0]
        x0 = mean(self.get('rm')[self.get(\
            'ixpt1')[0]+1:self.get('ixpt2')[0]+1, 0, 0])
        zm = self.get('zm')
        if flip is True:
            zm = -zm + self.disp 
        y0 = mean(zm[self.get('ixpt1')[0]+1:self.get('ixpt2')[0]+1, 0, 0])

        x = pol * self.eastnormaln[0] + rad * self.northnormaln[0]
        y = pol * self.eastnormaln[1] + rad * self.northnormaln[1]
        if flip is True:
            y *= -1
        x = sum(x)
        y = sum(y)
        norm = max(abs(x), abs(y))
        x /= norm
        y /= norm
        ax.arrow(x0, y0, x/3, y/3, width=width, color=color, **kwargs)


    def plot_streamline(self, varpol, varrad, s=None, surfnorm=True, **kwargs):
        """

        Arguments
        ---------

        Keyword arguments
        -----------------
            
        Modifies
        --------

        Returns
        -------
        """
        from numpy import zeros_like
        pol = self.get(varpol, s) / self.get('sx')**surfnorm
        rad = self.get(varrad, s) / self.get('sy')**surfnorm

        pol_use = zeros_like(pol)
        rad_use = zeros_like(rad)

        return self.streamline(pol, rad, **kwargs)

    def plot_2Dyldot(self, **kwargs):
        """Returns a series of figures to scroll through"""
        """

        Arguments
        ---------

        Keyword arguments
        -----------------
            
        Modifies
        --------

        Returns
        -------
        """
        from matplotlib.pyplot import subplots, ion, ioff
        from matplotlib.widgets import Slider
        from numpy import array, transpose

        ioff()
        f, ax = subplots(figsize=(7, 8))
        vararray = transpose(
            array((self.get("yldot") * self.get("sfscal"))[:-2])
            .reshape((self.get("ny") + 2, self.get("nx") + 2, self.get("numvar")))
            .T,
            (1, 2, 0),
        )
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
