# Plotting routines for case objects
from uetools.UePlot import Plot
from uetools.UeUtils import Tools


# TODO: implement profile plots
class Caseplot(Plot):
    def __init__(self, case, *args, **kwargs):
        # Couple get to case
        self.get = case.get
        self.info = case.info
        self.tools = Tools()
        if self.get('geometry')[0].decode('UTF-8').strip() \
            in ['uppersn', 'snull']:
            snull = True
            dnull = False
            if  self.get('geometry')[0].decode('UTF-8').strip() \
                    == 'uppersn':
                usn = True
            else:
                usn = False
        elif  self.get('geometry')[0].decode('UTF-8').strip() \
            == 'dnull':
            snull = False
            usn = False
            dnull = True
        elif  self.get('geometry')[0].decode('UTF-8').strip() \
            == 'dnbot':
            snull = True
            usn = False
            dnull = False
        else:
            raise TypeError("Geometry {} not recognized! Aborting.".format(\
                self.get('geometry')[0].decode('UTF-8').strip()))

        super().__init__(*args, snull=snull, dnull=dnull, usn=usn, **kwargs)
    
    def watermark(self, figure, bottom=0.15, top=0.95, left=0.09, right=0.98):
        """Adds metadata to figure"""
        from time import ctime

        label = '{}, case "{}"\n'.format(ctime(), self.info['casename'])
        label += 'UEDGE {} v{}, UETOOLS v{}, user "{}", hostname "{}"\n'.format(
            self.info['uedge_ver'].replace("$", r"\$"),
            self.info['pyver'],
            self.info['uetoolsversion'],
            self.info['user'],
            self.info['hostname'],
        )
        try:
            label += 'cwd "{}"'.format(self.info['location'])
        except:
            label += 'cwd "{}"'.format(self.info['casefname'])
        figure.subplots_adjust(bottom=bottom, top=top, left=left, right=right)
        figure.text(0.995, 0.005, label, fontsize=4, horizontalalignment="right")

        return


    def it(self, variable, ylabel=None, marksep=True, staggered=False, 
                xlim=(None, None), ylim=(None, None), primary=True, **kwargs
            ):
        """ Plots variable at inner plate as distance along the plate

        In case of double-null, plots values along primary target (bias 
        direction) unless kwarg primary=False is set. For balance double 
        null geometries, the lower target is assumed primary.

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
        trgt = ""
        x = self.get('yyrb')[1:-1]
        idx = 0
        if self.dnull == "upper":
            if not primary:
                idx = self.get("ixlb")[0] + (0 ** staggered)
                trgt = "lower"
                x = self.get('yylb')[1:-1,1]
            else:
                idx = self.get("ixrb")[0]
                trgt = "upper"
                x = self.get('yyrb')[1:-1,0]
        elif self.dnull in ["lower", "balanced"]:
            if not primary:
                idx = self.get("ixrb")[0]
                trgt = "upper"
                x = self.get('yyrb')[1:-1,0]
            else:
                idx = self.get("ixlb")[0] + (0 ** staggered)
                trgt = "lower"
                x = self.get('yyrb')[1:-1,0]
        fig = self.profile(x, variable[idx, 1:-1], **kwargs)
        # Add Sep location if requested
        if marksep is True:
            fig.get_axes()[0].axvline(0, color="grey", linewidth=1)
        fig.get_axes()[0].set_xlim(xlim)
        fig.get_axes()[0].set_ylim(ylim)
        fig.get_axes()[0].set_ylabel(ylabel)
        fig.get_axes()[0].set_xlabel(f"Distance from separatrix along {trgt} IT [m]")
        return fig

    def ot(self, variable, ylabel=None, marksep=True, xlim=(None, None), 
            ylim=(None, None), primary=True, staggered=False, **kwargs
            ):
        """ Plots variable at outer plate as distance along the plate

        In case of double-null, plots values along primary target (bias 
        direction) unless kwarg primary=False is set. For balance double 
        null geometries, the lower target is assumed primary.

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
        trgt = ""
        x = self.get('yyrb')[1:-1]
        idx = -2
        if self.dnull == "upper":
            if primary:
                idx = self.get("ixlb")[1] + (0 ** staggered)
                trgt = "upper"
                x = self.get('yylb')[1:-1,1]
            else:
                idx = self.get("ixrb")[1]
                trgt = "lower"
                x = self.get('yyrb')[1:-1,0]
        elif self.dnull in ["lower", "balanced"]:
            if primary:
                idx = self.get("ixrb")[1]
                trgt = "lower"
                x = self.get('yyrb')[1:-1,1]
            else:
                idx = self.get("ixlb")[1] + (0 ** staggered)
                trgt = "lower"
                x = self.get('yyrb')[1:-1,1]
        fig = self.profile(x, variable[idx, 1:-1], **kwargs)
        # Add Sep location if requested
        if marksep is True:
            fig.get_axes()[0].axvline(0, color="grey", linewidth=1)
        fig.get_axes()[0].set_ylabel(ylabel)
        fig.get_axes()[0].set_xlim(xlim)
        fig.get_axes()[0].set_ylim(ylim)
        fig.get_axes()[0].set_xlabel(f"Distance from separatrix along {trgt} OT [m]")
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
        fig = self.profile(self.get('yyc')[1:-1], variable[self.get('ixmp'), 1:-1], **kwargs)
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

        fig = self.profile(self.get('yyc')[1:-1], variable[ix, 1:-1], **kwargs)
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
        # TODO: how to deal with dnulls??

        from numpy import concatenate
        # Check if we are in the PFR
        ixpt1 = self.get('ixpt1')[0]
        ixpt2 = self.get('ixpt2')[0]
        if (self.get('ixp1')[ixpt1, iy] - ixpt1) == 1:
            x = self.lcon[1:-1, iy]
            fig = self.profile(x, variable[1:-1, iy], **kwargs)
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
            fig = self.profile(x, y, **kwargs)
            if markxpts is True:
                xpt = 0.5*(self.lcon[ixpt1,iy] + self.lcon[ixpt2+1, iy])
                fig.get_axes()[0].axvline(xpt, color="grey", linewidth=1)

        return fig

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
    def grid(self, **kwargs):
        return self.mesh(**kwargs)
    
    def gridue(self, file, flip=True, **kwargs):
        from h5py import is_hdf5

        if not is_hdf5(file):
            raise NotImplementedError("ASCII gridue support not implemented.")
        rm=self.tools.hdf5search(file,'rm')[1:-1,1:-1]
        zm=self.tools.hdf5search(file,'zm')[1:-1,1:-1]
        if (
            self.get("geometry")[0].strip().lower().decode("UTF-8") == "uppersn"
        ) and (flip is True):
            zm = self.disp-zm        
        return self.mesh( 
            lcfs=False,
            rm=rm,
            zm=zm,
            flip=flip,
            **kwargs
            )

    def heatmap(self, var, s=None, **kwargs):
        if isinstance(var, str):
            return self.mesh(self.get(var, s), **kwargs)
        else:
            return self.mesh(var, **kwargs)

    def selector2D(self, var, interactive, **kwargs):
        if interactive is True:
            return self.variablemesh(var, **kwargs)
        else:
            return self.mesh(var, **kwargs)

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
            return self.mesh(z, **kwargs)
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
        cbar, verts = self.mesh(
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

    def OMP_chunks(self, seed=1):
        from matplotlib.pyplot import subplots, get_cmap
        from matplotlib.patches import Polygon
        from numpy import linspace
        import random

        f, ax = subplots()
        ax.set_aspect('equal')
        cmap=get_cmap('jet')
        cols = linspace(0,1,self.get('Nchunks'))
        random.Random(seed).shuffle(cols)
        colors = iter(cmap(cols))
        nx = self.get('nx')
        ny = self.get('ny')
        iysptrx1 = self.get('iysptrx1')
        iysptrx2 = self.get('iysptrx2')
        ixpt1 = self.get('ixpt1')
        ixpt2 = self.get('ixpt2')
        ixrb = self.get('ixrb')
        rangexptchunk = self.get('rangexptchunk')
        nxpt = self.get('nxpt')

        # 
        ax.plot((0.5,0.5),(ny+.5,0.5),'k-', linewidth=1.5)
        ax.plot((nx+0.5,nx+0.5),(ny+.5,0.5),'k-', linewidth=1.5)
        ax.plot((0.5,nx+.5),(.5,0.5),'k-', linewidth=1.5)
        ax.plot((0.5,nx+0.5),(ny+.5,ny+0.5),'k-', linewidth=1.5)
        if nxpt > 1:
            ax.plot((ixrb[0]+1.5, ixrb[0]+1.5), (-0.5,ny+1.5), 'k:', linewidth=.5)
            ax.plot((ixrb[0]+0.5, ixrb[0]+0.5), (-0.5,ny+1.5), 'k-', linewidth=1.5)
            ax.plot((ixrb[0]+2.5, ixrb[0]+2.5), (-0.5,ny+1.5), 'k-', linewidth=1.5)

        ax.plot((-0.5,-0.5),(ny+1.5,-0.5),'k--', linewidth=.5)
        ax.plot((nx+1.5,nx+1.5),(ny+1.5,-0.5),'k--', linewidth=.5)
        ax.plot((-0.5,nx+1.5),(-.5,-0.5),'k--', linewidth=.5)
        ax.plot((-0.5,nx+1.5),(ny+1.5,ny+1.5),'k--', linewidth=.5)

        for ix in range(1, nx+1):
            ax.plot((ix+0.5,ix+0.5), (-0.5, ny+1.5), 'k:', linewidth=1)
        for iy in range(1, ny+1):
            ax.plot((-0.5,nx+1.5), (iy+0.5, iy+.5), 'k:', linewidth=1)

        for ixpt in range(nxpt):
            ax.plot((-0.5,nx+1.5),(iysptrx1[ixpt]+.5, iysptrx1[ixpt]+.5), 'k-', linewidth=1.5)
            ax.plot((ixpt1[ixpt]+0.5,ixpt1[ixpt]+0.5), (-0.5,iysptrx1[ixpt]+.5), 'k-', linewidth=1.5)
            ax.plot((ixpt2[ixpt]+0.5,ixpt2[ixpt]+0.5), (-0.5,iysptrx2[ixpt]+.5), 'k-', linewidth=1.5)

        for ixpt in range(nxpt):
            for icut in range(self.get('Nxptchunks')[ixpt]):
                for ii in range(2):
                    xpatch = rangexptchunk[ixpt, ii, icut]
                    ax.add_patch(Polygon( [(xpatch[0]-.5,xpatch[2]-.5),(xpatch[0]-.5,xpatch[3]+.5),(xpatch[1]+.5, xpatch[3]+.5),(xpatch[1]+.5, xpatch[2]-0.5)],
                    closed=True, edgecolor='k', facecolor='red'))

        for xpatch in self.get('rangechunk'):
            ax.add_patch(Polygon( [(xpatch[0]-.5,xpatch[2]-.5),(xpatch[0]-.5,xpatch[3]+.5),(xpatch[1]+.5, xpatch[3]+.5),(xpatch[1]+.5, xpatch[2]-0.5)],
                 closed=True, edgecolor='k', alpha=0.3, linestyle='--', linewidth=2, color=next(colors)))


        f.show()


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
        cbar, verts = self.mesh(
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

    def iongradB(self, ax=None, width=0.02, color='k', flip=True,
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
        if self.usn and (flip is True):
            zm = -zm + self.disp 
        y0 = mean(zm[self.get('ixpt1')[0]+1:self.get('ixpt2')[0]+1, 0, 0])

        x = pol * self.eastnormaln[0] + rad * self.northnormaln[0]
        y = pol * self.eastnormaln[1] + rad * self.northnormaln[1]
        if (flip is True) and self.usn:
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

        cbar, verts = self.mesh(
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
