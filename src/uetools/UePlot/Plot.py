# Object for plotting
from matplotlib.pyplot import ion

ion()

# TODO: implement divergence plotting/calculation

class Plot:
    def __init__(self, *args, rm=None, zm=None, snull=True, usn=False, dnull=False, **kwargs):
        """Constructs patches objects
        rm - UEDGE R-node object
        zm - UEDGE Z-node object
        """
        # Intialize empty
        if (rm is None) or (zm is None):
            try:
                rm = self.get("rm")
                zm = self.get("zm")
            except:
                return
        self.snull = snull
        self.dnull = dnull
        self.usn = usn

        if self.snull is False:
            sep = self.get("iysptrx1")
            if sep[0] == sep[1]:
                self.dnull = "balanced"
            elif sep[0] < sep[1]:
                self.dnull = "lower"
            else:
                self.dnull = "upper"


        self.createvertices(rm, zm)
        super().__init__(*args, **kwargs)
        return

    def getomit(self, var):
        """ Helper function to handled partial grids w/ omits """
        nyomit = self.get('nyomitmx', verbose=False)
        nxomit = self.get('nxomit', verbose = False)
        if nxomit is None:
            nxomit = 0
        if nyomit is None:
            nyomit = 0
        if isinstance(var, str):
            var = self.get(var)
        if nyomit > 0:
            var = var[:,:-nyomit]
        return var[nxomit:]

    def createvertices(self, rm, zm):
        from numpy import zeros, transpose, cross, sum
        # CREATE POLYGON COLLECTIONS TO USE
        self.vertices = self.createpolycollection(self.getomit(rm), self.getomit(zm))
        self.disp=0
        if self.get("geometry")[0].strip().lower().decode("UTF-8") == "uppersn":
            self.disp = 0
            if self.get("rmagx") + self.get("zmagx") == 0:
                # Normalizing to -min(-zm) results in "jumping" USN cases when
                # the core surfaces change: revert to using set 2.8m displacement
                self.disp = 2.8
            else:
                self.disp = 2 * self.get("zmagx")
            self.uppersnvertices = self.createpolycollection(
                self.getomit(rm), self.getomit(self.disp -zm), setparams=False
            )

        # CREATE NORMAL VECTOR IN LOCAL CELL COORDINATES
        nodes = zeros((self.get("nx") + 2, self.get("ny") + 2, 5, 2))
        nodes[:, :, :, 0] = self.getomit("rm")
        nodes[:, :, :, 1] = self.getomit("zm")
        nodes = transpose(nodes, (2, 3, 0, 1))

        # TODO: rather than align poloidal/radial in direction of cell,
        # evaluate face normal locally at cell face?
        # Find midpoints of y-faces
        self.symid = zeros((2, 2, self.get("nx") + 2, self.get("ny") + 2))
        self.symid[0] = (nodes[2] + nodes[1]) / 2  # Lower face center
        self.symid[1] = (nodes[4] + nodes[3]) / 2  # Upper face center
        # Find midpoints of x-faces
        self.sxmid = zeros((2, 2, self.get("nx") + 2, self.get("ny") + 2))
        self.sxmid[0] = (nodes[3] + nodes[1]) / 2  # Left face center
        self.sxmid[1] = (nodes[4] + nodes[2]) / 2  # Right face center

        # Find vectors of east faces
        eastface = zeros((3, self.get("nx") + 2, self.get("ny") + 2))
        eastface[:-1] = nodes[4] - nodes[2]
        # Find vectors of north faces
        northface = zeros((3, self.get("nx") + 2, self.get("ny") + 2))
        northface[:-1] = nodes[4] - nodes[3]
        # Find normals to faces
        toroidal = zeros((3, self.get("nx") + 2, self.get("ny") + 2))
        toroidal[-1] = 1
        eastnormal = cross(eastface, toroidal, axis=0)
        northnormal = cross(toroidal, northface, axis=0)

        self.northnormaln = zeros((2, self.get("nx") + 2, self.get("ny") + 2))
        for i in range(2):
            self.northnormaln[i] = northnormal[i] / (sum(northnormal**2, 
                axis=0) ** 0.5 + 1e-20)
        self.eastnormaln = zeros((2, self.get("nx") + 2, self.get("ny") + 2))
        for i in range(2):
            self.eastnormaln[i] = eastnormal[i] / (sum(eastnormal**2, 
                axis=0) ** 0.5 + 1e-20)

    def newplot(self, **kwargs):
        """ Creates a figure for 'dump' plots """
        from matplotlib.pyplot import subplots

        self.dumpfig, self.dumpax = subplots(**kwargs)


    def xy(self, 
            x=[], 
            y=[], 
            new=False,
            xlabel='', 
            ylabel='', 
            xlim=(None, None),
            ylim=(None, None),
            iax=0, 
            figsize = (7,5),
            nrows=1,
            ncols=1,
            color='k',
            plottype='plot',
            **kwargs):
        from matplotlib.pyplot import fignum_exists
        try:
            self.dumpfig.number
        except:
            self.newplot(nrows=nrows, ncols=ncols, figsize=figsize)
        if (not fignum_exists(self.dumpfig.number)) or new:
            self.newplot(nrows=nrows, ncols=ncols, figsize=figsize)

        getattr(self.dumpfig.get_axes()[iax], plottype)(
            x, y, color=color, **kwargs
        )
        self.dumpfig.get_axes()[iax].set_xlim(xlim)
        self.dumpfig.get_axes()[iax].set_ylim(ylim)
        self.dumpfig.get_axes()[iax].set_xlabel(xlabel)
        self.dumpfig.get_axes()[iax].set_ylabel(ylabel)

    def savefig(self, fname, **kwargs):
        from matplotlib.pyplot import fignum_exists
        try:
            if not fignum_exists(self.dumpfig.number):
                raise ValueError('No open figure found!')
        except:
            raise ValueError('No open figure found!')
        self.dumpfig.savefig(fname, **kwargs)

          

    def createpolycollection(self, rm, zm, margins=0.05, setparams=True):
        """Creates a poly collection and records boundaries"""
        from matplotlib.collections import PolyCollection
        from numpy import append
        try:
            from uedge import com
        except:
            pass
        from numpy import concatenate

        if setparams is True:
            '''
            try:
                self.nx
            except:
                self.nx = rm.shape[0] - 2
            try:
                self.ny
            except:
                self.ny = rm.shape[1] - 2
            '''
            self.sep = {}
            self.nx = self.get('nx')
            self.ny = self.get('ny')
            if self.snull:
                ixpt1 = self.get("ixpt1")[0]
                ixpt2 = self.get("ixpt2")[0]
                iysptrx = self.get("iysptrx")
                # Initialize separatrix lines to be drawn
                for line in ['ileg', 'oleg', 'core', 'iboundi', 'iboundo',
                    'obound', 'iplate', 'oplate', 'coresep', 'corepatch']:
                    self.sep[line] = {}
                # Define lines to be drawn
                for xy in ['r', 'z']:
                    self.sep['ileg'][xy] = locals()[f"{xy}m"][1:ixpt1+2, iysptrx+1, 1]
                    self.sep['oleg'][xy] = locals()[f"{xy}m"][ixpt2:-1, iysptrx+1, 2]
                    self.sep['core'][xy] = locals()[f"{xy}m"][ixpt1+1:ixpt2+1, 1, 1]
                    self.sep['corepatch'][xy] = locals()[f"{xy}m"][ixpt2, 1, 1:3]
                    self.sep['coresep'][xy] = locals()[f"{xy}m"][ixpt1:ixpt2+1, iysptrx+1, 2]
                    self.sep['obound'][xy] = self.getomit(f"{xy}m")[:, -1, 3]
                    self.sep['iboundi'][xy] = self.getomit(f"{xy}m")[:ixpt1+1, 1, 2]
                    self.sep['iboundo'][xy] = self.getomit(f"{xy}m")[ixpt2+1:, 1, 1]
                    self.sep['iplate'][xy] = self.getomit(f"{xy}m")[1, :, 1]
                    self.sep['oplate'][xy] = self.getomit(f"{xy}m")[-2, :, 2]
                    if self.get(f"{xy}bdry") is not None:
                        if 'efitsep' not in self.sep:
                            self.sep['efitsep'] = {}
                        self.sep['efitsep'][f'{xy}'] = self.get(f"{xy}bdry") 
            elif self.dnull is not False:
                iysptrx1 = self.get("iysptrx1")
                iysptrx2 = self.get("iysptrx2")
                ixpt1 = self.get("ixpt1")
                ixpt2 = self.get("ixpt2")
                ixlb = self.get("ixlb")
                ixrb = self.get("ixrb")
                for line in [   'ilegu', 'olegu', 'ilegl', 'olegl',
                                'corei', 'coreo', 'oboundi', 'oboundo',
                                'iboundli', 'iboundlo', 'iboundui', 'ibounduo',
                                'iplateu', 'oplateu', 'iplatel', 'oplatel',
                                    'isepi', 'isepo', 'osepi', 'osepo']:
                    self.sep[line] = {}
                for xy in ['r', 'z']:
                    self.sep['oplatel'][xy] = locals()[f"{xy}m"][
                        ixlb[1], :, 2
                    ]
                    self.sep['oplateu'][xy] = locals()[f"{xy}m"][
                        ixrb[1]+1, :, 1
                    ]
                    self.sep['iplatel'][xy] = locals()[f"{xy}m"][
                        ixlb[0], :, 2
                    ]
                    self.sep['iplateu'][xy] = locals()[f"{xy}m"][
                        ixrb[0]+1, :, 1
                    ]
                    # RADIAL OUTER BOUNDS
                    self.sep['iboundui'][xy] = locals()[f"{xy}m"][
                        ixpt2[0]+1:ixrb[0]+1, 1, 1
                    ]
                    self.sep['ibounduo'][xy] = locals()[f"{xy}m"][
                        ixlb[1]:ixpt1[1]+1, 1, 2
                    ]
                    self.sep['iboundli'][xy] = locals()[f"{xy}m"][
                        ixlb[0]:ixpt1[0]+1, 1, 2
                    ]
                    self.sep['iboundlo'][xy] = locals()[f"{xy}m"][
                        ixpt2[1]+1:ixrb[1]+2, 1, 1
                    ]
                    self.sep['oboundi'][xy] = locals()[f"{xy}m"][
                        ixlb[0]:ixrb[0]+2, -1, 1
                    ]
                    self.sep['oboundo'][xy] = locals()[f"{xy}m"][
                        ixlb[1]:ixrb[1]+2, -1, 1
                    ]
                    # CORE BOUNDS
                    corei = locals()[f"{xy}m"][
                        ixpt1[0]+1:ixpt2[0]+1, 1, 1
                    ]
                    self.sep['corei'][xy] = append(corei, locals()[f"{xy}m"][
                        ixpt2[0], 1, 2]
                    )
                    coreo = locals()[f"{xy}m"][\
                            ixpt1[1]+1:ixpt2[1]+1, 1, 1
                    ]
                    self.sep['coreo'][xy] = append(coreo, locals()[f"{xy}m"][\
                            ixpt2[1], 1, 2]
                    )
                    # OUTER SEPARATRIX
                    self.sep['osepo'][xy]  = locals()[f"{xy}m"][
                        ixlb[1]:ixpt2[1]+1, iysptrx1[0]+1, 2
                    ]
                    self.sep['osepi'][xy]  = locals()[f"{xy}m"][
                        ixpt1[0]+1:ixrb[0]+2, iysptrx1[0]+1, 1
                    ]
                    isepi = locals()[f"{xy}m"][
                        ixpt1[0]+1:ixpt2[0]+1, iysptrx2[0]+1, 1
                    ]
                    self.sep['isepi'][xy] = append(isepi, locals()[f"{xy}m"][
                        ixpt2[0], iysptrx2[0]+1, 2]
                    )
                    isepo = locals()[f"{xy}m"][\
                            ixpt1[1]+1:ixpt2[1]+1, iysptrx2[0]+1, 1
                    ]
                    self.sep['isepo'][xy] = append(isepo, locals()[f"{xy}m"][\
                            ixpt2[1], iysptrx2[0]+1, 2]
                    )
            if self.dnull=="upper":
                for xy in ['r', 'z']:
                    self.sep['ilegu'][xy] = locals()[f"{xy}m"][\
                            ixpt2[0]+1:ixrb[0]+2, iysptrx2[0]+1, 1
                    ]
                    self.sep['olegu'][xy] = \
                        locals()[f"{xy}m"][
                            ixlb[1]:ixpt1[1]+1, iysptrx2[0]+1, 2
                    ]
                    self.sep['olegl'][xy] = locals()[f"{xy}m"][\
                            ixpt2[1]+1:ixrb[1]+2, iysptrx1[0]+1, 1
                    ]
                    self.sep['ilegl'][xy] = \
                        locals()[f"{xy}m"][
                            ixlb[0]:ixpt1[0]+1, iysptrx1[0]+1, 2
                    ]
                
            elif self.dnull in ["lower", "balanced"]:
                for xy in ['r', 'z']:
                    # OUTER SEPARATRIX
                    self.sep['osepo'][xy]  = locals()[f"{xy}m"][
                        ixlb[1]:ixpt2[1]+1, iysptrx2[0]+1, 2
                    ]
                    self.sep['osepi'][xy]  = locals()[f"{xy}m"][
                        ixpt1[0]+1:ixrb[0]+2, iysptrx2[0]+1, 1
                    ]
                    isepi = locals()[f"{xy}m"][
                        ixpt1[0]+1:ixpt2[0]+1, iysptrx1[0]+1, 1
                    ]
                    self.sep['isepi'][xy] = append(isepi, locals()[f"{xy}m"][
                        ixpt2[0], iysptrx1[0]+1, 2]
                    )
                    isepo = locals()[f"{xy}m"][\
                            ixpt1[1]+1:ixpt2[1]+1, iysptrx1[0]+1, 1
                    ]
                    self.sep['isepo'][xy] = append(isepo, locals()[f"{xy}m"][\
                            ixpt2[1], iysptrx1[0]+1, 2]
                    )


                    self.sep['ilegu'][xy] = locals()[f"{xy}m"][\
                            :ixrb[0]+2, iysptrx2[0]+1, 1
                    ]
                    self.sep['olegu'][xy] = \
                        locals()[f"{xy}m"][
                            ixpt2[1]:, iysptrx2[0]+1, 2
                    ]
                    self.sep['olegl'][xy] = locals()[f"{xy}m"][\
                            ixpt2[1]:, iysptrx1[0]+1, 2
                    ]
                    self.sep['ilegl'][xy] = \
                        locals()[f"{xy}m"][
                            :ixpt1[0]+1, iysptrx1[0]+1, 2
                    ]
    
            if self.dnull=="balanced":
                del(self.sep['isepo'])
                del(self.sep['isepi'])

        vertices = []
        # Loop through all cells, omitting guard cells
        for i in range(1, len(rm) - 1):
            for j in range(1, len(rm[i]) - 1):
                vert = []
                for k in [1, 2, 4, 3]:
                    vert.append([rm[i, j, k], zm[i, j, k]])
                vertices.append(vert)
        self.nodes = vertices
        return PolyCollection(vertices)

    def checkusn(self, array, flip=True):
        if flip is False:
            return array
        elif self.get("geometry")[0].strip().lower().decode("UTF-8") == "uppersn":
            return -array + self.disp
        else:
            return array

    def profile(
        self,
        x,
        y,
        ax=None,
        xlim=(None, None),
        ylim=(0, None),
        figsize=(7, 5),
        xlabel=None,
        ylabel=None,
        title=None,
        logx=False,
        logy=False,
        color="k",
        watermark=True,
        **kwargs
    ):
        """Plots y as function of x"""
        from matplotlib.pyplot import figure, Axes, Figure

        if ax is None:  # Create figure if no axis supplied
            f = figure(title, figsize=figsize)
            ax = f.add_subplot()
        elif isinstance(ax, Figure):
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

    def mesh(
        self,
        z=None,
        rm=None,
        zm=None,
        ax=None,
        linewidth=0.05,
        linecolor="k",
        aspect="equal",
        figsize=(5, 7),
        cmap="magma",
        units="",
        xlim=(None, None),
        ylim=(None, None),
        zrange=(None, None),
        log=False,
        vessel=True,
        plates=True,
        platecolor=None,
        lcfs=True,
        lcfscolor='grey',
        title=None,
        grid=False,
        flip=True,
        watermark=True,
        mask=None,
        colorbar=False,
        interactive=False,
        **kwargs,
    ):
        """General plotting function
        z - values, if any. If None, plots grid
        rm, zm - radial and horizontal nodes
        """
        from matplotlib.pyplot import figure, Figure
        from h5py import File
        from matplotlib.colors import LogNorm
        from matplotlib.collections import PolyCollection
        from copy import deepcopy
        from numpy import array
        try:
            from uedge import com, bbb, grd
        except:
            pass
        try:
            self.vertices
        except:
            self.createvertices(self.get("rm"), self.get("zm"))
    
        if z is not None:
            # Assume this is a path to a grid file
            if isinstance(z, str):
                # Read data into local vars to be used to create vertices later
                with File(z) as f:
                    rm = self.getomit(f['grid/com/rm'][()])
                    zm = self.getomit(f['grid/com/zm'][()])
                    z = None
            else:
                if len(z.shape) != 2:
                    raise ValueError('Array to be plotted must be two-dimensional!')    

        if ax is None:
            f = figure(title, figsize=figsize)
            ax = f.add_subplot()
        if isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        if (rm is None) or (zm is None):
            # Use stored PolyCollection
            if (
                self.get("geometry")[0].strip().lower().decode("UTF-8") == "uppersn"
            ) and (flip is True):
                vertices = deepcopy(self.uppersnvertices)
            else:
                vertices = deepcopy(self.vertices)
        else:  # Create collection from data
            vertices = self.createpolycollection(rm, zm)
        if grid is False:
            vertices.set_linewidths(1)
            vertices.set_edgecolors("face")
        else:
            vertices.set_edgecolors("lightgrey")
            vertices.set_linewidths(0.08)
        if z is None:  # Plot grid
            vertices.set_facecolor((0, 0, 0, 0))
            vertices.set_edgecolors(linecolor)
            vertices.set_linewidths(linewidth)
        else:
            vertices.set_cmap(cmap)
            vertices.set_array(z[1:-1, 1:-1].reshape(self.nx * self.ny))
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
            self.lcfs(ax, flip, color=lcfscolor)
        if vessel is True:
            self.vessel(ax, flip)
        if plates is True:
            self.plates(ax, flip, color=platecolor)
        ax.autoscale_view()
        ax.set_title(title)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel("R [m]")
        ax.set_ylabel("Z [m]")
        ax.set_aspect(aspect)
        if (z is not None) and (colorbar is True):
            cbar = ax.get_figure().colorbar(vertices, ax=ax)
            cbar.ax.set_ylabel(units, va="bottom")

        if watermark is True:
            self.watermark(ax.get_figure(), bottom=0.1, left=0.02, right=0.95)

        self.Qvertices = vertices
        if interactive is True:
            return cbar, vertices
        else:
            return ax.get_figure()

    def lcfs(self, ax, flip=True, color="grey", linewidth=0.5, **kwargs):
        """Plots LCFS on ax"""
        try:
            from uedge import com, bbb, grd
        except:
            pass
        from numpy import int64
        for key, coords in self.sep.items():
            if not isinstance(coords['r'], (int, int64)):
                ax.plot(
                    coords['r'], 
                    self.checkusn(coords['z'], flip), 
                    color=color,
                    linewidth=linewidth,
                    label="lcfs",
                    **kwargs
                )
        return


        try:
            if not type(rbdry, int64):
                ax.plot(
                    rbdry,
                    self.checkusn(zbdry, flip),
                    color=color,
                    linewidth=linewidth,
                    label='lcfs'
                )
                ax.plot(
                    self.isepr,
                    self.checkusn(self.isepz, flip),
                    color=color,
                    linewidth=linewidth,
                    label='lcfs'
                )
                ax.plot(
                    self.osepr,
                    self.checkusn(self.osepz, flip),
                    color=color,
                    linewidth=linewidth,
                    label='lcfs'
                )
                ax.plot(
                    self.pfrboundr,
                    self.checkusn(self.pfrboundz, flip),
                    color=color,
                    linewidth=linewidth,
                    label='lcfs'
                )
                ax.plot(
                    self.otboundr,
                    self.checkusn(self.otboundz, flip),
                    color=color,
                    linewidth=linewidth,
                    label='lcfs'
                )
                ax.plot(
                    self.itboundr,
                    self.checkusn(self.itboundz, flip),
                    color=color,
                    linewidth=linewidth,
                    label='lcfs'
                )
                ax.plot(
                    self.solboundr,
                    self.checkusn(self.solboundz, flip),
                    color=color,
                    linewidth=linewidth,
                    label='lcfs'
                )
                plotted = True
        except:
            pass
        try:
            if not isinstance(self.get("rbdry"), int64):
                ax.plot(
                    self.get("rbdry"),
                    self.checkusn(self.get("zbdry"), flip),
                    color=color,
                    linewidth=linewidth,
                    label='lcfs'
                )
                ax.plot(
                    self.isepr,
                    self.checkusn(self.isepz, flip),
                    color=color,
                    linewidth=linewidth,
                    label='lcfs'
                )
                ax.plot(
                    self.osepr,
                    self.checkusn(self.osepz, flip),
                    color=color,
                    linewidth=linewidth,
                    label='lcfs'
                )
                ax.plot(
                    self.pfrboundr,
                    self.checkusn(self.pfrboundz, flip),
                    color=color,
                    linewidth=linewidth,
                    label='lcfs'
                )
                ax.plot(
                    self.otboundr,
                    self.checkusn(self.otboundz, flip),
                    color=color,
                    linewidth=linewidth,
                    label='lcfs'
                )
                ax.plot(
                    self.itboundr,
                    self.checkusn(self.itboundz, flip),
                    color=color,
                    linewidth=linewidth,
                    label='lcfs'
                )
                ax.plot(
                    self.solboundr,
                    self.checkusn(self.solboundz, flip),
                    color=color,
                    linewidth=linewidth,
                    label='lcfs'
                )
                plotted = True
        except:
            pass
        if plotted is False:
            ax.plot(
                self.sepcorer, 
                self.checkusn(self.sepcorez, flip),
                color=color,
                linewidth=linewidth,
                label='lcfs'
            )
            ax.plot(
                self.isepr,
                self.checkusn(self.isepz, flip),
                color=color,
                linewidth=linewidth,
                label='lcfs'
            )
            ax.plot(
                self.osepr,
                self.checkusn(self.osepz, flip),
                color=color,
                linewidth=linewidth,
                label='lcfs'
            )
            ax.plot(
                self.pfrboundr,
                self.checkusn(self.pfrboundz, flip),
                color=color,
                linewidth=linewidth,
                label='lcfs'
            )
            ax.plot(
                self.otboundr,
                self.checkusn(self.otboundz, flip),
                color=color,
                linewidth=linewidth,
                label='lcfs'
            )
            ax.plot(
                self.itboundr,
                self.checkusn(self.itboundz, flip),
                color=color,
                linewidth=linewidth,
                label='lcfs'
            )
            ax.plot(
                self.solboundr,
                self.checkusn(self.solboundz, flip),
                color=color,
                linewidth=linewidth,
                label='lcfs'
            )

    def vessel(self, ax, flip=True):
        """Plots vessel on ax"""
        try:
            from uedge import com, bbb, grd
        except:
            pass

        try:
            if self.get("xlim") is not None:
                for params in [['k-', 3], ['y-', 1]]:
                    ax.plot(
                        self.get("xlim"),
                        self.checkusn(self.get("ylim"), flip),
                        params[0],
                        linewidth=params[1],
                        label='vessel',
                    )
        except:
            pass


    def plates(self, ax, flip=True, color=None):
        """Plot plates on ax"""
        try:
            from uedge import com, bbb, grd
        except:
            pass

        if color is None:
            color = ['b', 'r']
        elif isinstance(color, str):
            p1c = color
            p2c = color
        
        try:
            for i in range(2):
                ax.plot(
                    self.get("rplate{}".format(i+1)),
                    self.checkusn(self.get("zplate{}".format(i+1)), flip),
                    "-",
                    color=color[i],
                    linewidth=1.5,
                    label='plate{}'.format(i+1),
                )
        except:
            pass

    def watermark(self, *args, **kwargs):
        pass

    def mesh_masked(self, z, zmask, maskvalues, figsize=(5,7), 
        **kwargs):
        from matplotlib.pyplot import subplots
        f, ax = subplots(figsize=figsize)       


        cbar, vertices = self.mesh(z, interactive=True, ax=ax, **kwargs)

        mask = zmask[1:-1,1:-1].reshape(self.nx*self.ny)
        vertices.set_alpha( [1*( (x<maskvalues[0]) or (x>maskvalues[1])) for x in mask])
        return f


    def quiver(
        self,
        pol, 
        rad,
        orthogonal=False,
        color='k',
        width=3e-3,
        headwidth=3,
        headlength=2,
        flip=True,
        plates=False,
        vessel=False,
        lcfs=True,
        linewidth=0.1,
        alpha=True,
        uniformsize=False,
        xlim=(None, None),
        ylim=(None, None),
        **kwargs
    ):

        f = self.mesh(plates=plates, lcfs=lcfs, vessel=vessel, 
            linewidth=linewidth)
        ax = f.get_axes()[0]
        # Check whether coords are given as poloidal or radial
        if orthogonal is False:
            x = pol * self.eastnormaln[0] + rad * self.northnormaln[0]
            y = pol * self.eastnormaln[1] + rad * self.northnormaln[1]
        else:
            x = pol
            y = rad
        magnitude = (pol**2 + rad**2)**0.5
        if uniformsize is True:
            x /= (magnitude + 1e-20)
            y /= (magnitude + 1e-20)

        magnitude=magnitude[1:-1, 1:-1]
        magnitude /= magnitude.max()
        if alpha is False:
            magnitude = magnitude ** 0        


        rm = self.get("rm")[:, :, 0]
        zm = self.get("zm")[:, :, 0]
        ax.quiver(
            rm[1:-1, 1:-1].ravel(), 
            zm[1:-1, 1:-1].ravel(), 
            x[1:-1,1:-1].ravel(), 
            y[1:-1, 1:-1].ravel(),
            pivot='mid',
            color=color,
            alpha=magnitude,
            width=width,
            headwidth=headwidth,
            headlength=headlength,
            headaxislength=headlength,
            **kwargs
        )
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel("R [m]")
        ax.set_ylabel("Z [m]")

        return f

        
    def streamline(
        self,
        pol,
        rad,
        resolution=(500j, 800j),
        linewidth="magnitude",
        broken_streamlines=False,
        color="k",
        maxlength=0.4,
        mask=True,
        density=2,
        xlim=(None, None),
        ylim=(None, None),
        **kwargs
    ):
        from numpy import zeros, sum, transpose, mgrid, nan, array, cross, nan_to_num
        from scipy.interpolate import griddata, bisplrep
        from matplotlib.patches import Polygon
        from copy import deepcopy

        rm = self.get("rm")
        zm = self.get("zm")
        nx = self.get("nx")
        ny = self.get("ny")
        # Create polygons for masking
        outerx = []
        outerx = outerx + list(rm[::-1][-self.get("ixpt1")[0] :, 0, 2])
        outerx = outerx + list(rm[0, :, 1])
        outerx = outerx + list(rm[:, -1, 3])
        outerx = outerx + list(rm[:, ::-1][-1, :, 4])
        outerx = outerx + list(rm[::-1][: nx - self.get("ixpt2")[0], 0, 1])
        outery = []
        outery = outery + list(zm[::-1][-self.get("ixpt1")[0] :, 0, 2])
        outery = outery + list(zm[0, :, 1])
        outery = outery + list(zm[:, -1, 3])
        outery = outery + list(zm[:, ::-1][-1, :, 4])
        outery = outery + list(zm[::-1][: nx - self.get("ixpt2")[0], 0, 1])

        innerx = rm[self.get("ixpt1")[0] + 1 : self.get("ixpt2")[0] + 1, 0, 1]
        innery = zm[self.get("ixpt1")[0] + 1 : self.get("ixpt2")[0] + 1, 0, 1]

        outer = Polygon(
            array([outerx, outery]).transpose(),
            closed=True,
            facecolor="white",
            edgecolor="none",
        )
        inner = Polygon(
            array([innerx, innery]).transpose(),
            closed=True,
            facecolor="white",
            edgecolor="none",
        )
        x = pol * self.eastnormaln[0] + rad * self.northnormaln[0]
        y = pol * self.eastnormaln[1] + rad * self.northnormaln[1]

        gx, gy = mgrid[
            rm.min() : rm.max() : resolution[0], zm.min() : zm.max() : resolution[1]
        ]

        xinterp = griddata(
            (self.sxmid[1, 0, 1:-1, 1:-1].ravel(), 
            self.sxmid[1, 1, 1:-1, 1:-1].ravel()), x[1:-1, 1:-1].ravel(),
            (gx, gy),
        )
        yinterp = griddata(
            (self.symid[1, 0, 1:-1, 1:-1].ravel(), 
            self.symid[1, 1, 1:-1, 1:-1].ravel()), y[1:-1, 1:-1].ravel(),
            (gx, gy),
        )

        if mask is True:
            for i in range(gx.shape[0]):
                for j in range(gx.shape[1]):
                    p = (gx[i, j], gy[i, j])
                    if (inner.contains_point(p)) or (not outer.contains_point(p)):
                        xinterp[i, j] = nan
                        yinterp[i, j] = nan

        f = self.mesh()
        if linewidth == "magnitude":
            linewidth = (xinterp**2 + yinterp**2) ** 0.5
            linewidth = linewidth.transpose()
            maxwidth = nan_to_num(deepcopy(linewidth)).max()
            linewidth /= maxwidth

        f.get_axes()[0].streamplot(
            gx.transpose(),
            gy.transpose(),
            xinterp.transpose(),
            yinterp.transpose(),
            linewidth=linewidth,
            broken_streamlines=broken_streamlines,
            color=color,
            maxlength=maxlength,
            density=density,
            **kwargs
        )

        ax = f.get_axes()[0]
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel("R [m]")
        ax.set_ylabel("Z [m]")
        return f



    def contour(
        self,
        var,
        resolution=(500j, 800j),
        color="k",
        flip=True,
        levels=14,
        linewidth=0.5,
        interplcfs=1,
        mask=True,
        ax=None,
        labels=True,
        vessel=False,
        plates=False,
        latecolor='k',
        lcfs=True,
        gridlinewidth=0.01,
        watermark=False,
        lcfscolor='grey',
        gridlinecolor='k',
        plotgrid=False,
        xlim=(None,None),
        ylim=(None,None),
        method='linear',
        **kwargs
    ):
        from numpy import zeros, sum, transpose, mgrid, nan, array, cross, nan_to_num, concatenate
        from scipy.interpolate import griddata, bisplrep
        from matplotlib.patches import Polygon
        from matplotlib.pyplot import subplots, Figure
        from copy import deepcopy

        if ax is None:
            f = self.mesh(linewidth=gridlinewidth, vessel=vessel, plates=plates,
                flip=flip, lcfs=lcfs, lcfscolor=lcfscolor, linecolor=gridlinecolor)
            ax = f.get_axes()[0]
        elif isinstance(ax, Figure):
            f = ax
            ax = f.get_axes()[0]
            
        else:
            f = ax.get_figure()
        if plotgrid is True:
            self.mesh(linewidth=gridlinewidth, vessel=vessel, plates=plates,
                flip=flip, lcfs=lcfs, lcfscolor=lcfscolor, linecolor=gridlinecolor, ax=ax, watermark=watermark)

        rm = self.get("rm")
        zm = self.get("zm")
        nx = self.get("nx")
        ny = self.get("ny")
        ixlb = self.get("ixlb")
        ixrb = self.get("ixrb")
        ixpt1 = self.get("ixpt1")
        ixpt2 = self.get("ixpt2")
        iysptrx1 = self.get("iysptrx1")
        iysptrx2 = self.get("iysptrx2")
        if self.get("geometry")[0].strip().lower().decode("UTF-8") == "uppersn":
            if flip is True:
                zm = -zm + self.disp

        # Create polygons for masking
        outerx = []
        outery = []
        innerx = []
        innery = []
        if self.snull:
            outerx = outerx + list(rm[::-1][-self.get("ixpt1")[0] :, 0, 2])
            outerx = outerx + list(rm[0, :, 1])
            outerx = outerx + list(rm[:, -1, 3])
            outerx = outerx + list(rm[:, ::-1][-1, :, 4])
            outerx = outerx + list(rm[::-1][: nx - self.get("ixpt2")[0], 0, 1])

            outery = outery + list(zm[::-1][-self.get("ixpt1")[0] :, 0, 2])
            outery = outery + list(zm[0, :, 1])
            outery = outery + list(zm[:, -1, 3])
            outery = outery + list(zm[:, ::-1][-1, :, 4])
            outery = outery + list(zm[::-1][: nx - self.get("ixpt2")[0], 0, 1])

            innerx = innerx + list(rm[self.get("ixpt1")[0] + 1 : self.get("ixpt2")[0] + 1, 0, 1])
            innery = innery + list(zm[self.get("ixpt1")[0] + 1 : self.get("ixpt2")[0] + 1, 0, 1])
        else:
            outerx = outerx + list(rm[::-1][-ixpt1[0]-1:, 0, 2])
            outerx = outerx + list(rm[0, :, 1])
            outerx = outerx + list(rm[:ixrb[0]+1, -1, 3])
            outerx = outerx + list(rm[:, ::-1][ixrb[0]+1, :, 4])
            outerx = outerx + list(rm[ixpt2[0]+1:ixrb[0]+1, 0, 1][::-1])

            outerx = outerx + list(rm[ixlb[1]+1:ixpt1[1]+1, 0, 1][::-1])
            outerx = outerx + list(rm[ixlb[1],:, 1])
            outerx = outerx + list(rm[ixlb[1]:, -1, 3])
            outerx = outerx + list(rm[-1, :, 4][::-1])
            outerx = outerx + list(rm[ixpt2[1]+1:, 0, 1][::-1])

            outery = outery + list(zm[::-1][-ixpt1[0]-1:, 0, 2])
            outery = outery + list(zm[0, :, 1])
            outery = outery + list(zm[:ixrb[0]+1, -1, 3])
            outery = outery + list(zm[:, ::-1][ixrb[0]+1, :, 4])
            outery = outery + list(zm[ixpt2[0]+1:ixrb[0]+1, 0, 1][::-1])

            outery = outery + list(zm[ixlb[1]+1:ixpt1[1]+1, 0, 1][::-1])
            outery = outery + list(zm[ixlb[1],:, 1])
            outery = outery + list(zm[ixlb[1]:, -1, 3])
            outery = outery + list(zm[-1, :, 4][::-1])
            outery = outery + list(zm[ixpt2[1]+1:, 0, 1][::-1])

            innerx = innerx + list(rm[ixpt1[0] + 1 : ixpt2[0] + 1, 0, 1])
            innerx = innerx + list(rm[ixpt1[1] + 1 : ixpt2[1] + 1, 0, 1])
            innery = innery + list(zm[ixpt1[0] + 1 : ixpt2[0] + 1, 0, 1])
            innery = innery + list(zm[ixpt1[1] + 1 : ixpt2[1] + 1, 0, 1])

        outer = Polygon(
            array([outerx, outery]).transpose(),
            closed=True,
            facecolor="white",
            edgecolor="none",
        )
        inner = Polygon(
            array([innerx, innery]).transpose(),
            closed=True,
            facecolor="white",
            edgecolor="none",
        )


        gx, gy = mgrid[
            rm.min() : rm.max() : resolution[0], zm.min() : zm.max() : resolution[1]
        ]

        
        if interplcfs == 1:
            # Fix for LCFS - no interpolation knows about the LCFS
            for i in range(len(ixpt1)):
                xtrax = 0.5*(rm[ixpt2[i], iysptrx2[i], 4] + rm[ixpt2[i], iysptrx2[i], 2])
                xtray = 0.5*(zm[ixpt2[i], iysptrx2[i], 4] + zm[ixpt2[i], iysptrx2[i], 2])
                d0 = (  (rm[ixpt2[i], iysptrx2[i], 0] - xtrax)**2 + 
                        (zm[ixpt2[i], iysptrx2[i], 0] - xtray)**2)**0.5
                d1 = (  (rm[ixpt1[i], iysptrx1[i], 0] - xtrax)**2 + 
                        (zm[ixpt1[i], iysptrx1[i], 0] - xtray)**2)**0.5
                varp1 = var[ixpt1[i], iysptrx1[i]]
                xtraz = (d0*var[ixpt2[i], iysptrx2[i]] + d1*varp1)/(d0 + d1)
                orig = (concatenate((rm[:, :, 0].ravel(), array((xtrax,)))),
                        concatenate((zm[:, :, 0].ravel(), array((xtray,)))))
                fullvar = concatenate((var.ravel(), array((xtraz,))))

        elif interplcfs == 2:
            raise NotImplementedError("Option interplcfs=2 not implemented/verified")
            # Fix for LCFS - no interpolation knows about the LCFS
            xtrax, xtray, xtraz = [], [], []
            for i in range(ixpt1+1, ixpt2+1):
                xtrax.append(rm[i, iysptrx, 4])
                xtray.append(zm[i, iysptrx, 4])
                d0 = (  (rm[i, iysptrx, 0] - xtrax[-1])**2 + 
                        (zm[i, iysptrx, 0] - xtray[-1])**2)**0.5
                if i == ixpt2:
                    d1 = (  (rm[ixpt1, iysptrx, 0] - xtrax[-1])**2 + 
                            (zm[ixpt1, iysptrx, 0] - xtray[-1])**2)**0.5
                    varp1 = var[ixpt1, iysptrx]
                else:
                    d1 = (  (rm[i+1, iysptrx, 0] - xtrax[-1])**2 + 
                            (zm[i+1, iysptrx, 0] - xtray[-1])**2)**0.5
                    varp1 = var[i+1, iysptrx]
                xtraz.append( (d0*var[i, iysptrx] + d1*varp1)/\
                        (d0 + d1))
            orig = (concatenate((rm[:, :, 0].ravel(), array(xtrax))),
                    concatenate((zm[:, :, 0].ravel(), array(xtray))))
            fullvar = concatenate((var.ravel(), array(xtraz)))

        else:
            orig = (rm[:, :, 0].ravel(), zm[:, :, 0].ravel())
            fullvar = var.ravel()
            

        # Previous implementation, which (incorrectly) assumed values to
        # be given for cell-centers
        #        xinterp = griddata( (rm[1:-1,1:-1,0].ravel(), zm[1:-1,1:-1,0].ravel()),
        #            x[1:-1,1:-1].ravel(), (gx, gy))
        #        yinterp = griddata( (rm[1:-1,1:-1,0].ravel(), zm[1:-1,1:-1,0].ravel()),
        #            y[1:-1,1:-1].ravel(), (gx, gy))
        varinterp = griddata( orig,
            fullvar, (gx, gy), method=method
#            (rm[1:-1, 1:-1, 0].ravel(), zm[1:-1, 1:-1, 0].ravel()), 
#            var[1:-1, 1:-1].ravel(), (gx, gy)
        )

        if mask is True:
            for i in range(gx.shape[0]):
                for j in range(gx.shape[1]):
                    p = (gx[i, j], gy[i, j])
                    if (inner.contains_point(p)) or (not outer.contains_point(p)):
                        varinterp[i, j] = nan

        CS = ax.contour(
                gx, gy, varinterp, levels=levels, linewidths=linewidth, 
                colors=color, **kwargs
        )
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel("R [m]")
        ax.set_ylabel("Z [m]")
        if labels is True:
            ax.clabel(CS, fontsize=9, inline=False, fmt='% 1.2e')

        return f


    def LFSleg_polprof(self, var, ax=None, irad=3, wrad =2, figsize=(7,5), color='r', 
        ylim=None, ylabel='', alpha=0.05, plottype='plot', **kwargs):
        """ Plots the poloidal profile of var between X-point and target """
        from matplotlib.pyplot import subplots, Figure
        from numpy import amin, amax

        if ax is None:  # Create figure if no axis supplied
            f, ax = subplots(figsize=figsize)
        elif ax is Figure:
            ax = ax.get_axes()[0]


        ixpt2 = self.get('ixpt2')[0]
        xcs = self.get('xcs')
        iysptrx = self.get('iysptrx')
        ny = self.get('ny')
        # Poloidal distance of X-point
        xcsXpt = 0.5*sum(xcs[ixpt2:ixpt2+2])
        
        # Shade the region that spans wrad flux-tubes to each side of irad
        ax.fill_between(
            xcs[ixpt2:-1] - xcsXpt, 
            amin( var[ixpt2:-1, iysptrx+irad-wrad : iysptrx + irad+wrad], axis=1),
            amax( var[ixpt2:-1, iysptrx+irad-wrad : iysptrx + irad+wrad], axis=1),
            color = color,
            alpha = alpha
        )

        # Plot irad
        getattr(ax, plottype)( 
                xcs[ixpt2:-1] - xcsXpt, 
                var[ixpt2:-1, iysptrx+irad], 
                color=color, **kwargs
        )

        ax.set_xlim(0, 1.04*(xcs[-2] - xcsXpt))
        if ylim is None:
            ax.autoscale(axis='y')
        else:
            ax.set_ylim(ylim)
        ax.set_ylabel(ylabel)
        ax.set_xlabel('Poloidal distance from X-point along LFS leg')


        return ax.get_figure()

































