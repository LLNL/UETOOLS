# Object for plotting
from matplotlib.pyplot import ion

ion()

# TODO: implement divergence plotting/calculation

class Plot:
    def __init__(self, *args, rm=None, zm=None, **kwargs):
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

        self.createvertices(rm, zm)
        super().__init__(*args, **kwargs)
        return

    def createvertices(self, rm, zm):
        from numpy import zeros, transpose, cross, sum
        # CREATE POLYGON COLLECTIONS TO USE
        self.vertices = self.createpolycollection(rm, zm)
        if self.get("geometry")[0].strip().lower().decode("UTF-8") == "uppersn":
            self.disp = 0
            if self.get("rmagx") + self.get("zmagx") == 0:
                self.disp = -(-zm).min()
            else:
                self.disp = 2 * self.get("zmagx")
            self.uppersnvertices = self.createpolycollection(
                rm, -zm + self.disp, setparams=False
            )

        # CREATE NORMAL VECTOR IN LOCAL CELL COORDINATES
        nodes = zeros((self.get("nx") + 2, self.get("ny") + 2, 5, 2))
        nodes[:, :, :, 0] = self.get("rm")
        nodes[:, :, :, 1] = self.get("zm")
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


    def plot(self, 
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
        try:
            from uedge import com
        except:
            pass
        from numpy import concatenate

        if setparams is True:
            try:
                self.nx
            except:
                self.nx = rm.shape[0] - 2
            try:
                self.ny
            except:
                self.ny = rm.shape[1] - 2
            ixpt1 = self.get("ixpt1")[0]
            ixpt2 = self.get("ixpt2")[0]
            iysptrx = self.get("iysptrx")
            # TODO: where to pass/find xpoint indices?
            self.isepr = rm[1 : ixpt1 + 2, iysptrx + 1, 1]
            self.isepz = zm[1 : ixpt1 + 2, iysptrx + 1, 1]
            self.osepr = rm[ixpt2 : -1, iysptrx + 1, 2]
            self.osepz = zm[ixpt2 : -1, iysptrx + 1, 2]

            self.sepcorer = self.get("rm")[ixpt1 : ixpt2 + 1, iysptrx + 1, 2]
            self.sepcorez = self.get("zm")[ixpt1 : ixpt2 + 1, iysptrx + 1, 2]
            self.solboundr = self.get("rm")[:, -2, 3:]
            self.solboundz = self.get("zm")[:, -2, 3:]
            self.pfrboundr = concatenate((
                                self.get("rm")[: ixpt1+1, 1, 2], 
                                self.get("rm")[ixpt2+1 :, 1, 1]
                            ))
            self.pfrboundz = concatenate((
                                self.get("zm")[: ixpt1+1, 1, 2],
                                self.get("zm")[ixpt2+1 :, 1, 1] 
                            ))
            self.itboundr = self.get("rm")[1, :, 1]
            self.itboundz = self.get("zm")[1, :, 1]
            self.otboundr = self.get("rm")[-2, :, 2]
            self.otboundz = self.get("zm")[-2, :, 2]

        vertices = []
        # Loop through all cells, omitting guard cells
        for i in range(1, len(rm) - 1):
            for j in range(1, len(rm[i]) - 1):
                vert = []
                for k in [1, 2, 4, 3]:
                    vert.append([rm[i, j, k], zm[i, j, k]])
                vertices.append(vert)
        return PolyCollection(vertices)

    def checkusn(self, array, flip=False):
        if flip is False:
            return array
        elif self.get("geometry")[0].strip().lower().decode("UTF-8") == "uppersn":
            return -array + self.disp

    def plotprofile(
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

    def plotmesh(
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
        flip=False,
        watermark=True,
        mask=None,
        colorbar=True,
        interactive=False,
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
                    rm = f['grid/com/rm'][()]
                    zm = f['grid/com/zm'][()]
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
            vertices.set_edgecolors(linecolor)
            vertices.set_linewidths(linewidth)
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
            self.plotlcfs(ax, flip, color=lcfscolor)
        if vessel is True:
            self.plotvessel(ax, flip)
        if plates is True:
            self.plotplates(ax, flip, color=platecolor)
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

        if interactive is True:
            return cbar, vertices
        else:
            return ax.get_figure()

    def plotlcfs(self, ax, flip=False, color="grey", linewidth=0.5, **kwargs):
        """Plots LCFS on ax"""
        try:
            from uedge import com, bbb, grd
        except:
            pass
        from numpy import int64
        plotted = False
        try:
            if not type(com.rbdry, int64):
                ax.plot(
                    com.rbdry,
                    self.checkusn(com.zbdry, flip),
                    color=color,
                    linewidth=linewidth,
                )
                ax.plot(
                    self.isepr,
                    self.checkusn(self.isepz, flip),
                    color=color,
                    linewidth=linewidth,
                )
                ax.plot(
                    self.osepr,
                    self.checkusn(self.osepz, flip),
                    color=color,
                    linewidth=linewidth,
                )
                ax.plot(
                    self.pfrboundr,
                    self.checkusn(self.pfrboundz, flip),
                    color=color,
                    linewidth=linewidth,
                )
                ax.plot(
                    self.otboundr,
                    self.checkusn(self.otboundz, flip),
                    color=color,
                    linewidth=linewidth,
                )
                ax.plot(
                    self.itboundr,
                    self.checkusn(self.itboundz, flip),
                    color=color,
                    linewidth=linewidth,
                )
                ax.plot(
                    self.solboundr,
                    self.checkusn(self.solboundz, flip),
                    color=color,
                    linewidth=linewidth,
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
                )
                ax.plot(
                    self.isepr,
                    self.checkusn(self.isepz, flip),
                    color=color,
                    linewidth=linewidth,
                )
                ax.plot(
                    self.osepr,
                    self.checkusn(self.osepz, flip),
                    color=color,
                    linewidth=linewidth,
                )
                ax.plot(
                    self.pfrboundr,
                    self.checkusn(self.pfrboundz, flip),
                    color=color,
                    linewidth=linewidth,
                )
                ax.plot(
                    self.otboundr,
                    self.checkusn(self.otboundz, flip),
                    color=color,
                    linewidth=linewidth,
                )
                ax.plot(
                    self.itboundr,
                    self.checkusn(self.itboundz, flip),
                    color=color,
                    linewidth=linewidth,
                )
                ax.plot(
                    self.solboundr,
                    self.checkusn(self.solboundz, flip),
                    color=color,
                    linewidth=linewidth,
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
            )
            ax.plot(
                self.isepr,
                self.checkusn(self.isepz, flip),
                color=color,
                linewidth=linewidth,
            )
            ax.plot(
                self.osepr,
                self.checkusn(self.osepz, flip),
                color=color,
                linewidth=linewidth,
            )
            ax.plot(
                self.pfrboundr,
                self.checkusn(self.pfrboundz, flip),
                color=color,
                linewidth=linewidth,
            )
            ax.plot(
                self.otboundr,
                self.checkusn(self.otboundz, flip),
                color=color,
                linewidth=linewidth,
            )
            ax.plot(
                self.itboundr,
                self.checkusn(self.itboundz, flip),
                color=color,
                linewidth=linewidth,
            )
            ax.plot(
                self.solboundr,
                self.checkusn(self.solboundz, flip),
                color=color,
                linewidth=linewidth,
            )

    def plotvessel(self, ax, flip=False):
        """Plots vessel on ax"""
        try:
            from uedge import com, bbb, grd
        except:
            pass

        try:
            ax.plot(com.xlim, self.checkusn(com.ylim, flip), "k-", linewidth=3)
            ax.plot(com.xlim, self.checkusn(com.ylim, flip), "y-", linewidth=1)
        except:
            pass
        try:
            if self.get("xlim") is not None:
                ax.plot(
                    self.get("xlim"),
                    self.checkusn(self.get("ylim"), flip),
                    "k-",
                    linewidth=3,
                )
                ax.plot(
                    self.get("xlim"),
                    self.checkusn(self.get("ylim"), flip),
                    "y-",
                    linewidth=1,
                )
        except:
            pass


    def plotplates(self, ax, flip=False, color=None):
        """Plot plates on ax"""
        try:
            from uedge import com, bbb, grd
        except:
            pass

        try:
            ax.plot(grd.rplate1, self.checkusn(grd.zplate1, flip), "b-", linewidth=1.5)
            ax.plot(grd.rplate2, self.checkusn(grd.zplate2, flip), "r-", linewidth=1.5)
        except:
            pass
        if color is None:
            p1c = 'b'
            p2c = 'r'
        else:
            p1c = color
            p2c = color
        
        try:
            if self.get("rplate1") is not None:
                ax.plot(
                    self.get("rplate1"),
                    self.checkusn(self.get("zplate1"), flip),
                    "-",
                    color=p1c,
                    linewidth=1.5,
                )
                ax.plot(
                    self.get("rplate2"),
                    self.checkusn(self.get("zplate2"), flip),
                    "-",
                    color=p2c,
                    linewidth=1.5,
                )
        except:
            pass

    def watermark(self, figure, bottom=0.15, top=0.95, left=0.09, right=0.98):
        """Adds metadata to figure"""
        from time import ctime

        label = '{}, case "{}"\n'.format(ctime(), self.casename)
        label += 'UEDGE {} v{}, UETOOLS v{}, user "{}", hostname "{}"\n'.format(
            self.uedge_ver.replace("$", "\$"),
            self.pyver,
            self.uetoolsversion,
            self.user,
            self.hostname,
        )
        try:
            label += 'cwd "{}"'.format(self.location)
        except:
            label += 'cwd "{}"'.format(self.casefname)
        figure.subplots_adjust(bottom=bottom, top=top, left=left, right=right)
        figure.text(0.995, 0.005, label, fontsize=4, horizontalalignment="right")

        return

    def plotmesh_masked(self, z, zmask, maskvalues, figsize=(5,7), 
        **kwargs):
        from matplotlib.pyplot import subplots
        f, ax = subplots(figsize=figsize)       


        cbar, vertices = self.plotmesh(z, interactive=True, ax=ax, **kwargs)

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
        flip=False,
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

        f = self.plotmesh(plates=plates, lcfs=lcfs, vessel=vessel, 
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

        f = self.plotmesh()
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
        flip=False,
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
            f = self.plotmesh(linewidth=gridlinewidth, vessel=vessel, plates=plates,
                flip=flip, lcfs=lcfs, lcfscolor=lcfscolor, linecolor=gridlinecolor)
            ax = f.get_axes()[0]
        elif isinstance(ax, Figure):
            f = ax
            ax = f.get_axes()[0]
            
        else:
            f = ax.get_figure()
        if plotgrid is True:
            self.plotmesh(linewidth=gridlinewidth, vessel=vessel, plates=plates,
                flip=flip, lcfs=lcfs, lcfscolor=lcfscolor, linecolor=gridlinecolor, ax=ax)

        rm = self.get("rm")
        zm = self.get("zm")
        nx = self.get("nx")
        ny = self.get("ny")
        ixpt1 = self.get("ixpt1")[0]
        ixpt2 = self.get("ixpt2")[0]
        iysptrx = self.get("iysptrx")
        if self.get("geometry")[0].strip().lower().decode("UTF-8") == "uppersn":
            if flip is True:
                zm = -zm + self.disp

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

        gx, gy = mgrid[
            rm.min() : rm.max() : resolution[0], zm.min() : zm.max() : resolution[1]
        ]

        if interplcfs == 1:
            # Fix for LCFS - no interpolation knows about the LCFS
            xtrax = 0.5*(rm[ixpt2, iysptrx, 4] + rm[ixpt2, iysptrx, 2])
            xtray = 0.5*(zm[ixpt2, iysptrx, 4] + zm[ixpt2, iysptrx, 2])
            d0 = (  (rm[ixpt2, iysptrx, 0] - xtrax)**2 + 
                    (zm[ixpt2, iysptrx, 0] - xtray)**2)**0.5
            d1 = (  (rm[ixpt1, iysptrx, 0] - xtrax)**2 + 
                    (zm[ixpt1, iysptrx, 0] - xtray)**2)**0.5
            varp1 = var[ixpt1, iysptrx]
            xtraz = (d0*var[ixpt2, iysptrx] + d1*varp1)/(d0 + d1)
            orig = (concatenate((rm[:, :, 0].ravel(), array((xtrax,)))),
                    concatenate((zm[:, :, 0].ravel(), array((xtray,)))))
            fullvar = concatenate((var.ravel(), array((xtraz,))))

        elif interplcfs == 2:
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

































