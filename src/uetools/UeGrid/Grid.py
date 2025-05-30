class Grid:
    """Object for interacting with the UEDGE GRIDUE grid generator

    This is still a stub, currently only offering some limited
    plotting routines

    Inherits GridPlot

    Attributes
    ----------
    plot: GridPlot object, providing plotting functions

    Methods
    -------
    pick_aeqdskdata(geqdsk, ncontour=250, interpres=2000, **kwargs)
        Lets the user define the magnetic axes and X-point locations
    """

    def __init__(self, case):
        """Links class to uetools.Case functions"""
        self.get = case.getue
        self.setue = case.setue
        self.plot = GridPlot(case)

    def pick_aeqdskdata(self, geqdsk, ncontour=250, interpres=2000, 
        colormesh=False, **kwargs):
        """Tool for manually defining aeqdsk data

        Passes **kwargs to self.plot.efit

        Arguments
        ---------
        geqdsk - path to geqdsk file to be read

        Keyword arguments
        -----------------
        ncontour : int (default = 250)
            number of contours to plot for EFIT psis
        interpres : int (default = 2000)
            resolution of interpolation used to identify magnetic
            axis and X-points
        colormesh : bool (default = False)
            switch whether to overlay equilibrium on colored mesh

        Return
        ------
        None

        Modifies
        --------
        Nothing (yet)
        """
        from uedge import com
        from scipy.interpolate import griddata
        from scipy.optimize import fmin
        from numpy import linspace, array, gradient, sum, meshgrid, where
        from matplotlib.pyplot import subplots, ginput, waitforbuttonpress
        from shapely import Polygon, Point, LineString

        f, ax = subplots(figsize=(7, 9))
        ax.plot([], [], "ko")
        ax.plot([], [], "ro")
        ax.plot([], [], "bo")
        self.plot.efit(geqdsk, ax=ax, ncontour=ncontour, **kwargs)

        rorig, zorig = self.get("nxefit"), self.get("nyefit")
        xf = lambda r: linspace(0, self.get("xdim"), r) + self.get("rgrid1")
        yf = lambda r: linspace(0, self.get("zdim"), r) - (
            self.get("zdim") * 0.5 - self.get("zmid")
        )
        xo, yo = meshgrid(xf(rorig), yf(zorig))
        fold = self.get("fold").transpose()
        grad = gradient(fold, xo[0, 1] - xo[0, 0], yo[1, 0] - yo[0, 0])

        xi, yi = meshgrid(xf(interpres), yf(interpres))
        interp = griddata((xo.ravel(), yo.ravel()), fold.ravel(), (xi, yi))
        gradinterp = griddata(
            (xo.ravel(), yo.ravel()),
            array((grad[0] ** 2 + grad[1] ** 2) ** 0.5).ravel(),
            (xi, yi),
        )

        if colormesh is True:
            c = ax.pcolormesh(xi, yi, gradinterp, cmap='hot_r', vmax=0.1)#, vmin=-0.2, vmax=0.2)

        print("Manually identify the following points in order")
        print("(Choose by clicking in the figure)")  # , undo by right-clicking)')
        textbox = f.get_axes()[0].text(0.87, 1.5, "", backgroundcolor="w", zorder=10)
        boxtext = "Please click on the:\n{}\nPress return to accepts surface\nMouse-click to clear and pick again"
        xlim = self.get("xlim")
        ylim = self.get("ylim")
        vessel = Polygon(zip(xlim, ylim))

        pts = []
        gradmin = []
        psimin = []
        contours = []
        crosses = []
        strike_points = []
        ptlabels = [" - Magnetic axis", " - Lower X-point", " - Upper X-point"]
        colors = ["b", "r", "m"]
        # TODO: Fix message
        # TODO: Figure out why double-clicks are needed to clear+redo
        for i in range(3):
            print(ptlabels[i])
            while True:
                textbox.set_text(boxtext.format(ptlabels[i]))
                pt = ginput(1, 0, False)
                if len(pt) > 0:
                    nearest = [
                        abs(xf(interpres) - pt[0][0]).argmin(),
                        abs(yf(interpres) - pt[0][1]).argmin(),
                    ]
                    gradmin_nearest = gradinterp[nearest[1], nearest[0]]
                    psimin_nearest = interp[nearest[1], nearest[0]]
                    gradmin_optimum = gradinterp[
                        nearest[1]
                        - int(interpres / 50) : nearest[1]
                        + int(interpres / 50),
                        nearest[0]
                        - int(interpres / 50) : nearest[0]
                        + int(interpres / 50),
                    ].min()
                    y, x = where(gradinterp == gradmin_optimum)
                    psimin_optimum = interp[y[0], x[0]]
                    if psimin_optimum < psimin_nearest:
                        gradmin_use = gradmin_optimum
                        psimin_use = psimin_optimum
                        pts_use = [xf(interpres)[x[0]], yf(interpres)[y[0]]]
                    else:
                        gradmin_use = gradmin_nearest
                        psimin_use = psimin_nearest
                        pts_use = [xf(interpres)[nearest[0]], yf(interpres)[nearest[1]]]
                    crosses.append(
                        ax.plot(*pts_use, "x", color=colors[i], markersize=12)
                    )
                    contours.append(
                        ax.contour(
                            xf(rorig),
                            yf(zorig),
                            fold,
                            [psimin_use],
                            colors=colors[i],
                            linewidths=1.5,
                            linestyles="-",
                        )
                    )
                    # Identify strike-point locations
                    if i > 0:
                        sps = []
                        # Store line (x,y)'s
                        for line in contours[-1].collections:
                            lines = []
                            for v in line.get_paths():
                                lines.append([v.vertices[:, 0], v.vertices[:, 1]])
                        # Find contour intersects w/ vessel
                        intersects = []
                        for segment in lines:
                            # Create LineString from object
                            seg = LineString(zip(segment[0], segment[1]))
                            if vessel.exterior.intersects(seg):
                                for point in vessel.exterior.intersection(seg).geoms:
                                    intersects.append(point)
                        dists = [Point(pts_use).distance(x) for x in intersects]
                        strike_points.append(
                            [x for _, x in sorted(zip(dists, intersects))][:2]
                        )
                        for point in strike_points[-1]:
                            sps.append(ax.plot(*point.xy, "o", color=colors[i]))
                    textbox.set_text(boxtext.format(ptlabels[i]))
                    f.canvas.draw()
                    # Check whether line is accepted or rejected
                    if waitforbuttonpress():
                        try:
                            textbox.set_text(boxtext.format(ptlabels[i + 1]))
                        except:
                            textbox.set_text("All points defined!")
                        # Store lines for finding strike-points

                        # TODO: Find and store intersects here
                        f.canvas.draw()
                        gradmin.append(gradmin_use)
                        psimin.append(psimin_use)
                        pts.append(pts_use)
                        break
                    else:
                        for line in contours[i].collections:
                            line.remove()
                        del contours[i]
                        l = crosses[i].pop(0)
                        l.remove()
                        del crosses[i]
                        if i > 0:
                            for sp in sps:
                                o = sp.pop(0)
                                o.remove()
                        f.canvas.draw()

        print(psimin)
        # TODO: Set UEDGE variables based on detected variable
        # ( rseps, zseps, rseps, zseps2,) rvsin, rvsout, zvsin, zvsout


class GridPlot:
    """Class providing grid generation plotting routines

    Methods
    -------
    flx(ax=None, surfaces=None)
        Plots the flux-surfaces constructed by UEDGE
    efit(geqdsk, aeqdsk=None, ax=None, ncontour=80, color='grey',
            linestyle='solid', sepcolor='k', linewidth=0.5)
        Plots the efit equilibrium
    """
    def __init__(self, case):
        self.get = case.get
        self.setue = case.setue
        self.reload = case.reload

    def flx(self, ax=None, surfaces=None):
        """Plots flux surfaces from UEDGE memory

        Based on plotflx.bas by Rensink, Rognlien & Porter

        Keyword arguments
        -----------------
        ax : matplotlib.pyplot.Figure or Axes (default = None)
            axis to plot on. If None, creates new figure
        surfaces : range (default = None)
            range of surfaces to plot. If none, plots all

        Returns:
        ========
        matplotlib.pyplot.Figure
        """
        from Forthon import packageobject
        from matplotlib.pyplot import subplots, Figure, Axes

        # Execute flxrun
        packageobject("flx").__getattribute__("flxrun")
        # Validate ax
        if ax is None:
            f, ax = subplots(figsize=(7, 9))
        elif isinstance(ax, Axes):
            pass
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        else:
            raise TypeError("Axes type {} not compatible".format(type(ax)))
        # Get the number of surfaces available
        mxs = (
            2 * (self.get("nycore")[0] + self.get("nysol")[0] + self.get("nyout")[0])
            + 2
        )
        # Validate surface plot request
        if surfaces is None:
            rng = range(mxs + 1)
        elif not isinstance(surfaces, range):
            raise TypeError("surfaces must be type 'range'")
        else:
            rng = surfaces
            if rng.stop > mxs:
                rng = range(rng.start, mxs)
            if rng.start < 0:
                rng = range(0, rng.stop)
        # Plot vessel if present in EQDSK files
        if self.get("nlim") > 0:
            ax.plot(self.get("xlim"), self.get("ylim")+self.get('zshift'), "k-", linewidth=3)
            ax.plot(
                self.get("xlim"), self.get("ylim")+self.get('zshift'), "-", linewidth=1.5, color="yellow"
            )
        # Plot target plates, if specified
        try:
            ax.plot(self.get("rplate1"), self.get("zplate1"), "ro-")
        except:
            pass
        try:
            ax.plot(self.get("rplate2"), self.get("zplate2"), "b.-")
        except:
            pass
        # Helpers
        jmin = self.get("jmin")
        jmax = self.get("jmax")
        jsptrx = self.get("jsptrx")
        xcurve = self.get("xcurve")
        ycurve = self.get("ycurve")
        ijumpf = self.get("ijumpf")
        # Plot the flux surfaces within the sepcified range
        for i in rng:
            # Plot SOL flux surfaces for each half-mesh
            if ((i >= jmin[0] - 1) and (i < jsptrx[0])) or (
                (i >= jsptrx[1]) and (i <= jmax[1])
            ):
                ax.plot(
                    xcurve[:, i][abs(xcurve[:, i]) > 0],
                    ycurve[:, i][abs(ycurve[:, i]) > 0],
                    "k-",
                    linewidth=0.3,
                )
            # Plot CORE/PFR flux surfaces for each half-mesh
            elif ((i >= jsptrx[0]) and (i <= jmax[0])) or (
                (i >= jmin[1] - 1) and (i <= jsptrx[1] + 1)
            ):
                ax.plot(
                    xcurve[:, i][abs(xcurve[:, i]) > 0][: ijumpf[i]],
                    ycurve[:, i][abs(ycurve[:, i]) > 0][: ijumpf[i]],
                    "k-",
                    linewidth=0.3,
                )
                ax.plot(
                    xcurve[:, i][abs(xcurve[:, i]) > 0][ijumpf[i] :],
                    ycurve[:, i][abs(ycurve[:, i]) > 0][ijumpf[i] :],
                    "k-",
                    linewidth=0.3,
                )
            for i in [jsptrx[0] - 1, jsptrx[1] - 1]:
                ax.plot(
                    xcurve[:, i][abs(xcurve[:, i]) > 0],
                    ycurve[:, i][abs(ycurve[:, i]) > 0],
                    "r-",
                    linewidth=0.5,
                )
        ax.set_aspect("equal")
        ax.set_xlabel("Horizontal position [m]")
        ax.set_ylabel("Vertical position [m]")
        return ax.get_figure()

    def efit(
        self,
        geqdsk,
        aeqdsk=None,
        ax=None,
        ncontour=80,
        color="grey",
        linestyle="solid",
        sepcolor="k",
        linewidth=0.5,
        labels=False,
    ):
        """Function to plot EFIT contours

        Arguments
        ---------
        geqdsk - path to EFIT file to read

        Keyword arguments
        -----------------
        aeqdsk : str (default = None)
            Path to aeqdsk file. If None, does not read any aeqdsk
        ax : matplotlib.pyplot.Figure or Axes (default = None)
            axis to plot on. If None, creates new figure
        ncontour : int (default = 80)
            number of contours to plot for EFIT psis
        color : str (default = 'grey')
            color of plot contours
        linestyle : str (default = 'solid')
            line style for plot contours
        sepcolor : str (default = 'k')
            color of separatrices
        linewidth : float (default = 0.5)
            line width of contours

        Returns
        -------
        matplotlib.pyplot.Figure
        """
        from matplotlib.pyplot import subplots, Figure, Axes
        from copy import deepcopy
        from numpy import linspace
        from uedge import com, flx
        from os.path import exists
        from scipy.interpolate import interp2d
        from Forthon import packageobject

        # Backup original pointers
        oldaeqdskfname = deepcopy(self.get("aeqdskfname"))
        oldgeqdskfname = deepcopy(self.get("geqdskfname"))
        # Set new file paths
        self.setue("geqdskfname", geqdsk)

        # Check whether the aeqdsk file can be located: if not, do not execute aeqdsk()
        if aeqdsk is not None:
            if exists(aeqdsk):
                self.setue("aeqdskfname", aeqdsk)
                packageobject("flx").__getattribute__("aeqdsk")()

        if exists(geqdsk):
            packageobject("flx").__getattribute__("neqdsk")()
        else:
            raise FileNotFoundError(
                'EFIT geqdsk file "{}" not found.\nAborting...'.format(geqdsk)
            )
        self.reload()

        if ax is None:
            f, ax = subplots(figsize=(7, 9))
        elif isinstance(ax, Axes):
            pass
        elif isinstance(ax, Figure):
            ax = ax.get_axes()[0]
        else:
            raise TypeError("Axes type {} not compatible".format(type(ax)))

        # Reconstruct EFIT grid
        x = linspace(0, self.get("xdim"), self.get("nxefit")) + self.get("rgrid1")
        y = linspace(0, self.get("zdim"), self.get("nyefit")) - (
            self.get("zdim") * 0.5 - self.get("zmid")
        )

        fold = self.get("fold").transpose()
        # Interpolate on EFIT grid to find X-points
        interp = interp2d(x, y, fold)

        CS = ax.contour(
            x,
            y,
            fold,
            ncontour,
            colors=color,
            linewidths=linewidth,
            linestyles=linestyle,
        )
        if labels:
            ax.clabel(CS, CS.levels, inline=1, fontsize=10)

        rseps2 = self.get("rseps2")
        zseps2 = self.get("zseps2")
        # Check whether the upper X-point exists
        if (x.min() <= rseps2 <= x.max()) and (y.min() <= zseps2 <= y.max()):
            upperxpoint = interp(rseps2, zseps2)
            ax.contour(
                x,
                y,
                fold,
                [upperxpoint],
                colors=sepcolor,
                linewidths=1,
                linestyles="solid",
            )

        rseps = self.get("rseps")
        zseps = self.get("zseps")
        # Check whether the lower X-point exists
        if (x.min() <= rseps <= x.max()) and (y.min() <= zseps <= y.max()):
            lowerxpoint = interp(rseps, zseps)
            ax.contour(
                x,
                y,
                fold,
                [lowerxpoint],
                colors=sepcolor,
                linewidths=1,
                linestyles="solid",
            )

        ax.plot(self.get("xlim"), self.get("ylim"), "k-", linewidth=2)
        ax.set_aspect("equal")
        ax.set_xlabel("Horizontal position [m]")
        ax.set_ylabel("Vertical position [m]")
        ax.set_title(self.get("runid")[0].decode("UTF-8").strip())

        # Restore original pointers
        self.setue("aeqdskfname", oldaeqdskfname)
        self.setue("geqdskfname", oldgeqdskfname)

        return ax.get_figure()
