

class  DB_2DPlots:

    def ng_2Dseries(self, species, **kwargs):
        return self.plot_2Dseries(self.get("ng")[:,:,:,species], **kwargs)

    def tg_2Dseries(self, species, **kwargs):
        return self.plot_2Dseries(self.get("tg")[:,:,:,species] / 1.602e-19, **kwargs)

    def ni_2Dseries(self, species, **kwargs):
        return self.plot_2Dseries(self.get("ni")[:,:,:,species], **kwargs)

    def ti_2Dseries(self, **kwargs):
        return self.plot_2Dseries(self.get("ti") / 1.602e-19, **kwargs)

    def ne_2Dseries(self, **kwargs):
        return self.plot_2Dseries(self.get("ne"), **kwargs)

    def te_2Dseries(self, **kwargs):
        return self.plot_2Dseries(self.get("te") / 1.602e-19, **kwargs)

    def plot_2Dseries(self, 
        vararray, 
        cmap='magma', 
        grid=False,
        linewidth=0.05, 
        linecolor='k',
        flip=True, 
        lcfscolor='grey', 
        lcfs=True,
        platecolor='r', 
        plates=True, 
        vessel=True,
        **kwargs):
        """Returns a series of figures to scroll through"""
        from matplotlib.pyplot import subplots, ion, ioff
        from matplotlib.widgets import Slider, RangeSlider
        from copy import deepcopy

        # Return to avoid garbage collection!
        return  InteractivePlot(
                vararray,
                self,
                watermark=False,
                interactive=True, 
                cmap=cmap,
                grid=grid,
                linewidth=linewidth,
                linecolor=linecolor,
                flip=flip,
                lcfscolor=lcfscolor,
                lcfs=lcfs,
                platecolor=platecolor,
                plates=plates,
                vessel=vessel,
                **kwargs

        )


class InteractivePlot():
    def __init__(self, vararray, db, flip=True, xlim=None, ylim=None, **kwargs):
        """Returns a series of figures to scroll through"""
        from matplotlib.pyplot import subplots, ion, ioff
        from matplotlib.widgets import Slider, RangeSlider
        from copy import deepcopy
        
        self.vararray = vararray
        self.db = db
        for key, value in kwargs.items():
            self.__setattr__(key, value)
        self.flip=flip

        ioff()
        self.f, self.ax = subplots(figsize=(7, 8))

        try:
            kwargs["zrange"]
            origrange = kwargs["zrange"]
        except:
            kwargs["zrange"] = (
                vararray[:,1:-1, 1:-1].min(),
                vararray[:,1:-1, 1:-1].max(),
            )
            origrange = kwargs["zrange"]

        self.cbar, self.verts = db.getcase(0).plot.mesh(
                vararray[0], 
                flip=self.flip,
                ax=self.ax, 
                **kwargs
        )
        self.f.axes[0].set_position([0.125, 0.13, 0.55, 0.85])
        self.f.axes[1].set_position([0.7, 0.13, 0.82, 0.85])
        slice_position = self.f.add_axes([0.1, 0.02, 0.65, 0.04])
        self.slice_slider = Slider(
            slice_position,
            self.db.sortvar,
            self.db.sortvalues.min(),
            self.db.sortvalues.max(),
            valstep=self.db.sortvalues,
        )
        zrange_position = self.f.add_axes([0.85, 0.13, 0.04, 0.85])
        self.zrange_slider = RangeSlider(
            zrange_position,
            "",
            vararray[:,1:-1, 1:-1].min(),
            vararray[:,1:-1, 1:-1].max(),
            valinit=(origrange),
            orientation="vertical",
        )
        if xlim is None:
            xlim = [self.db.get('rm').min(), self.db.get('rm').max()]
            for var in ['xlim', 'rplate1', 'rplate2']:
                xlim[0] = min(self.db.get(var, ravel=True).min(), xlim[0])-0.03
                xlim[1] = max(self.db.get(var, ravel=True).max(), xlim[1])+0.03
        if ylim is None:
            ylim = [self.db.get('zm').min(), self.db.get('zm').max()]
            for var in ['ylim', 'zplate1', 'zplate2']:
                ylim[0] = min(self.db.get(var, ravel=True).min(), ylim[0])-0.03
                ylim[1] = max(self.db.get(var, ravel=True).max(), ylim[1])+0.03
        if (
            self.db.getcase(0).get("geometry")[0].strip().lower().decode("UTF-8") == "uppersn"
        ) and (self.flip is True):
            ylim = [db.getcase(0).disp - x for x in ylim[::-1]]
        # TODO:
        # Normalize view to the X-point location?


        self.ax.set_xlim((xlim[0], xlim[1]))
        self.ax.set_ylim((ylim[0], ylim[1]))

        self.lcfs = False
        self.plates = False
        self.vessel = False
        def update(val):
            from numpy import where, array

            if self.slice_slider.val != self.slce:
                self.slce = self.slice_slider.val
                index = where(self.db.sortvalues == self.slce)[0][0]
                c = self.db.getcase(index)
                nodes = array(c.plot.nodes)
                if not self.flip:
                    nodes[:,:,1] = c.disp-nodes[:,:,1]
                self.verts.set_verts(nodes)
                self.verts.set_array(
                    vararray[index, 1:-1, 1:-1].reshape(c.get("nx") * c.get("ny"))
                )
                for line in self.ax.lines:
                    if line.get_label() == 'lcfs':
                        self.lcfs = True
                    elif 'plate' in line.get_label():
                        self.plates = True
                    elif line.get_label() == 'vessel':
                        self.vessel = True
                    line.remove()
                if self.lcfs:
                    c.plot.lcfs(self.ax, flip=self.flip, color="grey", linewidth=0.5)
                if self.vessel:
                    c.plot.vessel(self.ax, flip=self.flip)
                if self.plates:
                    c.plot.plates(self.ax, flip=self.flip)

            self.verts.set_clim(self.zrange_slider.val)
            self.verts.set_cmap(self.cmap)
            return
        self.slce = self.slice_slider.val
        self.slice_slider.on_changed(update)
        self.zrange_slider.on_changed(update)
        

        self.f.show()
        ion()


