

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
        flip=False, 
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
    def __init__(self, vararray, db, **kwargs):
        """Returns a series of figures to scroll through"""
        from matplotlib.pyplot import subplots, ion, ioff
        from matplotlib.widgets import Slider, RangeSlider
        from copy import deepcopy
        
        self.vararray = vararray
        self.db = db
        for key, value in kwargs.items():
            self.__setattr__(key, value)
        
        ioff()
        self.f, self.ax = subplots(figsize=(7, 8))

        try:
            kwargs["zrange"]
            origrange = kwargs["zrange"]
        except:
            kwargs["zrange"] = (
                vararray[1:-1, 1:-1, :].min(),
                vararray[1:-1, 1:-1, :].max(),
            )
            origrange = kwargs["zrange"]

        self.cbar, self.verts = db.getcase(0).plotmesh(
                vararray[0], 
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
            vararray[1:-1, 1:-1, :].min(),
            vararray[1:-1, 1:-1, :].max(),
            valinit=(origrange),
            orientation="vertical",
        )
        xlim = [self.db.get('rm').min(), self.db.get('rm').max()]
        ylim = [self.db.get('zm').min(), self.db.get('zm').max()]
        for var in ['xlim', 'rplate1', 'rplate2']:
            xlim[0] = min(self.db.get(var).min(), xlim[0])
            xlim[1] = max(self.db.get(var).max(), xlim[1])
        for var in ['ylim', 'zplate1', 'zplate2']:
            ylim[0] = min(self.db.get(var).min(), ylim[0])
            ylim[1] = max(self.db.get(var).max(), ylim[1])
        if self.flip is True:
            ylim = [db.getcase(0).disp - x for x in ylim[::-1]]


        self.ax.set_xlim((xlim[0]-0.03, xlim[1]+0.03))
        self.ax.set_ylim((ylim[0]-0.03, ylim[1]+0.03))

        def update(val):
            from numpy import where

            slce = self.slice_slider.val
            index = where(self.db.sortvalues == slce)[0][0]
            buffcase = self.db.getcase(index)
            self.verts.remove()
            if (
                buffcase.get("geometry")[0].strip().lower().decode("UTF-8") == "uppersn"
            ) and (self.flip is True):
                self.verts = deepcopy(buffcase.uppersnvertices)
            else:
                self.verts = deepcopy(buffcase.vertices)
            self.ax.add_collection(self.verts)
            
            if self.grid is False:
                self.verts.set_linewidths(1)
                self.verts.set_edgecolors("face")
            else:
                self.verts.set_edgecolors(linecolor)
                self.verts.set_linewidths(linewidth)
            self.verts.set_cmap(self.cmap)
            
            self.verts.set_array(
                vararray[index, 1:-1, 1:-1].reshape(
                    buffcase.get("nx") * buffcase.get("ny")
                )
            )

            for line in self.ax.get_lines():
                line.remove()
            if self.lcfs:
                buffcase.plotlcfs(self.ax, self.flip, 'grey')#flip, color=lcfscolor)
            if self.plates:
                buffcase.plotplates(self.ax, self.flip, 'r')#flip, color=platecolor)
            if self.vessel:
                buffcase.plotvessel(self.ax, self.flip)#flip)
            self.verts.set_clim(self.zrange_slider.val)

        self.slice_slider.on_changed(update)
        self.zrange_slider.on_changed(update)

        self.f.show()
        ion()


