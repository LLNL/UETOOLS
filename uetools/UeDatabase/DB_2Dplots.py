

class  DB_2DPlots:

    def ng_2Dseries(self, species, **kwargs):
        return self.plot_2Dseries(self.get("ng")[:,:,:,species], **kwargs)

    def tg_2Dseries(self, species, **kwargs):
        self.plot_2Dseries(self.get("tg")[:,:,:,species] / 1.602e-19, **kwargs)

    def ni_2Dseries(self, species, **kwargs):
        return self.plot_2Dseries(self.get("ni")[:,:,:,species], **kwargs)

    def ti_2Dseries(self, **kwargs):
        self.plot_2Dseries(self.get("ti") / 1.602e-19, **kwargs)

    def ne_2Dseries(self, **kwargs):
        return self.plot_2Dseries(self.get("ne"), **kwargs)

    def te_2Dseries(self, **kwargs):
        self.plot_2Dseries(self.get("te") / 1.602e-19, **kwargs)

    def plot_2Dseries(self, vararray, **kwargs):
        """Returns a series of figures to scroll through"""
        from matplotlib.pyplot import subplots, ion, ioff
        from matplotlib.widgets import Slider, RangeSlider

        ioff()
        f, ax = subplots(figsize=(7, 8))

        try:
            kwargs["zrange"]
            origrange = kwargs["zrange"]
        except:
            kwargs["zrange"] = (
                vararray[1:-1, 1:-1, :].min(),
                vararray[1:-1, 1:-1, :].max(),
            )
            origrange = kwargs["zrange"]

        c = self.getcase(0)
        cbar, verts = c.plotmesh(
            vararray[0], ax=ax, watermark=False, interactive=True, **kwargs
        )
        f.axes[0].set_position([0.125, 0.13, 0.55, 0.85])
        f.axes[1].set_position([0.7, 0.13, 0.82, 0.85])
        slice_position = f.add_axes([0.1, 0.02, 0.65, 0.04])
        slice_slider = Slider(
            slice_position,
            self.sortvar,
            self.sortvalues.min(),
            self.sortvalues.max(),
            valstep=self.sortvalues,
        )
        zrange_position = f.add_axes([0.85, 0.13, 0.04, 0.85])
        zrange_slider = RangeSlider(
            zrange_position,
            "",
            vararray[1:-1, 1:-1, :].min(),
            vararray[1:-1, 1:-1, :].max(),
            valinit=(origrange),
            orientation="vertical",
        )

        def update(val):
            from numpy import where

            slce = slice_slider.val
            index = where(self.sortvalues == slce)[0][0]
            verts.set_array(
                vararray[index, 1:-1, 1:-1].reshape(
                    self.getcase(index).get("nx") * self.getcase(index).get("ny")
                )
            )
            verts.set_clim(zrange_slider.val)

        slice_slider.on_changed(update)
        zrange_slider.on_changed(update)

        f.show()
        ion()
        return f, slice_slider, zrange_slider


