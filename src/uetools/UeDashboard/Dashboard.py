
class Dashboard():
    def dashboard(self, flip=False, **kwargs):
        from matplotlib.pyplot import ion, ioff, subplots
        from matplotlib.widgets import Slider, RangeSlider, Button
        from matplotlib import is_interactive
        from copy import deepcopy

        # ioff()

        z = self.get("ne")
        f, ax = subplots(figsize=(15, 8))

        try:
            self.vertices
        except:
            self.createvertices(self.get("rm"), self.get("zm"))
        if (self.get("geometry")[0].strip().lower().decode("UTF-8") == "uppersn") and (
            flip is True
        ):
            vertices = deepcopy(self.uppersnvertices)
        else:
            vertices = deepcopy(self.vertices)

        vertices.set_cmap("magma")

        self.plotlcfs(ax, flip)
        self.plotvessel(ax, flip)
        self.plotplates(ax, flip)

        ax.set_aspect("equal")
        ax.autoscale_view()
        ax.set_xlabel("R [m]")
        ax.set_ylabel("Z [m]")
        ax.add_collection(vertices)
        vertices.set_array(z[1:-1, 1:-1].reshape(self.nx * self.ny))
        vertices.set_edgecolors("face")
        cbar = f.colorbar(vertices, ax=ax)
        cbar.ax.set_ylabel(r"$\rm n_e$ [$\rm m^{-3}$]", va="bottom")

        f.axes[0].set_position([0.06, 0.1, 0.27, 0.85])
        f.axes[1].set_position([0.35, 0.1, 0.82, 0.85])

        zrange_position = f.add_axes([0.4, 0.1, 0.02, 0.85])
        zrange_slider = RangeSlider(
            zrange_position,
            "",
            z[1:-1, 1:-1].min(),
            z[1:-1, 1:-1].max(),
            valinit=((z[1:-1, 1:-1].min(), z[1:-1, 1:-1].max())),
            orientation="vertical",
        )

        te_position = f.add_axes([0.5, 0.9, 0.12, 0.03])
        te_button = Button(te_position, "Electron temperature")

        #        vertices.set_array(z[1:-1,1:-1].reshape(self.nx*self.ny))
        #        vertices.set_clim(*zrange)
        #        if log is True:
        #            vertices.set_norm(LogNorm())
        #            vertices.set_clim(*zrange)

        # TODO: devise scheme to look for variables in memory, from
        try:
            kwargs["zrange"]
            origrange = kwargs["zrange"]
        except:
            kwargs["zrange"] = (z[1:-1, 1:-1].min(), z[1:-1, 1:-1].max())
            origrange = kwargs["zrange"]
        #        cbar, verts = self.plotmesh(z, ax=ax, watermark=False,
        #            interactive=True, **kwargs)
        #        reset_position = f.add_axes([0.02, 0.02, 0.1, 0.04])
        #        reset_button = Button(reset_position, 'Reset')
        #            , color='gold',hovercolor='skyblue')

        def reset_zrange(event):
            zrange_slider.reset()

        def update(val):
            from numpy import floor, ceil, linspace

            vertices.set_clim(zrange_slider.val)

        def plot_te(event):
            var = self.get("te")[1:-1, 1:-1].reshape(self.nx * self.ny) / 1.602e-19
            zrange_slider.set_val((var.min(), var.max()))
            zrange_slider.set_min(var.min())
            zrange_slider.set_max(var.max())
            vertices.set_clim(var.min(), var.max())
            vertices.set_array(var)

        #        reset_button.on_clicked(reset_zrange)
        zrange_slider.on_changed(update)
        te_button.on_clicked(plot_te)
        f.show()

        ion()

        return f, cbar, vertices, (zrange_slider, te_button)

    def plot_ne(self, vertices):
        z = self.get("ne")[1:-1, 1:-1].reshape(self.nx * self.ny)
        vertices.set_edgecolors("face")
        vertices.set_facecolors((0, 0, 0, 0))
        vertices.set_array(z)
        vertices.set_clim(z.min(), z.max())

    def plot_grid(self, vertices):
        vertices.set_facecolors("none")
        vertices.set_linewidths(0.1)
        vertices.set_edgecolors("k")
        return
