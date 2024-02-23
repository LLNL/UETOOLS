#import matplotlib
#matplotlib.use("Qt5Agg")

class MainMenu():
    def __init__(self, casepath=None):
        from matplotlib.pyplot import figure
        self.f = figure(figsize=(4, 3))
        if casepath is None:
            self.f.suptitle("No UEDGE case loaded")
                
        

class CaseDashboard2D():
    def __init__(self, 
        case, 
        **kwargs):
        """Returns a series of figures to scroll through"""
        from matplotlib.pyplot import subplots, ion, ioff, figure
        from matplotlib.widgets import Slider, RangeSlider, Button, \
            TextBox, RadioButtons
        from copy import deepcopy
        from os.path import expanduser
        
        self.case = case
        self.var = self.case.get('te')/1.602e-19
        self.varname = "Te [eV]"
        self.case.set_speciesarrays()
        self.multispecies = 0
        self.s = 0
        self.ns = 0
        self.smax = self.case.get("nisp")
        self.nsmax = self.case.get("ngsp")
        self.fromccmap = False
        default_cmaps = ('magma', 'bwr', 'gist_heat', 'jet', 'gist_rainbow','')
        for key, value in kwargs.items():
            self.__setattr__(key, value)
        ioff()
        self.board = figure(layout='constrained', figsize=(14,8))
        [self.f, self.controls] = self.board.subfigures(1,2, wspace=0,
                                                width_ratios=[1,1])
        self.ax = self.f.subplots(1,1)
        [self.buttonarea, self.controlarea] = self.controls.subfigures(2, 1, 
                                                hspace=0, height_ratios=[8,2])

        [self._varbuttonarea, self.radioarea, self._switcharea] = \
                self.buttonarea.subfigures(
                            1, 
                            3, 
                            width_ratios=[1,1,1], 
                            wspace=0
                )
        try:
            kwargs["zrange"]
            origrange = kwargs["zrange"]
        except:
            kwargs["zrange"] = (
                self.var[1:-1, 1:-1].min(),
                self.var[1:-1, 1:-1].max(),
            )
            origrange = kwargs["zrange"]
        self.cbar, self.verts = case.plotmesh(
                self.var, 
                ax=self.ax, 
                interactive=True,
                watermark=False,
                **kwargs
        )
        # TODO: Dropdown with parameters?
        # TODO: Button to flip?
        self.f.get_axes()[0].set_title(self.varname)
        # Set figure location
        self.f.get_axes()[0].set_position([0.05, 0.08, 0.6, 0.88])
        # Set colorbar location
        self.f.get_axes()[1].set_position([0.72, 0.08, 0.06, 0.88])
        # Set slider location
        zrange_position = self.f.add_axes([0.88, 0.1, 0.06, 0.84])
        self.zrange_slider = RangeSlider(
            zrange_position,
            "",
            self.var[1:-1, 1:-1].min(),
            self.var[1:-1, 1:-1].max(),
            valinit=(origrange),
            orientation="vertical",
        )
        self.zrange_slider.valtext.set_visible(False)
        axulim = self.f.add_axes([0.86, 0.96, 0.11, 0.025])
        self.tbulim = TextBox(axulim, "", textalignment='center')
        self.tbulim.on_submit(self.update_upper)
        self.display_value(self.tbulim, kwargs['zrange'][1])

        axllim = self.f.add_axes([0.86, 0.02, 0.11, 0.025])
        self.tbllim = TextBox(axllim, "", textalignment='center')
        self.tbllim.on_submit(self.update_lower)
        self.display_value(self.tbllim, kwargs['zrange'][0])
        maxbuttons = 18
        """ Set up variable button options """
        self.varbuttons = {
            "Electron temperature": {'func': self.te},
            "Ion temperature": {'func': self.ti},
            "Electron density": {'func': self.ne},
            "Ion density": {'func': self.ni},
            "Gas density": {'func': self.ng},
        }
        [self.varbuttonarea, a] = self._varbuttonarea.subfigures(
                    2, 
                    1, 
                    height_ratios=[
                            len(self.varbuttons)/maxbuttons,
                            1-(len(self.varbuttons)/maxbuttons)
                    ]
        )
        buttonaxes = self.varbuttonarea.subplots(\
                    max(len(self.varbuttons), 0), 
                    1
        )
        i = 0
        for key, setup in self.varbuttons.items():
            self.varbuttons[key]['button'] = Button(buttonaxes[i], key, 
                color=(key != "Electron temperature")*'light'+'grey')
            self.varbuttons[key]['button'].on_clicked(setup['func'])
            i += 1
        for j in range(i, len(buttonaxes)):
            buttonaxes[j].set_visible(False)
        self.varbuttons['Electron temperature']['button'].color = 'grey'
        """ Set up the switches """
        self.switches = {
            "Log": {'func': self.toggle_log, 'default': False},
            "Vessel": {'func': self.toggle_vessel, 'default': True},
            "Plates": {'func': self.toggle_plates, 'default': True},
            "Separatrix": {'func': self.toggle_sep, 'default': True},
            "Grid": {'func': self.toggle_grid, 'default': False},
        }
        [self.switcharea, b] = self._switcharea.subfigures(
                    2, 
                    1,
                    height_ratios=[
                            len(self.switches)/maxbuttons,
                            1-(len(self.switches)/maxbuttons)
                    ]
        )
        switchaxes = self.switcharea.subplots(max(len(self.switches), 0), 1)
        i = 0
        for key, setup in self.switches.items():
            if setup['default']:
                color = 'grey'
                setattr(self, key, True)
            else:
                color = 'lightgrey'
                setattr(self, key, False)
            self.switches[key]['button'] = Button(switchaxes[i], 
                        key, 
                        color=color
            )
            self.switches[key]['button'].on_clicked(self.switches[key]['func'])
            i += 1
        for j in range(i, len(switchaxes)):
            switchaxes[j].set_visible(False)
        """ Set up radio buttons """
        plupps = self.nsmax + self.smax + len(default_cmaps)
        radioaxes = self.radioarea.subplots(
                3, 1, 
                height_ratios=[
                    self.smax/plupps, 
                    self.nsmax/plupps, 
                    len(default_cmaps)/plupps
        ]) 
        self.radios = {
            'ion': {
                'labels': tuple(self.case.ionarray), 
                'func': self.update_ionspecies
            },
            'gas': {
                'labels': tuple(self.case.gasarray), 
                'func': self.update_gasspecies
            },
            'cmap': {
                'labels': default_cmaps,
                'func': self.update_cmap
            },
        }
        i = 0
        for key, setup in self.radios.items():
            self.radios[key]['button'] = RadioButtons(radioaxes[i],
                        self.radios[key]['labels'], radio_props={'s': 42})
            self.radios[key]['button'].on_clicked(self.radios[key]['func'])
            i += 1
        x = radioaxes[-1].get_position().x0
        y = radioaxes[-1].get_position().y0
        axccmap = self.radioarea.add_axes([0.245, 0.6/plupps, 0.7, 0.034])
        self.tbccmap = TextBox(axccmap, "")
        self.tbccmap.on_submit(self.custom_cmap)
        self.tbccmap.color='lightgrey'
        self.customcmap=False

        """ SET TITLES """
        radioaxes[0].set_title("Ion array species")
        radioaxes[1].set_title("Gas array species")
        radioaxes[2].set_title("Colormap")
        
        self.varbuttonarea.suptitle("Plot variables")
        self.switcharea.suptitle("Plot switches")

        [_, self.varbox, _, self.savebox, self.menu] = \
                self.controlarea.subfigures(\
                        5,
                        1, 
                        wspace=0.1, 
                        height_ratios=[0.5, 1, 0.5, 1,2
        ])
        axep = self.savebox.subplots(1, 1)
        self.tbe = TextBox(axep, "Export path:".rjust(25))
        self.tbe.on_submit(self.exportpath)
        self.tbe.color='lightgrey'
        self.exportpath = expanduser("~/fig.png")
        self.tbe.text_disp.set_text('~/fig.png')
        
        axprompt = self.varbox.subplots(1,1)
        self.tbprompt = TextBox(axprompt, "Plot expression:".rjust(25))
        self.tbprompt.on_submit(self.prompt)
        self.tbprompt.color='lightgrey'

        self.menubuttons = {
            "Export": {'func': self.export},
            "Close": {'func': self.close},
        }
        maxmenus = 10
        [_, self.active_menu] = self.menu.subfigures(1, 2,
            width_ratios=[
                    1-(len(self.menubuttons)/maxmenus),
                    len(self.menubuttons)/maxmenus])
        menuaxes = self.active_menu.subplots(1, len(self.menubuttons))
        i = 0
        for key, setup in self.menubuttons.items():
            self.menubuttons[key]['button'] = \
                Button(menuaxes[i], key, color='lightgrey')
            self.menubuttons[key]['button'].on_clicked(setup['func'])
            i += 1

        for line in self.ax.lines:
            if "child" in line.get_label():
                line.remove()

        xlim = [self.case.get('rm').min(), self.case.get('rm').max()]
        ylim = [self.case.get('zm').min(), self.case.get('zm').max()]
        for var in ['xlim', 'rplate1', 'rplate2']:
            xlim[0] = min(self.case.get(var).min(), xlim[0])
            xlim[1] = max(self.case.get(var).max(), xlim[1])
        for var in ['ylim', 'zplate1', 'zplate2']:
            ylim[0] = min(self.case.get(var).min(), ylim[0])
            ylim[1] = max(self.case.get(var).max(), ylim[1])
        if self.flip is True:
            ylim = [self.case.disp - x for x in ylim[::-1]]


        self.ax.set_xlim((xlim[0]-0.03, xlim[1]+0.03))
        self.ax.set_ylim((ylim[0]-0.03, ylim[1]+0.03))
        self.board.show()
        self.zrange_slider.on_changed(self.update)
        ion()

    def close(self, event):
        from matplotlib.pyplot import close
        close(self.board)
        del self

    def exportpath(self, val):
        from os.path import expanduser
        self.exportpath = expanduser(val)

    def export(self, event):
        from matplotlib.transforms import Bbox
        from matplotlib.pyplot import close
        from copy import deepcopy
        if self.multispecies == 0:
            var = self.var
        elif self.multispecies == 'ion':
            var = self.var[:,:,self.s]
        elif self.multispecies == 'gas':
            var = self.var[:,:,self.ns]
        figsize=self.board.get_size_inches()
        figsize[0] *= 0.42
        bb=Bbox([[0,0],[*figsize]])
        self.board.savefig(self.exportpath, dpi=300, bbox_inches=bb)
        """
        # TODO: Figure out why closing the figure causes
        # a segfault??!
        self.expfig=self.case.plotmesh(
                deepcopy(var),
                cmap=self.cmap,
                figsize=figsize,
                xlim=self.ax.get_xlim(),
                ylim=self.ax.get_ylim(),
                zrange=self.zrange_slider.val,
                log=self.Log,
                vessel=self.Vessel,
                plates=self.Plates,
                lcfs=self.Separatrix,
                title=self.varname,
                grid=self.Grid,
                flip=self.flip
        )
        self.expfig.savefig(self.exportpath, dpi=300)       
#       Closing the new figure segfaults for some reason?!
        close(self.expfig)
        """
        print(f"Figure saved successfully to: {self.exportpath}")
        

    """ RADIO BUTTONS """
        
    def update_ionspecies(self, label):
        self.s = self.case.ionarray.index(label)
        if self.multispecies == 'ion':
            self.update_slider()
            self.update(self.zrange_slider.val)
            self.reset_buttons()

    def update_gasspecies(self, label):
        self.ns = self.case.gasarray.index(label)
        if self.multispecies == 'gas':
            self.update_slider()
            self.update(self.zrange_slider.val)
            self.reset_buttons()

    def update_cmap(self, label):
        self.cmap = label
        if label != '':
            self.verts.set_cmap(self.cmap)
            self.customcmap=False
        else:
            self.customcmap=True
            if not self.fromccmap:
                try:
                    self.custom_cmap(self.tbccmap.text)
                except:
                    pass
        self.f.canvas.draw()


    def custom_cmap(self, val):
        if not self.customcmap:
            self.fromccmap = True
            self.radios['cmap']['button'].set_active(len(\
                self.radios['cmap']['button'].labels)-1)
            self.fromccmap = False
        self.cmap = val
        try:
            self.verts.set_cmap(self.cmap)
        except:
            self.tbccmap.color='salmon'
            return
        
    """ TOGGLES """

    def toggle(self, label):
        status = self.__getattribute__(label)
        self.switches[label]['button'].color = ('light'*status) + 'grey'
        setattr(self, label, not status)
        self.f.canvas.draw()

    def toggle_plates(self, event):
        if self.Plates is True:
            for line in self.ax.lines:
                if 'plate' in line.get_label():
                    line.set_visible(False)
        else:
            for line in self.ax.lines:
                if 'plate' in line.get_label():
                    line.set_visible(True)
        self.toggle('Plates')

    def toggle_vessel(self, event):
        from matplotlib.pyplot import show
        if self.Vessel is True:
            for line in self.ax.lines:
                if line.get_label() == 'vessel':
                    line.set_visible(False)
        else:
            for line in self.ax.lines:
                if line.get_label() == 'vessel':
                    line.set_visible(True)
        self.toggle('Vessel')

    def toggle_sep(self, event):
        if self.Separatrix is True:
            for line in self.ax.lines:
                if line.get_label() == 'lcfs':
                    line.set_visible(False)
        else:
            for line in self.ax.lines:
                if line.get_label() == 'lcfs':
                    line.set_visible(True)
                    line.set_alpha(1)
        self.toggle('Separatrix')

    def toggle_log(self, event):
        from matplotlib.colors import LogNorm, Normalize
        from numpy import log
        lims = self.get_lims()
        for x in lims:
            if x<=0:
                raise ValueError("Cannot plot negative values on log scale!")
        if self.Log is True:
            self.verts.set_norm(Normalize(vmin=lims[0], vmax=lims[1]))
        else:
            self.verts.set_norm(LogNorm(*lims))
            self.verts.set_clim(*lims)
        self.toggle('Log')

    def toggle_grid(self, event):
        if self.Grid is True:
            self.verts.set_linewidths(1)
            self.verts.set_edgecolors("face")
        else:
            self.verts.set_edgecolors("lightgrey")
            self.verts.set_linewidths(0.08)
        self.toggle('Grid')

    def update_lower(self, val):
        self.zrange_slider.set_min(float(val))
        self.update_slider()

    def update_upper(self, val):
        self.zrange_slider.set_max(float(val))
        self.update_slider()

    def display_value(self, box, val):
        valstr = "{:.3g}".format(val)
        box.text_disp.set_text(valstr)



    def update_slider(self):
        from copy import deepcopy
        lims = self.get_lims()
        self.zrange_slider.valmin = lims[0]
        self.zrange_slider.valmax = lims[1]
        vals = deepcopy(self.zrange_slider.val)
        valinit = deepcopy(self.zrange_slider.valinit)
        self.zrange_slider.valinit = vals
        # TODO: only update new_valint when zrange_slider set?
        new_valinit = list(deepcopy(vals))
        # If values within 1/1000th of the resolution, rescale
        if (vals[0] < lims[0]) or (vals[0] > lims[1]) or \
            (vals[0] > lims[0]*1e4):
            self.zrange_slider.set_min(lims[0])
            new_valinit[0]=lims[0]
        if (vals[1] > lims[1]) or (vals[1] < lims[0]) or \
            (vals[1] < lims[1]/1e4):
            self.zrange_slider.set_max(lims[1])
            new_valinit[1]=lims[1]
        self.zrange_slider.valinit = (new_valinit)
        self.zrange_slider.ax.set_ylim(lims)
        self.zrange_slider.set_val(self.zrange_slider.val)
        self.display_value(self.tbllim, self.zrange_slider.val[0])  
        self.display_value(self.tbulim, self.zrange_slider.val[1])  
        self.zrange_slider.reset()

    def get_lims(self):
        if self.multispecies == 0:
            var = self.var
        elif self.multispecies == 'ion':
            var = self.var[:,:,self.s]
        elif self.multispecies == 'gas':
            var = self.var[:,:,self.ns]
        return (var[1:-1,1:-1].min(), 
                var[1:-1,1:-1].max())
 
    def reset_buttons(self):
        for key, setup in self.varbuttons.items():
            setup['button'].color="lightgrey"
        self.tbprompt.color='lightgrey'
        self.tbprompt.text_disp.set_text('')

    def execute_button(self, title, label, multispecies=0):
        self.multispecies = multispecies
        self.varname = title
        self.update_slider()
        self.update(self.zrange_slider.val)
        self.reset_buttons()
        self.varbuttons[label]['button'].color = 'grey'

    def ti(self, event):
        self.var = self.case.get("ti")/1.602e-19
        self.execute_button("Ti [eV]", "Ion temperature")

    def te(self, event):
        self.var = self.case.get("te")/1.602e-19
        self.execute_button("Te [eV]", "Electron temperature")
        return

    def ne(self, event):
        self.var = self.case.get("ne")
        self.execute_button("ne [m**-3]", "Electron density")
        return

    def ni(self, event):
        self.var = self.case.get("ni")
        self.execute_button("ni [m**-3], {}".format(\
                self.case.ionarray[self.s]),
                "Ion density", 
                "ion"
        )
        return

    def ng(self, event):
        self.var = self.case.get("ng")
        self.execute_button("ng [m**-3], {}".format(\
                self.case.gasarray[self.ns]), 
                "Gas density", 
                "gas"
        )
        return

    def prompt(self, command):
        from uedge import bbb, com, grd, flx, aph, api
        # TODO: add log-checks here
        get = self.case.get
        self.multispecies = 0
        self.varname = command
        for substr in ['get', '(', ')', '"', "'", ".", "bbb", "com", 
            "grd", "flx", "aph", "api"]:
            self.varname = self.varname.replace(substr, '')
        try:
            exec( "self.var = " + command)
        except:
            self.tbprompt.color='salmon'
            return
        self.reset_buttons()
        self.tbprompt.text_disp.set_text(command)
        try:
            self.verts.set_array(
                self.var[1:-1, 1:-1].reshape(
                    self.case.get("nx") * self.case.get("ny")
                )
            )
        except:
            self.tbprompt.color='salmon'
            return
        self.update_slider()
        self.update(self.zrange_slider.val)
        self.tbprompt.color='grey'
        # TODO: Use radio buttons for custom prompts?
        # No, species indices may not be straightforward

    def update(self, val):
        from copy import deepcopy
  
        self.verts.set_cmap(self.cmap)
        if self.multispecies == 0:
            var = self.var[1:-1, 1:-1]
        elif self.multispecies == 'ion':
            var = self.var[1:-1, 1:-1, self.s]
        elif self.multispecies == 'gas':
            var = self.var[1:-1, 1:-1, self.ns]
        else:
            raise IndexError("Shape of var could not be determined")
        self.verts.set_array(
            var.reshape(
                self.case.get("nx") * self.case.get("ny")
            )
        )
        self.f.get_axes()[0].set_title(self.varname)
        self.verts.set_clim(self.zrange_slider.val)
        self.display_value(self.tbllim, self.zrange_slider.val[0])  
        self.display_value(self.tbulim, self.zrange_slider.val[1])  


