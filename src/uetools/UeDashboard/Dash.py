import sys
from uetools import Case
import matplotlib
matplotlib.use('Qt5Agg')

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib import rcParams

#from range_slider import RangeSlider
from PyQt5.QtCore import Qt, QMargins
from PyQt5.QtWidgets import QMenu
from PyQt5.QtGui import QFont, QDoubleValidator
from range_slider import RangeSlider
from PyQt5.QtWidgets import QAction, QWidget
from functools import partial
from PyQt5.QtWidgets import (
    QApplication, 
    QLabel, 
    QMainWindow,
    QMenu,
    QAction,
    QSlider,
    QVBoxLayout,
    QHBoxLayout,
    QGridLayout,
    QPushButton,
    QLineEdit,
    QRadioButton,
    QButtonGroup,
    QCheckBox,
    QFrame,
    QFormLayout,
    QComboBox,
)


# TODO: Text boxes to set limits
# TODO: Open dialogue
    # update permanent bar


# TODO: Store previously opened files in separate yaml under .uedgerc
# TODO: Make dialogue work with active UEDGE run
# TODO: make raise_message bolded red in bar


# TODO: How to deal with CX/psorx?
# TODO: Set zeros to 1e-10

class CaseDashboard(QMainWindow):
    """Main Window."""
    def __init__(self, case, parent=None):
        """Initializer."""
        super().__init__(parent)
        self.file = case.filename
        self.case = case        

        self.setFocus()
        self.setWindowTitle("UETOOLS Case Heatmap")
        self.resize(1300, 900)

        # TODO Initialize to dummy file to get prompt open?
        self.layout = QGridLayout()

        self._createCanvas()
        self._createVarButtons()
        self._createSlider()
        self._createSwitches()
        self._createRadios()
    
        self._createActions()
        self._createMenuBar()
        self._connectActions()
        self._createStatusBar()
        self._createMiscOptions()
        self._createSettings()

        self.layout.addLayout(self.buttons['layout'],    
                0, 0, 5, 2
        )
        self.layout.addLayout(self.radios['layout'],
                0, 2, 5, 1
        )
        self.layout.addLayout(self.canvaslayout, 0, 5, 17, 18)
        self.layout.addLayout(self.settings, 7, 0, 7, 3)
        self.layout.addWidget(self.slider['items']['ulim'],
                0, 23, 1, 3
        )
        self.layout.addWidget(self.slider['items']['slider'],        
                1, 24, 15, 1, Qt.AlignHCenter
        )
        self.layout.addWidget(self.slider['items']['llim'],
                16, 23, 1, 3
        )

        self.centralWidget = QWidget(self)
        self.centralWidget.setLayout(self.layout)

        self.setCentralWidget(self.centralWidget)


    def _createSettings(self):
        self.settings = QVBoxLayout()
        decks = QHBoxLayout()
        self.settings.addLayout(decks)
        decks.addWidget(self.misc['frame'])
        decks.addWidget(self.cmap_radio['frame'])
        decks.addWidget(self.switches['frame'])
        

    def _createCanvas(self):
        self.canvas = MplCanvas(self, width=5, height=4, dpi=100)
        self.case.te2D(ax=self.canvas.axes)#, flip=True)
        self.case.set_speciesarrays()
        self.verts = self.case.Qvertices
        self.suptitle = 'Electron temperature [eV]'
        self.canvas.fig.suptitle(self.suptitle)
        
        self.canvastoolbar = NavigationToolbar2QT(self.canvas, self)
        self.canvaslayout = QVBoxLayout()
        self.canvaslayout.addWidget(self.canvastoolbar)
        self.canvaslayout.addWidget(self.canvas)

        self.gasspecies = 0
        self.ionspecies = 0
        self.cmap = 'magma'
        self.grid = False
        self.log = False
        self.abs = False
        self.varscale = 1
        self.var = self.case.get("te")/1.602e-19
        self.multispecies = False
        self.xy = self.case.get("nx") * self.case.get("ny")
        for line in self.canvas.axes.lines:
            if "child" in line.get_label():
                line.remove()


    def _createMiscOptions(self):
        self.misc = {
            'frame': QFrame(),
            'items': {
                'abs': QCheckBox(),
                'varscale': MyLineEdit("1"),
                'ulim': MyClearLineEdit(),
                'llim': MyClearLineEdit(),
            }
        }
        tmp = self.misc['items']

        form =  QFormLayout()
        form.addRow("Variable absolute value", tmp["abs"])
        tmp['abs'].toggled.connect(self.toggle_abs)
        tmp['abs'].setToolTip("Plots the absolute value of variables.")

        tmp['varscale'].setValidator(QDoubleValidator())
        tmp['varscale'].returnPressed.connect(self.change_varscale)
        tmp['varscale'].setFixedWidth(50)
        form.addRow("Variable scaling factor", tmp["varscale"])
        tmp['varscale'].mousePressEvent = lambda _ : tmp['varscale'].selectAll()
        tmp['varscale'].setToolTip("Hit return to apply.")

        tmp['ulim'].setValidator(QDoubleValidator())
        tmp['ulim'].returnPressed.connect(self.change_ulim)
        tmp['ulim'].setFixedWidth(50)
        tmp['ulim'].setToolTip("Hit return to apply.")
        form.addRow("Set upper plot range", tmp ['ulim'])

        tmp['llim'].setValidator(QDoubleValidator())
        tmp['llim'].returnPressed.connect(self.change_llim)
        tmp['llim'].setFixedWidth(50)
        tmp['llim'].setToolTip("Hit return to apply.")
        form.addRow("Set lower plot range", tmp ['llim'])

        self.misc['frame'].setFrameStyle(QFrame.Panel | QFrame.Raised)
        layout = QVBoxLayout()
        layout.addWidget(QLabel(self.title("Plot options")))
        layout.addLayout(form)
        layout.addStretch()
        self.misc['frame'].setLayout(layout)


    def set_custom_cmap(self):
        cmap = self.cmap_radio['items']['cmap'].text()
        try:
            self.verts.set_cmap(cmap)
            self.cmap = cmap
            self.cmap_radio['group'].setExclusive(False)
            for key, item in self.cmap_radio['items'].items():
                try:
                    item.setChecked(False)
                except:
                    pass
            self.cmap_radio['group'].setExclusive(True)
        except:
            self.raise_message(f"Colormap {cmap} not available!")
            self.cmap_radio['items']['cmap'].clear()
        self.canvas.draw()
        self.cmap_radio['items']['cmap'].clearFocus()


    def toggle_abs(self):
        self.abs = not self.abs
        self.update_var()

    def change_varscale(self):
        self.varscale = float(self.misc['items']['varscale'].text())
        self.update_var()
        self.misc['items']['varscale'].clearFocus()

    def change_ulim(self):
        ulim = float(self.misc['items']['ulim'].text())
        clim = self.verts.get_clim()
        lims = self.get_lims()
        if ulim <= clim[0]:
            self.raise_message("Requested upper range {:.3g} is".format(ulim)+
                " below the current lower limit {:.3g}".format(clim[0]))
            return
        ulim = min(ulim, lims[1])
        self.verts.set_clim((clim[0], ulim))
        self.slider['items']['ulim'].setText(\
                self.title("{:.3g}".format(ulim)
        ))
        self.slider['items']['slider'].setHigh(
           int( self.value2position(ulim)
        ))
        self.update_var()
        self.misc['items']['ulim'].clear()
        self.misc['items']['ulim'].clearFocus()

    def change_llim(self):
        llim = float(self.misc['items']['llim'].text())
        clim = self.verts.get_clim()
        lims = self.get_lims()
        if llim >= clim[1]:
            self.raise_message("Requested lower range {:.3g} is".format(llim)+
                " above the current upper limit {:.3g}".format(clim[1]))
            return
        llim = max(llim, lims[0])
        self.verts.set_clim((llim, clim[1]))
        self.slider['items']['llim'].setText(\
                self.title("{:.3g}".format(llim)
        ))
        self.slider['items']['slider'].setLow(
           int( self.value2position(llim)
        ))
        self.update_var()
        self.misc['items']['llim'].clear()
        self.misc['items']['llim'].clearFocus()





    def _createRadios(self):
        self._createIonRadio()
        self._createGasRadio()
        self._createCmapRadio()
        self.radios = {'layout': QVBoxLayout()}

        self.radios['layout'].addWidget(self.ion_radio['frame'])
        self.radios['layout'].addWidget(self.gas_radio['frame'])
        self.radios['layout'].addStretch()

    def _createIonRadio(self):
        self.ion_radio = {
            'frame': QFrame(),
            'layout': QVBoxLayout(),
            'group': QButtonGroup(),
            'title': QLabel(self.title("Ion species")),
            'labels': self.case.ionarray,
            'items': {}
        }
        self.ion_radio['title'].setAlignment(\
                    Qt.AlignHCenter |\
                    Qt.AlignVCenter
        )
        self.ion_radio['layout'].addWidget(self.ion_radio['title'])
        ind = 0
        for label in self.ion_radio['labels']:
            self.ion_radio['items'][label] = QRadioButton(label)
            self.ion_radio['layout'].addWidget(\
                    self.ion_radio['items'][label])
            self.ion_radio['group'].addButton(\
                    self.ion_radio['items'][label], ind)
            self.ion_radio['items'][label].clicked.connect(self.ion_radio_clicked)
            ind += 1
        self.ion_radio['frame'].setFrameStyle(QFrame.Panel | QFrame.Raised)
        self.ion_radio['frame'].setLayout(self.ion_radio['layout'])
        self.ion_radio['items'][self.case.ionarray[0]].setChecked(True)

    def _createGasRadio(self):
        self.gas_radio = {
            'frame': QFrame(),
            'layout': QVBoxLayout(),
            'group': QButtonGroup(),
            'title': QLabel(self.title("Gas species")),
            'labels': self.case.gasarray,
            'items': {}
        }
        self.gas_radio['title'].setAlignment(\
                    Qt.AlignHCenter |\
                    Qt.AlignVCenter
        )
        self.gas_radio['layout'].addWidget(self.gas_radio['title'])
        ind = 0
        for label in self.gas_radio['labels']:
            self.gas_radio['items'][label] = QRadioButton(label)
            self.gas_radio['layout'].addWidget(\
                    self.gas_radio['items'][label])
            self.gas_radio['group'].addButton(\
                    self.gas_radio['items'][label], ind)
            self.gas_radio['items'][label].clicked.connect(self.gas_radio_clicked)
            ind += 1
        self.gas_radio['frame'].setFrameStyle(QFrame.Panel | QFrame.Raised)
        self.gas_radio['frame'].setLayout(self.gas_radio['layout'])
        self.gas_radio['items'][self.case.gasarray[0]].setChecked(True)


    def _createCmapRadio(self):
        self.cmap_radio = {
            'frame': QFrame(),
            'layout': QVBoxLayout(),
            'group': QButtonGroup(),
            'title': QLabel(self.title("Colormaps")),
            'labels': [
                'magma',
                'bwr',
                'cividis',
                'hsv',
                'Dark2',
                'gnuplot2',
                'hot',
                'jet',
            ],
            'items': {
                'cmap': MyLineEdit()
            }
        }
        self.cmap_radio['layout'].addWidget(self.cmap_radio['title'])
        ind = 0
        for label in self.cmap_radio['labels']:
            self.cmap_radio['items'][label] = QRadioButton(label)
            self.cmap_radio['layout'].addWidget(\
                    self.cmap_radio['items'][label])
            self.cmap_radio['group'].addButton(\
                    self.cmap_radio['items'][label], ind)
            self.cmap_radio['items'][label].clicked.connect(self.cmap_radio_clicked)
            ind += 1
        self.cmap_radio['frame'].setFrameStyle(QFrame.Panel | QFrame.Raised)
        self.cmap_radio['frame'].setLayout(self.cmap_radio['layout'])
        self.cmap_radio['items'][self.cmap_radio['labels'][0]].setChecked(True)

        self.cmap_radio['layout'].addWidget(self.cmap_radio['items']["cmap"])
        self.cmap_radio['items']['cmap'].returnPressed.connect(self.set_custom_cmap)
        self.cmap_radio['items']['cmap'].setFixedWidth(100)
        self.cmap_radio['items']['cmap'].mousePressEvent = \
                lambda _ : self.cmap_radio['items']['cmap'].selectAll()
        self.cmap_radio['items']['cmap'].setToolTip("Hit return to apply.")
        self.cmap_radio['layout'].addStretch()


    def cmap_radio_clicked(self):
        self.cmap = self.cmap_radio['group'].checkedButton().text()
        self.verts.set_cmap(self.cmap)
        self.cmap_radio['items']['cmap'].clear()
        self.canvas.draw()

    def ion_radio_clicked(self):
        self.ionspecies = self.ion_radio['group'].checkedId()
        title = self.suptitle
        if self.multispecies == 'ion':
            if " for " in title:
                title = title.split("for")[0] + "for {}".format(
                    self.case.ionarray[self.ionspecies])
            self.suptitle = title
            self.update_var()
        
    def gas_radio_clicked(self):
        self.gasspecies = self.gas_radio['group'].checkedId()
        if self.multispecies == 'gas':
            title = self.canvas.fig._suptitle.get_text()
            if " for " in title:
                title = title.split("for")[0] + "for {}".format(
                    self.case.gasarray[self.gasspecies])
            self.canvas.fig.suptitle(title)
            self.update_var()

    def _createStatusBar(self):
        self.statusbar = self.statusBar()


    def _createSlider(self):
        lims = self.get_lims()
        self.slider = {
            'layout': QVBoxLayout(),
            'items': {
                'ulim': QLabel(self.title("{:.3g}".format(lims[1]))),
                'slider': RangeSlider(Qt.Vertical),
                'llim': QLabel(self.title("{:.3g}".format(lims[0]))),
            }
        }

        self.slider['items']['llim'].setAlignment(\
                        Qt.AlignHCenter | \
                        Qt.AlignVCenter
        )
        self.slider['items']['ulim'].setAlignment(\
                        Qt.AlignHCenter | \
                        Qt.AlignVCenter
        )

        self.slider['items']['slider'].setMinimum(0)
        self.slider['items']['slider'].setMaximum(1000)
        self.slider['items']['slider'].setLow(0)
        self.slider['items']['slider'].setHigh(1000)
        self.slider['items']['slider'].setMinimumHeight(30)
        self.slider['items']['slider'].setValue(1000)
        self.slider['items']['slider'].sliderMoved.connect(\
                self.update_plot_range)


    def _createSwitches(self):
        self.switches = {
            'frame': QFrame(),
            'layout': QVBoxLayout(),
            'title': QLabel(self.title("Switches")),
            'items': {}
        }
        
        grid = QHBoxLayout()
        col = [QVBoxLayout(), QVBoxLayout(), QVBoxLayout()]
        grid.addLayout(col[0])
        grid.addLayout(col[1])
        grid.addLayout(col[2])
        col[0].setContentsMargins(10,0,10,5)
        col[1].setContentsMargins(10,0,10,5)

        self.switches['layout'].addWidget(self.switches['title'])
        self.switches['layout'].addLayout(grid)
        self.switches['layout'].setSpacing(5)
        grid.setSpacing(5)
        self.switches['frame'].setLayout(self.switches['layout'])
        self.switches['frame'].setFrameStyle(QFrame.Panel | QFrame.Raised)

        i = 0
        for switch in ["Vessel", "Plates", "Separatrix", "Grid", "Scale"]:
            if switch != "Scale":
                on, off = "On", "Off"
            else:
                on, off = "Lin", "Log"
            self.switches['items'][switch.lower()] = {
                    'group': QButtonGroup(),
                    on.lower(): QRadioButton(on),
                    off.lower(): QRadioButton(off),
                    'layout': QVBoxLayout()
            }
            tmpsw = self.switches['items'][switch.lower()]
            tmpsw[on.lower()].toggled.connect(self.__getattribute__(\
                "{}_switch_clicked".format(switch.lower())))
            tmpsw['group'].addButton(tmpsw[on.lower()], 1)
            tmpsw['group'].addButton(tmpsw[off.lower()], 0)
            tmpsw['layout'].addWidget(QLabel(f"<b>{switch}</b>"))
            tmpsw['layout'].addWidget(tmpsw[on.lower()])
            tmpsw['layout'].addWidget(tmpsw[off.lower()])
            tmpsw['layout'].setContentsMargins(0,0,0,0)
            tmpsw['layout'].setSpacing(0)
            col[int(i/5)].addLayout(tmpsw['layout'])
            i += 1
        col[-1].addStretch()
        self.switches['items']['vessel']['on'].setChecked(True)
        self.switches['items']['plates']['on'].setChecked(True)
        self.switches['items']['separatrix']['on'].setChecked(True)
        self.switches['items']['grid']['off'].setChecked(True)
        self.switches['items']['scale']['lin'].setChecked(True)
        self.plates_switch_clicked()
        self.vessel_switch_clicked()
        self.separatrix_switch_clicked()
        self.scale_switch_clicked()


    def plates_switch_clicked(self):
        for line in self.canvas.axes.lines:
            if 'plate' in line.get_label():
                line.set_visible(not line.get_visible())
        self.canvas.draw()

    def grid_switch_clicked(self):
        if self.grid:
            self.verts.set_edgecolor('face')
            self.verts.set_linewidths(1)
        else:
            self.verts.set_edgecolor('lightgrey')
            self.verts.set_linewidths(0.08)
        self.grid = not self.grid
        self.canvas.draw()

    def separatrix_switch_clicked(self):
        for line in self.canvas.axes.lines:
            if line.get_label() == "lcfs":
                line.set_visible(not line.get_visible())
        self.canvas.draw()

    def vessel_switch_clicked(self):
        for line in self.canvas.axes.lines:
            if line.get_label() == 'vessel':
                line.set_visible(not line.get_visible())
        self.canvas.draw()

    def scale_switch_clicked(self):
        from matplotlib.colors import LogNorm, Normalize
        from numpy import log
        lims = self.get_lims()

        '''
        self.log = not self.log
        if self.log:
            self.verts.set_norm(Normalize(*lims))#vmin=lims[0], vmax=lims[1])
        else:
            self.verts.set_norm(LogNorm(*lims))
        self.update_var()  

        return
        '''
        if (lims[0] <= 0) and (self.log is False):
            self.raise_message("Negative values in var: masking array!")
            var = self.verts.get_array()
            var.mask=(var<=0)
            self.verts.set_array(var)
            lims = self.get_lims()
        else:
            var = self.verts.get_array()
            var.mask=False
            self.verts.set_array(var)
        if self.log:
            self.verts.set_norm(Normalize(*lims))#vmin=lims[0], vmax=lims[1])
        else:
            self.verts.set_norm(LogNorm(*lims))
        vals = self.get_slider_values()
        self.log = not self.log 
        newpos = [int(self.value2position(x)) for x in vals]
        self.slider['items']['slider'].setLow(newpos[0])
        self.slider['items']['slider'].setHigh(newpos[1])
        clim = list(self.verts.get_clim())
        if self.log:
            if clim[0] < 0:
                clim[0]=1e-100
            if (clim[1] <0) or (clim[1]<clim[0]):
                clim[1]=1e100
        self.verts.set_clim(clim)
        self.canvas.draw()


    def _createVarButtons(self):
        self.buttons = {
            'layout': QGridLayout(),
            'title': QLabel(self.title("Plots")),
            'items': {
                'te': 'Electron temperature',
                'ti': 'Ion temperature',
                'tg': 'Gas temperature',
                'ne': 'Electron density',
                'ni': 'Ion density',
                'ng': 'Gas density',
                'phi': 'Potential',
                'prad': 'Total radiated power',
                'pradhyd': 'Hyd. radiated power',
                'pradimp': 'Imp. radiated power',
                'psorc': "Ion ioniz. source",
                'psorgc': "Gas ioniz. sink",
                'psorxrc': "Ion rec. + CX sink",
                'psorrgc': "Gas rec. source",
                'psorcxg': "Gas CX source",
            },
        }
        self.buttons["title"].setAlignment(Qt.AlignHCenter)
        self.buttons['layout'].addWidget(self.buttons['title'],0,0,1,2)
        self.buttons['layout'].setSpacing(0)

        cols = [QVBoxLayout()]
        col = cols[-1]
        i = 10
        j = 0
        for key, setup in self.buttons['items'].items():
            if i > 7:   
                i = 0
                cols.append(QVBoxLayout())
                col = cols[-1]
                self.buttons['layout'].addLayout(col, 1, j, 1, 1)
                j +=1
            self.buttons['items'][key] = QPushButton(setup)
            self.buttons['items'][key].clicked.connect(
                self.__getattribute__(f"plot_{key}"))
            col.addWidget(self.buttons['items'][key])
            i += 1
        # Add drop-down
        self.buttons['items']['dropdown'] = QComboBox()
        varlist = list(self.buttons['items'].keys())
        self.buttons['items']['dropdown'].addItem("", 0)
        for vartype in ['centered', 'staggered']:
            if self.case.inplace:
                for var, path in self.case.vars.items():
                    if (vartype in path) and (var not in varlist):
                        self.buttons['items']['dropdown'].addItem(var, i)
        cols[-1].addWidget(self.buttons['items']['dropdown']) 
        self.buttons['items']['dropdown'].activated[str].connect(self.plot_dropdown)
        for col in cols:
            col.addStretch()

        self.buttons['layout'].addWidget(QLabel("Custom formula"), 2, 0, 1, 2)
        self.buttons['items']['custom'] = MyLineEdit()
        self.buttons['items']['custom'].mousePressEvent = \
            lambda _ : self.buttons['items']['custom'].selectAll()
        self.buttons['items']['custom'].returnPressed.connect(self.plot_custom)
        self.buttons['items']['custom'].setToolTip("Plots custom formula. "+
            "Access var by get('var'). Python syntax applies. Expects output"+
            "of shape ({},{}). Hit return to apply. ".format(self.case.nx, 
            self.case.ny))
        self.buttons['layout'].addWidget(self.buttons['items']['custom'], 
            3, 0, 1, 2)

    def plot_custom(self):
        try:
            from uedge import com, bbb, grd, flx, aph, api
        except:
            pass
        command = self.buttons['items']['custom'].text()
        self.buttons['items']['dropdown'].setCurrentIndex(0)
        get = self.case.get
        try:
            exec("locals()['var'] =" + command)
            locals()['var'].shape # Fails non-ndarray results
        except:
            self.raise_message(f"{command} not recognized!")
            return
        if len(locals()['var'].shape) > 2:
            self.raise_message("Command variable has wrong shape: "+
                "{}, expected ({},{})".format(\
                    locals()['var'].shape,
                    self.case.nx+2,
                    self.case.ny+2
                )
            )
        else:
            self.var = locals()['var']
            for substr in ['get', '(', ')', '"', "'", ".", "bbb", "com",
                "grd", "flx", "aph", "api", "self"]:
                command = command.replace(substr, '')
            self.suptitle = command
            self.buttons['items']['custom'].clearFocus()
            self.multispecies = False
            self.update_var()


    def plot_dropdown(self, text):
        self.raise_message(text)
        var = self.case.get(text)
        if len(var.shape) == 3:
            if var.shape[-1] == len(self.case.ionarray):
                self.multispecies = 'ion'
                self.suptitle = text + " for {}".format(\
                    self.case.ionarray[self.ionspecies])
                self.var = var
            elif var.shape[-1] == len(self.case.gasarray):
                self.multispecies = 'gas'
                self.suptitle = text + " for {}".format(\
                    self.case.gasarray[self.gasspecies])
                self.var = var
            else:
                self.raise_message("Could not determine species"+
                    f" index of {text}")
                return
        elif len(var.shape) == 2:
            self.var = var
            self.suptitle = text
        else:
            self.raise_message(f"Could not determine shape of {text}")
            return
            
        self.buttons['items']['custom'].clear()
        self.update_var()

    def update_var(self):
        from copy import deepcopy
        from numpy import sum
        if self.multispecies == 'ion':
            var = self.var[1:-1,1:-1,self.ionspecies]
        elif self.multispecies == 'gas':
            var = self.var[1:-1,1:-1,self.gasspecies]
        else:
            var = self.var[1:-1,1:-1]
        var *= self.varscale
        if sum(var) == 0:
            self.raise_message("Requested variable unpopulated! "+
                "Select another variable or species to plot.")
            return
        if self.abs:
            var = abs(var)
        if (var.min() < 0) and self.log:
            self.raise_message("Negative values in var: masking array!")
            var.mask=(var<=0)
        # Record old values slider values
        [llim, ulim] = deepcopy(self.verts.get_clim())
        self.verts.set_array(var.reshape(self.xy))
        clim = self.verts.get_clim()
        lims = self.get_lims()
        if ulim <= lims[0]:
            ulim = lims[1]
        ulim = min(ulim, lims[1])
        if llim >= lims[1]:
            llim = lims[0]
        llim = max(llim, lims[0])
        # Enseure some separation between sliders when var is constant
        if lims[1]==lims[0]:
            llim, ulim = 1, 1e100
        # Enseure some separation between sliders when bounds are maintained
        if abs(ulim-llim)/max(1e-10, abs(lims[1]-lims[0])) < 1e-2:
            llim, ulim = lims[0], lims[1]
        # NOTE
        # I HAVE NO CLUE WHY THE LOWER LIMIT DOES NOT REGISTER
        # AT THE FIRST CET_CLIM CALL, AND DONT ASK ME HOW LONG
        # IT TOOK TO FIGURE IT OUT!!!
        self.verts.set_clim((llim, ulim))
        self.verts.set_clim((llim, ulim))
        self.slider['items']['ulim'].setText(\
                self.title("{:.3g}".format(ulim)
        ))
        self.slider['items']['llim'].setText(\
                self.title("{:.3g}".format(llim)
        ))
        self.slider['items']['slider'].setHigh(
            int( self.value2position(ulim)
        ))
        self.slider['items']['slider'].setLow(
            int( self.value2position(llim)
        ))
        self.verts.set_cmap(self.verts.get_cmap())
        self.canvas.fig.suptitle(self.suptitle)
        self.canvas.draw()

    def plot_driver(self, var, title, multispecies=False):
        from numpy.ma import masked_array
        # Save variable as masked_array to mask for log
        self.var = masked_array(var, mask=False)
        self.multispecies = multispecies
        species=''
        if multispecies == 'gas':
            species = self.case.gasarray[self.gasspecies]
        elif multispecies == 'ion':
            species = self.case.ionarray[self.ionspecies]
        self.suptitle = title + f" for {species}"*(multispecies is not False)
        self.update_var()
        self.buttons['items']['dropdown'].setCurrentIndex(0)
        command = self.buttons['items']['custom'].clear()
        
    def plot_te(self):
        self.plot_driver(
            self.case.get('te')/1.602e-19,
            'Electron temperature [eV]',
        )

    def plot_ti(self):
        self.plot_driver(
            self.case.get('ti')/1.602e-19,
            'Ion temperature [eV]',
        )

    def plot_tg(self):
        self.plot_driver(
            self.case.get('tg')/1.602e-19,
            'Gas temperature [eV]',
            'gas'
        )

    def plot_ne(self):
        self.plot_driver(
            self.case.get('ne'),
            r'Electron density [m$\mathrm{{}^{-3}}$]',
        )

    def plot_ni(self):
        self.plot_driver(
            self.case.get('ni'),
            r'Ion density [m$\mathrm{{}^{-3}}$]',
            'ion'
        )

    def plot_ng(self):
        self.plot_driver(
            self.case.get('ng'),
            r'Gas density [m$\mathrm{{}^{-3}}$]',
            'gas'
        )

    def plot_phi(self):
        self.plot_driver(
            self.case.get('phi'),
            'Electrical potential [V]',
        )

    def plot_prad(self):
        self.plot_driver(
            self.case.get('prad')+self.case.get('pradhyd'),
            'Total radiated power [W/m$\mathrm{{}^{-3}}$]',
        )

    def plot_pradhyd(self):
        self.plot_driver(
            self.case.get('pradhyd'),
            'Hydrogenic radiated power [W/m$\mathrm{{}^{-3}}$]',
        )

    def plot_pradimp(self):
        self.plot_driver(
            self.case.get('prad'),
            'Impurity radiated power [W/m$\mathrm{{}^{-3}}$]',
        )

    def plot_psorc(self):
        self.plot_driver(
            self.case.get('psorc'),
            'Ion ionization source [parts/s]',
            'ion'
        )

    def plot_psorgc(self):
        self.plot_driver(
            self.case.get('psorgc'),
            'Gas ionization sink [parts/s]',
            'gas'
        )

    def plot_psorxrc(self):
        self.plot_driver(
            self.case.get('psorxrc'),
            'Ion recombination + CX sink [parts/s]',
            'ion'
        )

    def plot_psorrgc(self):
        self.plot_driver(
            self.case.get('psorrgc'),
            'Gas recombination source [parts/s]',
            'gas'
        )

    def plot_psorcxg(self):
        self.plot_driver(
            self.case.get('psorcxg'),
            'Gas CX source [parts/s]',
            'gas'
        )


    def _createMenuBar(self):
        menuBar = self.menuBar()
        menuBar.setNativeMenuBar(False)
        # File menu
        fileMenu = QMenu("&File", self)
        menuBar.addMenu(fileMenu)
#        fileMenu.addAction(self.openAction)
#        openMenu.addAction(self.openYamlAction)
#        openMenu.addAction(self.openHdf5Action)
#        openMenu.addAction(self.openSaveAction)
#        self.openRecentMenu = fileMenu.addMenu("Open Recent")
        fileMenu.addSeparator()
        fileMenu.addAction(self.exitAction)
#        editMenu = menuBar.addMenu("&Edit")
        helpMenu = menuBar.addMenu("&Help")
        helpMenu.addAction(self.helpContentAction)
        helpMenu.addAction(self.aboutAction)

    def _createActions(self):
        # Creating actions using the second constructor
        self.openAction = QAction("Open HDF5 save file...", self)
        self.exitAction = QAction("&Exit", self)
        self.helpContentAction = QAction("&Help Content", self)
        self.aboutAction = QAction("&About", self)
        
    def _createStatusBar(self):
        self.statusbar = self.statusBar()
        self.wcLabel = QLabel(f"{self.getCurrentFile()}")
        self.statusbar.addPermanentWidget(self.wcLabel)

    def _connectActions(self):
        # Connect File actions
        self.openAction.triggered.connect(self.openFile)
        self.exitAction.triggered.connect(self.close)
#        self.openRecentMenu.aboutToShow.connect(self.populateOpenRecent)

    """ ===========================
               SLIDER HELPERS 
        ==========================="""
    def get_lims(self):
        var = self.verts.get_array()
        return (var.min(), var.max())

    def get_slider_position(self):
        return (self.slider_lower_position(), self.slider_higher_position())

    def get_slider_values(self):
        pos = self.get_slider_position()
        return  (max(min(self.position2value(pos[0]), 1000), 0),
                 max(min(self.position2value(pos[1]), 1000),0))

    def slider_higher_position(self):
        return self.slider['items']['slider'].high()

    def slider_lower_position(self):
        return self.slider['items']['slider'].low()

    def position2value(self, pos):
        lims = self.get_lims()
        pos /= 1000
        if self.log:
            return (lims[1]**pos)*(lims[0]**(1-pos))
        else:
            return pos*lims[1] + (1-pos)*lims[0]

    def value2position(self, val):
        from numpy import log
        from numpy.ma import is_masked
        lims = self.get_lims()
        if abs(lims[0]-lims[1])<1e-20:
            return 0
        for lim in lims:
            if is_masked(lim):
                return 0
        if self.log:
            return 1000/(1+(log(lims[1])-log(val))/max(1e-20, 
                        (log(val)-log(lims[0]))))
        else:
            return 1000/(1 + ((lims[1]-val)/max(1e-10,(val-lims[0]))))
        
    def update_plot_range(self, l, u):
        from numpy.ma import is_masked
        nl = self.position2value(l)
        nu = self.position2value(u)
        for var in [nl, nu]:
            if is_masked(var):
                return 0
        self.verts.set_clim((nl, nu))
        self.verts.set_cmap(self.verts.get_cmap())
        self.slider['items']['llim'].setText(self.title("{:.3g}".format(nl)))
        self.slider['items']['ulim'].setText(self.title("{:.3g}".format(nu)))
        self.canvas.draw()


    def getCurrentFile(self):
        return self.file

    def populateOpenRecent(self):
        # Step 1. Remove the old options from the menu
        self.openRecentMenu.clear()
        # Step 2. Dynamically create the actions
        actions = []
        filenames = [f"File-{n}" for n in range(5)]
        for filename in filenames:
            action = QAction(filename, self)
            action.triggered.connect(partial(self.openRecentFile, filename))
            actions.append(action)
        # Step 3. Add the actions to the menu
        self.openRecentMenu.addActions(actions)

    def openRecentFile(self, filename):
        # Logic for opening a recent file goes here...
        self.promptbox.setText(f"<b>{filename}</b> opened")

    def openFile(self):
        # Logic for opening an existing file goes here...
        self.raise_message("<b>File > Open YAML input...</b> clicked")

    def openYamlFile(self):
        # Logic for opening an existing file goes here...
        self.promptbox.setText("<b>File > Open YAML input...</b> clicked")

    def openHdf5File(self):
        # Logic for opening an existing file goes here...
        self.promptbox.setText("<b>File > Open HDF5 input...</b> clicked")

    def openSaveFile(self):
        # Logic for opening an existing file goes here...
        self.promptbox.setText("<b>File > Open Save file...</b> clicked")

    """ ===========================
                FORMATTERS 
        ==========================="""
    def title(self, text):
        return f"<b><font size=4>{text}</font></b>"

    def show_message(self, message, time=5000):    
        self.statusbar.showMessage(message, time)

    def raise_message(self, message, time=5000):
        self.statusbar.showMessage(message, time)

    def mousePressEvent(self, event):
        focused_widget = QApplication.focusWidget()
        if isinstance(focused_widget, MyLineEdit):
            focused_widget.clearFocus()
        if isinstance(focused_widget, MyClearLineEdit):
            focused_widget.clear()
            focused_widget.clearFocus()
        QMainWindow.mousePressEvent(self, event)


class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=2, height=4, dpi=300):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MplCanvas, self).__init__(self.fig)

class MyLineEdit(QLineEdit):
        pass

class MyClearLineEdit(QLineEdit):
        pass





if __name__ == "__main__":
    app = QApplication(sys.argv)
    try:
        win = CaseDashboard(Case(sys.argv[1], inplace=True))
    except:
#        raise ValueError("No file supplied. Aborting!")
        file = "/Users/holm10/Documents/fusion/uedge/src/"+\
                "UETOOLS/jupyter/testcase_hires/dashtest.hdf5"
        print("USING TUTORIAL CASE")
        win = CaseDashboard(Case(file, inplace=True))
    win.show()
    sys.exit(app.exec_())
else:
    app = QApplication(sys.argv)
    win = CaseDashboard(self)
    
     
    win.show()
    sys.exit(app.exec_())
