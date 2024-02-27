import sys
from uetools import Case
import matplotlib
matplotlib.use('Qt5Agg')

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib import rcParams

#from range_slider import RangeSlider
from PyQt5.QtCore import (
    Qt, 
    QMargins,
    QTimer,
    QSize,
)
from PyQt5.QtWidgets import QMenu
from PyQt5.QtGui import (
    QFont, 
    QDoubleValidator,
    QStatusTipEvent,
)
from range_slider import RangeSlider
from functools import partial
from PyQt5.QtWidgets import (
    QAction,
    QWidget,
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
    QTabWidget,
    QRadioButton,
    QButtonGroup,
    QCheckBox,
    QFrame,
    QFormLayout,
    QSizePolicy,
    QComboBox,
    QFileDialog,
    QDialog,
    QDialogButtonBox,
    QInputDialog,
)


    # update permanent bar
# TODO: Display file name somewhere else


# TODO: Store previously opened files in separate yaml under .uedgerc
# TODO: Make dialogue work with active UEDGE run
# TODO: make raise_message bolded red in bar


# TODO: How to deal with CX/psorx?
# TODO: Set zeros to 1e-10


class CaseDashboard(QWidget):
    """Main Window."""
    def __init__(self, case, parent=None):
        """Initializer."""
        super().__init__(parent)
        self.file = case.filename

        self.caseplot = HeatmapInteractiveFigure(case)
        # Cannibalize Case functions
        self.inplace = case.inplace
        self.casevars = case.vars
        self.get = case.get
        self.nx = case.nx
        self.ny = case.ny
        self.ionarray = case.ionarray
        self.gasarray = case.gasarray
        self.casename = case.casename
        self.ionspecies = 0
        self.gaspecies = 0
        self.log = False


        self._createVarButtons()
        self._createSwitches()
        self._createRadios()
        self._createMiscOptions()

        self.plot_te()
        self.caseplot.setTitle("Electron temperature [eV]")

        self.settings = QGridLayout()
        self.settings.addWidget(self.misc['frame'], 0, 0, 1, 1)
        self.settings.addWidget(self.cmap_radio['frame'], 1, 0, 1, 1)
        self.settings.addWidget(self.switches['frame'], 0, 1, 2, 1)

        self.menu = QGridLayout()
        self.menu.addLayout(self.buttons['layout'],    
                0, 0, 1, 2
        )
        self.menu.addLayout(self.radios['layout'],
                0, 2, 1, 1
        )
        self.menu.addLayout(self.settings, 1, 0, 1, 3)

        '''
        self.layout.addLayout(self.canvaslayout, 0, 5, 17, 30)
        self.layout.addWidget(self.slider['items']['ulim'],
                0, 35, 1, 3
        )
        self.layout.addWidget(self.slider['items']['slider'],        
                1, 35, 15, 1, Qt.AlignHCenter
        )
        self.layout.addWidget(self.slider['items']['llim'],
                16, 35, 1, 3
        )

        self.centralWidget = QWidget(self)
        self.centralWidget.setLayout(self.layout)
        '''
        self.layout = QHBoxLayout()
        self.layout.addLayout(self.menu)
        self.layout.addWidget(self.caseplot)
    
        print(self.caseplot)
        self.centralWidget = QWidget(self)
        self.centralWidget.setLayout(self.layout)

        
    """ ===========================
               WIDGET CREATION
        ==========================="""

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
        self.buttons['layout'].setSpacing(0)

        maxbuttons = 8
        cols = [QVBoxLayout() for x in range(int(len(self.buttons['items'])/maxbuttons+1))]
        i = 0
        for key, setup in self.buttons['items'].items():
            self.buttons['items'][key] = QPushButton(setup)
            self.buttons['items'][key].clicked.connect(
                self.__getattribute__(f"plot_{key}"))
            cols[int(i/maxbuttons)].addWidget(
                self.buttons['items'][key], 
            )
            i += 1
        for vbox in cols:
            vbox.addStretch()
        # Add drop-down
        self.buttons['items']['dropdown'] = QComboBox()
        varlist = list(self.buttons['items'].keys())
        self.buttons['items']['dropdown'].addItem("", 0)
        # TODO: Check that can be plotted --> len(shape)
        for vartype in ['centered', 'staggered']:
            if self.inplace:
                for var, path in self.casevars.items():
                    if (vartype in path) and (var not in varlist):
                        if (len(self.get(var).shape) <= 3) and \
                            (len(self.get(var).shape) > 1):
                            self.buttons['items']['dropdown'].addItem(var, i)

        cols[-1].addWidget(self.buttons['items']['dropdown'])
        self.buttons['items']['dropdown'].activated[str].connect(self.plot_dropdown)

        custom = QVBoxLayout()
        custom.addWidget(QLabel("Custom formula"))
        self.buttons['items']['custom'] = MyLineEdit()
        self.buttons['items']['custom'].mousePressEvent = \
            lambda _ : self.buttons['items']['custom'].selectAll()
        self.buttons['items']['custom'].returnPressed.connect(self.plot_custom)
        self.buttons['items']['custom'].setToolTip("Plots custom formula. "+
            "Access var by get('var'). Python syntax applies. Expects output"+
            "of shape ({},{}). Hit return to apply. ".format(self.nx, self.ny))
        custom.addWidget(self.buttons['items']['custom'])
        custom.addStretch()

        self.buttons['layout'].addWidget(self.buttons['title'],0,0,1, len(cols))
        for i in range(len(cols)):
            self.buttons['layout'].addLayout(cols[i], 1, i, 1, 1)
        self.buttons['layout'].addLayout(custom, 2,0,1,len(cols))


    def _createRadios(self):
        self._createIonRadio()
        self._createGasRadio()
        self._createCmapRadio()
        self.radios = {'layout': QGridLayout()}

        self.radios['layout'].addWidget(self.ion_radio['frame'], 0, 0, 1, 1)
        self.radios['layout'].addWidget(self.gas_radio['frame'], 1, 0, 1, 1)

    def _createIonRadio(self):
        self.ion_radio = {
            'frame': QFrame(),
            'layout': QVBoxLayout(),
            'group': QButtonGroup(),
            'title': QLabel(self.title("Ion species")),
            'labels': self.ionarray,
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
        self.ion_radio['items'][self.ionarray[0]].setChecked(True)

    def _createGasRadio(self):
        self.gas_radio = {
            'frame': QFrame(),
            'layout': QVBoxLayout(),
            'group': QButtonGroup(),
            'title': QLabel(self.title("Gas species")),
            'labels': self.gasarray,
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
        self.gas_radio['layout'].addStretch()
        self.gas_radio['frame'].setFrameStyle(QFrame.Panel | QFrame.Raised)
        self.gas_radio['frame'].setLayout(self.gas_radio['layout'])
        self.gas_radio['items'][self.gasarray[0]].setChecked(True)


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



    """ ===========================
                 ACTIONS
        ==========================="""

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
        self.caseplot.toggleAbs()

    def change_varscale(self):
        self.caseplpt.setVarscale(float(self.misc['items']['varscale'].text()))
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


    def cmap_radio_clicked(self):
        self.cmap = self.cmap_radio['group'].checkedButton().text()
        self.verts.set_cmap(self.cmap)
        self.cmap_radio['items']['cmap'].clear()
        self.canvas.draw()

    def ion_radio_clicked(self):
        self.ionspecies = self.ion_radio['group'].checkedId()
        title = self.caseplot.getTitle()
        if self.multispecies == 'ion':
            if " for " in title:
                title = title.split("for")[0] + "for {}".format(
                    self.ionarray[self.ionspecies])
            self.caseplot.setTitle(title)
            self.update_var()
        
    def gas_radio_clicked(self):
        self.gasspecies = self.gas_radio['group'].checkedId()
        if self.multispecies == 'gas':
            title = self.caseplot.getTitle()
            if " for " in title:
                title = title.split("for")[0] + "for {}".format(
                    self.gasarray[self.gasspecies])
            self.caseplot.setTitle(title)
            self.update_var()


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
        self.caseplot.toggleGrid()

    def grid_switch_clicked(self):
        self.caseplot.toggleGrid()

    def separatrix_switch_clicked(self):
        self.caseplot.toggleSeparatrix()

    def vessel_switch_clicked(self):
        self.caseplot.toggleVessel()

    def scale_switch_clicked(self):
        from matplotlib.colors import LogNorm, Normalize
        from numpy import log
        lims = self.caseplot.get_lims()
        if (lims[0] <= 0) and (self.log is False):
            self.raise_message("Negative values in var: masking array!")
        self.caseplot.toggleLog()
        self.log = not self.log
        # TODO: detect log from figure?

    """ ===========================
                PLOT ACTIONS
        ==========================="""
    def plot_driver(self, var, title, multispecies=False):
        from numpy.ma import masked_array
        # Save variable as masked_array to mask for log
        self.var = masked_array(var, mask=False)
        self.multispecies = multispecies
        species=''
        if (var.min() < 0) and self.log:
            self.raise_message("Negative values in var: masking array!")
        if multispecies == 'gas':
            species = self.gasarray[self.gasspecies]
            var = var[1:-1, 1:-1, self.gasspecies]
        elif multispecies == 'ion':
            species = self.ionarray[self.ionspecies]
            var = var[1:-1, 1:-1, self.ionspecies]
        else:
            var = var[1:-1, 1:-1]
        self.caseplot.setTitle(title + \
                f" for {species}"*(multispecies is not False))
        self.caseplot.setVar(var)
        self.buttons['items']['dropdown'].setCurrentIndex(0)
        command = self.buttons['items']['custom'].clear()
 

    def plot_dropdown(self, text):
        # TODO: attempt auto-detecting shape
        self.raise_message(text)
        var = self.get(text)
        if len(var.shape) == 3:
            if var.shape[-1] == len(self.ionarray):
                self.multispecies = 'ion'
                self.caseplot.setTitle(text + " for {}".format(\
                    self.ionarray[self.ionspecies]))
                self.var = var
            elif var.shape[-1] == len(self.gasarray):
                self.multispecies = 'gas'
                self.caseplot.setTitle(text + " for {}".format(\
                    self.gasarray[self.gasspecies]))
                self.var = var
            else:
                self.raise_message("Could not determine species"+
                    f" index of {text}")
                return
        elif len(var.shape) == 2:
            self.var = var
            self.caseplot.setTitle(uptitle)
        else:
            self.raise_message(f"Could not determine shape of {text}")
            return
        if (var.min() < 0) and self.log:
            self.raise_message("Negative values in var: masking array!")

        self.buttons['items']['custom'].clear()
        self.update_var()

    def plot_custom(self):
        try:
            from uedge import com, bbb, grd, flx, aph, api
        except:
            pass
        command = self.buttons['items']['custom'].text()
        self.buttons['items']['dropdown'].setCurrentIndex(0)
        get = self.get
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
                    self.nx+2,
                    self.ny+2
                )
            )
        else:
            self.var = locals()['var']
            for substr in ['get', '(', ')', '"', "'", ".", "bbb", "com",
                "grd", "flx", "aph", "api", "self"]:
                command = command.replace(substr, '')
            self.caseplot.setTitle(command)
            self.buttons['items']['custom'].clearFocus()
            self.multispecies = False
            self.update_var()


       

    def plot_te(self):
        self.plot_driver(
            self.get('te')/1.602e-19,
            'Electron temperature [eV]',
        )

    def plot_ti(self):
        self.plot_driver(
            self.get('ti')/1.602e-19,
            'Ion temperature [eV]',
        )

    def plot_tg(self):
        self.plot_driver(
            self.get('tg')/1.602e-19,
            'Gas temperature [eV]',
            'gas'
        )

    def plot_ne(self):
        self.plot_driver(
            self.get('ne'),
            r'Electron density [m$\mathrm{{}^{-3}}$]',
        )

    def plot_ni(self):
        self.plot_driver(
            self.get('ni'),
            r'Ion density [m$\mathrm{{}^{-3}}$]',
            'ion'
        )

    def plot_ng(self):
        self.plot_driver(
            self.get('ng'),
            r'Gas density [m$\mathrm{{}^{-3}}$]',
            'gas'
        )

    def plot_phi(self):
        self.plot_driver(
            self.get('phi'),
            'Electrical potential [V]',
        )

    def plot_prad(self):
        self.plot_driver(
            self.get('prad')+self.get('pradhyd'),
            'Total radiated power [W/m$\mathrm{{}^{-3}}$]',
        )

    def plot_pradhyd(self):
        self.plot_driver(
            self.get('pradhyd'),
            'Hydrogenic radiated power [W/m$\mathrm{{}^{-3}}$]',
        )

    def plot_pradimp(self):
        self.plot_driver(
            self.get('prad'),
            'Impurity radiated power [W/m$\mathrm{{}^{-3}}$]',
        )

    def plot_psorc(self):
        self.plot_driver(
            self.get('psorc'),
            'Ion ionization source [parts/s]',
            'ion'
        )

    def plot_psorgc(self):
        self.plot_driver(
            self.get('psorgc'),
            'Gas ionization sink [parts/s]',
            'gas'
        )

    def plot_psorxrc(self):
        self.plot_driver(
            self.get('psorxrc'),
            'Ion recombination + CX sink [parts/s]',
            'ion'
        )

    def plot_psorrgc(self):
        self.plot_driver(
            self.get('psorrgc'),
            'Gas recombination source [parts/s]',
            'gas'
        )

    def plot_psorcxg(self):
        self.plot_driver(
            self.get('psorcxg'),
            'Gas CX source [parts/s]',
            'gas'
        )

    """ ===========================
                FORMATTERS 
        ==========================="""
    def title(self, text):
        return f"<b><font size=4>{text}</font></b>"

    def raise_message(self, message, time=5000):
        QApplication.sendEvent(self, QStatusTipEvent(message))
        self.timer = QTimer()
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(self.hide_message)
        self.timer.start(time)

    def hide_message(self):
        QApplication.sendEvent(self, QStatusTipEvent(''))

    def mousePressEvent(self, event):
        focused_widget = QApplication.focusWidget()
        if isinstance(focused_widget, MyLineEdit):
            focused_widget.clearFocus()
        if isinstance(focused_widget, MyClearLineEdit):
            focused_widget.clear()
            focused_widget.clearFocus()
        QMainWindow.mousePressEvent(self, event)


class HeatmapInteractiveFigure(QWidget):
    """Figure Widget with slider."""
    def __init__(self, case, parent=None):
        """Initializer."""
        from numpy import zeros
        super().__init__(parent)

        # TODO: Figure out what is needed to initialize a figure!
        self.canvas = MplCanvas(self, width=5, height=4, dpi=100)
        self.canvas.setMinimumWidth(650)
        self.canvas.setMinimumHeight(750)
        self.case = case        
        self.case.plotmesh(ax=self.canvas.axes, colorbar=True)
        self.case.set_speciesarrays()
        self.verts = self.case.Qvertices
        
        self.canvastoolbar = NavigationToolbar2QT(self.canvas, self)
        self.canvaslayout = QVBoxLayout()
        self.canvaslayout.addWidget(self.canvastoolbar)
        self.canvaslayout.addWidget(self.canvas)
        (ax, cax) = self.canvas.fig.get_axes()
        old_cax = cax.get_position()
        old_ax = ax.get_position()
        dx = 0.055
        cax.set_position([old_cax.x0+dx, old_cax.y0, old_cax.x1+dx, 0.85])
        ax.set_anchor("C")
        ax.set_position([0.1, old_ax.y0, 0.75, 0.85])

        self.casename = self.case.casename
        self.gasspecies = 0
        self.ionspecies = 0
        self.cmap = 'magma'
        self.grid = False
        self.log = False
        self.abs = False
        self.varscale = 1
        self.suptitle = ""
        self.verts.set_cmap(self.cmap)
        self.verts.set_edgecolor('face')
        self.verts.set_linewidths(1)
        self.xy = self.case.get("nx") * self.case.get("ny")
        self.var = zeros((self.case.get("nx"), self.case.get("ny")))
        for line in self.canvas.axes.lines:
            if "child" in line.get_label():
                line.remove()

        self._createSlider()

        self.sliderBar = QHBoxLayout()
        self.sliderBar.addStretch()
        self.sliderBar.addWidget(self.slider['items']['slider'])
        self.sliderBar.addStretch()

        self.sliderControl = QVBoxLayout()
        self.sliderControl.addWidget(self.slider['items']['ulim'])
        self.sliderControl.addLayout(self.sliderBar)
        self.sliderControl.addWidget(self.slider['items']['llim'])

        self.layout = QHBoxLayout()
        self.layout.addLayout(self.canvaslayout)
        self.layout.addLayout(self.sliderControl)

        self.centralWidget = QWidget(self)
        self.centralWidget.setLayout(self.layout)
        self.canvas.draw()

    def setTitle(self, title):
        self.suptitle = title
        self.canvas.fig.suptitle(self.suptitle)
        self.updatePlot()

    def getTitle(self):
        return self.suptitle

    def setVar(self, var):
        self.var = var
        self.updatePlot()

    def getVar(self):
        return self.var

    def setCmap(self, cmap):
        self.cmap = cmap
        vertices.set_cmap(cmap)
        self.updatePlot()

    def toggleAbs(self):
        self.abs = not self.abs
        self.updatePlot()

    def toggleSeparatrix(self):
        for line in self.canvas.axes.lines:
            if line.get_label() == "lcfs":
                line.set_visible(not line.get_visible())
        self.canvas.draw()

    def setVarscale(self, val):
        self.varscale = val
        self.update_var()

    def toggleVessel(self):
        for line in self.canvas.axes.lines:
            if line.get_label() == 'vessel':
                line.set_visible(not line.get_visible())
        self.canvas.draw()

    def toggleLog(self):
        from matplotlib.colors import LogNorm, Normalize
        from numpy import log
        lims = self.get_lims()

        if (lims[0] <= 0) and (self.log is False):
            var = self.verts.get_array()
            if var is not None:
                var.mask=(var<=0)
                self.verts.set_array(var)
                lims = self.get_lims()
        else:
            var = self.verts.get_array()
            if var is not None:
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


    def togglePlates(self):
        for line in self.canvas.axes.lines:
            if 'plate' in line.get_label():
                line.set_visible(not line.get_visible())
        self.canvas.draw()

    def toggleGrid(self):
        if self.grid:
            self.verts.set_edgecolor('face')
            self.verts.set_linewidths(1)
        else:
            self.verts.set_edgecolor('lightgrey')
            self.verts.set_linewidths(0.08)
        self.grid = not self.grid
        self.canvas.draw()


    def updatePlot(self):
        from copy import deepcopy
        from numpy import sum
        self.suptitle = self.suptitle.split("(")[0] \
            + (self.varscale != 1)*" (x{:.3g})".format(self.varscale)
        var = self.varscale*self.var
#        if sum(var) == 0:
#            self.raise_message("Requested variable unpopulated! "+
#                "Select another variable or species to plot.")
#            return
        if self.abs:
            var = abs(var)
#        if (var.min() < 0) and self.log:
#            self.raise_message("Negative values in var: masking array!")
#            var.mask=(var<=0)
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
        self.slider['items']['ulim'].setFixedWidth(75)
        self.slider['items']['llim'].setFixedWidth(75)

        self.slider['items']['slider'].setMinimum(0)
        self.slider['items']['slider'].setMaximum(1000)
        self.slider['items']['slider'].setLow(0)
        self.slider['items']['slider'].setHigh(1000)
        self.slider['items']['slider'].setMinimumHeight(30)
        self.slider['items']['slider'].setValue(1000)
        self.slider['items']['slider'].sliderMoved.connect(\
                self.update_plot_range)


    """ ===========================
               SLIDER HELPERS 
        ==========================="""
    def get_lims(self):
        var = self.verts.get_array()
        if var is None:
            return (0, 1)
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


    def title(self, text):
        return f"<b><font size=4>{text}</font></b>"



class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=2, height=4, dpi=300):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MplCanvas, self).__init__(self.fig)

class MyLineEdit(QLineEdit):
        pass

class MyClearLineEdit(QLineEdit):
        pass


class MainMenu(QMainWindow):
    """Main Window."""
    def __init__(self, parent=None, case=None):
        """Initializer."""
        from os import getcwd
        super().__init__(parent)
        
        self.tablist = []
        self.file = "No file loaded!"
        self.lastpath = getcwd()

        self.setFocus()
        self.setWindowTitle("UETOOLS Main Menu")
        self.resize(1300, 900)

        self._createActions()
        self._createMenuBar()
        self._connectActions()
        self._createStatusBar()

#        self.centralWidget = QTabWidget()#QLabel("Hello, World")
#        self.centralWidget.setMinimumWidth(1300)
#        self.centralWidget.setMinimumHeight(900)
#        self.setCentralWidget(self.centralWidget)

        if case is not None:
            case =  CaseDashboard(Case(file, inplace=True))
            self.tablist.append(self.centralWidget.addTab(case,
                f"{case.casename} heatmap")
            )


    def openFile(self):
        # Logic for opening an existing file goes here...
        # TODO: Move focus to opened file

        if 1==0:
            file = QFileDialog.getOpenFileName(self, 'Open UETOOLS save', 
            self.lastpath, "All files (*.*)")[0]
            self.lastpath = "/".join(file.split("/")[:-1])
        else:
            file = "/Users/holm10/Documents/fusion/uedge/src/"+\
                    "UETOOLS/jupyter/testcase_hires/dashtest.hdf5"
            print("USING TUTORIAL CASE")
        if len(file.strip()) == 0:
            return

        """
        case = Case(file, inplace=True)
        te = case.get('te')[1:-1,1:-1]/1.602e-19
        case =  HeatmapInteractiveFigure(case)
        case.setVar(te)
#        case =  HeatmapInteractiveFigure(case, case.get("te")/1.602e-19)
        self.tablist.append(self.centralWidget.addTab(case, 'test'))
        return
        """

#        case =  CaseDashboard(Case(file, inplace=True))
        case = Case(file, inplace=True)
        case =  HeatmapInteractiveFigure(case)
        self.setCentralWidget(case)#self.centralWidget)
        return
        case.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        # TODO: Figure out how to resize Widget with Window?!
        self.tablist.append(self.centralWidget.addTab(case,
            f"{case.casename} heatmap")
        )

        self.raise_message(f"File > Opened {file}.")


    def _createMenuBar(self):
        menuBar = self.menuBar()
        menuBar.setNativeMenuBar(False)
        # File menu
        fileMenu = QMenu("&File", self)
        menuBar.addMenu(fileMenu)
        fileMenu.addAction(self.openAction)
        fileMenu.addSeparator()
        fileMenu.addAction(self.exitAction)
        editMenu = menuBar.addMenu("&Tab")
        editMenu.addAction(self.renameAction)
        editMenu.addAction(self.closeTabAction)
        helpMenu = menuBar.addMenu("&Help")
        helpMenu.addAction(self.helpContentAction)
        helpMenu.addAction(self.aboutAction)



    def _createActions(self):
        # Creating actions using the second constructor
        self.openAction = QAction("Open heatmap from HDF5 save...", self)
        self.renameAction = QAction("Rename tab", self)
        self.closeTabAction = QAction("Close tab", self)
        self.exitAction = QAction("&Exit", self)
        self.helpContentAction = QAction("&Help Content", self)
        self.aboutAction = QAction("&About", self)
        
    def _createStatusBar(self):
        self.statusbar = self.statusBar()
#        self.wcLabel = QLabel(f"{self.file}")
#        self.statusbar.addPermanentWidget(self.wcLabel)

    def _connectActions(self):
        # Connect File actions
        self.openAction.triggered.connect(self.openFile)
        self.renameAction.triggered.connect(self.editTab)
        self.closeTabAction.triggered.connect(self.closeTab)
        self.exitAction.triggered.connect(self.close)
#        self.openRecentMenu.aboutToShow.connect(self.populateOpenRecent)

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

    def editTab(self):
        # Logic for opening an existing file goes here...

        newname, ok = QInputDialog().getText(self, "Edit tab name",
                "New tab name:", QLineEdit.Normal,
                self.centralWidget.tabText(self.centralWidget.currentIndex()))
        if ok and newname:
            self.centralWidget.setTabText(
                self.centralWidget.currentIndex(), 
                newname
            )


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


    def closeTab(self):
        # TODO Get a handle on the tab naming
        self.centralWidget.removeTab(self.centralWidget.currentIndex())

    def raise_message(self, message, time=5000):
        self.statusbar.showMessage(message, time)

if __name__ == "__main__":
    app = QApplication(sys.argv)

    win = MainMenu()
    win.show()
    win.openFile()
    sys.exit(app.exec_())
else:
    app = QApplication(sys.argv)
    win = CaseDashboard(self)
    
     
    win.show()
    sys.exit(app.exec_())
