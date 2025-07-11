from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
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
try:
    from .range_slider import RangeSlider
except:
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

# TODO: Add "Grid" plot, plotting a B&W grid only
# TODO: Implement diagnostics GUI 
# TODO: Store previously opened files in separate yaml under .uedgerc
# TODO: Implement opening an additional window for secondary parameter

class DatabaseDashboard(QWidget):
    """ Main Window """
    def __init__(self, database, parent=None):
        """Initializer."""
        super().__init__(parent)
        self.ind = 0
        self.db = database
        self.values = self.db.sortvalues
        self.db.sortlabel = self.db.sortlabel.replace(" ", ",")
        self.range = (self.values.min(), self.values.max())

        self.dbtools = QGridLayout()

        self.sortbutton = QPushButton("Sort")
        self.sortbutton.clicked.connect(self.sort)
        
        self.animate = QPushButton("Create animation")
        self.animate.setEnabled(False)

        self.makecase = QPushButton("Open Case in new tab")
        self.makecase.setEnabled(False)

        self.slider = QSlider(Qt.Horizontal)
        self.label = QLabel("")
        self.updateLabel()
        self.label.setFixedWidth(220)

        self.slider.setMinimum(0)
        self.slider.setMaximum(1000)
        self.slider.setValue(0)
        self.slider.valueChanged.connect(\
                self.update_case)
        
        self.buttons = QHBoxLayout()
        self.buttons.addWidget(self.sortbutton)
        self.buttons.addWidget(self.animate)
        self.buttons.addWidget(self.makecase)

        self.sliderbar = QHBoxLayout()
        self.sliderbar.addWidget(self.label)
        self.sliderbar.addWidget(self.slider)

        self.dbtools.addLayout(self.buttons, 0, 0, 1, 1)
        self.dbtools.addLayout(self.sliderbar, 0, 1, 1, 2)

        self.layout = QVBoxLayout()
        self.dash = CaseDashboard(self.db.getcase(0))
        self.case = self.db.getcase(0)
        self.layout.addWidget(self.dash)
        self.layout.addLayout(self.dbtools)
        self.centralWidget = QWidget(self)

        self.setLayout(self.layout)

    def updateLabel(self):
        self.label.setText("<b><font size=4>{} = {:.3g}</font></b>".format(
            self.db.sortlabel, self.values[self.ind]
        ))


    def update_case(self, val):
        if self.range[0] != self.range[1]:
            val = self.range[0] + 0.001*val*(self.range[1]-self.range[0])
            newind = abs(self.values - val).argmin()
        else:
            self.dash.raise_message("Sorting variable constant! Change variable/location.")
            return
        if self.ind != newind:
            self.case = self.db.getcase(self.ind)
            self.ind = newind
            self.updateLabel()
            self.dash.case = self.case
            self.dash.get = self.case.get
            self.dash.caseplot.verts.set_verts(self.case.plot.nodes)
            self.dash.caseplot.updatePlot()
            self.dash.updateVar()
            lines = {}
            for line in self.dash.caseplot.canvas.axes.lines:
                lines[line.get_label()] = line.get_visible()
                line.remove()
            self.case.plot.lcfs(ax=self.dash.caseplot.canvas.axes, flip=True)
            if not lines['lcfs']:
                self.dash.caseplot.toggleSeparatrix()
            self.case.plot.plates(ax=self.dash.caseplot.canvas.axes, flip=True)
            if not lines['plate1']:
                self.dash.caseplot.togglePlates()
            self.case.plot.vessel(ax=self.dash.caseplot.canvas.axes, flip=True)
            if not lines['vessel']:
                self.dash.caseplot.toggleVessel()
            self.dash.caseplot.canvas.draw()

    def sort(self):
        dlg = DatabaseSort(self.db)
        if dlg.exec():
            self.values = self.db.sortvalues
            self.range = (self.values.min(), self.values.max())
            self.ind = -1
            self.slider.setValue(0)
            self.update_case(self.range[0])
            # TODO: Updating stuff goes here!

class DatabaseSort(QDialog):
    def __init__(self, db, parent=None):
        super().__init__(parent)
        self.setFocus()
        self.complete = False
        self.db = db
        self.setWindowTitle("Sort Database")

        self.setMinimumSize(400,500)
        self.layout = QVBoxLayout()
        message = QLabel("Select a variable to sort the Database by:")
        self.layout.addWidget(message)

        self.type = QComboBox()
        self.type.addItem("Predefined", 0)
        self.type.addItem("Variable", 1)
        self.type.addItem("Custom", 2)
        self.type.activated[str].connect(self.sort_var)
        self.layout.addWidget(self.type)

        self.selector = QVBoxLayout()
        self.layout.addLayout(self.selector)
        self.layout.addStretch()

        self.layout.addWidget(QLabel("Scaling factor"))
        self.scaling = MyLineEdit("1")
        self.layout.addWidget(self.scaling)

        self.layout.addWidget(QLabel("Label"))
        self.label = MyLineEdit(self.db.sortlabel)
        self.label.mousePressEvent = lambda _ : self.label.selectAll()
        self.layout.addWidget(self.label)

        QBtn =  QDialogButtonBox.Ok | \
                QDialogButtonBox.Cancel
        self.buttonBox = QDialogButtonBox(QBtn)
        self.buttonBox.accepted.connect(self.sortdb)
        self.buttonBox.rejected.connect(self.reject)

        self.message = QLabel()
        self.layout.addWidget(self.message)    

        self.layout.addWidget(self.buttonBox)
        self.setLayout(self.layout)

        self.omitvar = [
            "equationkey",
        ]
        self.predef_loc = {
            "OMP separatrix": [0, "OMP-sep"],
            "LFS target separatrix": [11, "LFSt-sep"],
            "HFS target separatrix": [12, "HFSt-sep"],
            "LFS target max": [21, "LFSt-max"],
            "HFS target max": [22, "HFSt-max"],
            "LFS target min": [31, "LFSt-min"],
            "HFS target min": [32, "HFSt-max"],
        }
        self.predef_var = {
            "Electron density": ["ne", 1],
            "Electron temperature": ["te", 1/1.602e-19],
            "Ion temperature": ["ti", 1/1.602e-19],
            "Ion density": ["ni", 1],
        }
        self.sort_var("Predefined")


    def sort_var(self, text):
        self.complete=False
        self.setFocus()
        for i in reversed(range(self.selector.count())): 
            try:
                self.selector.itemAt(i).widget().setParent(None)
            except:
                layout_item = self.selector.itemAt(i)
                if isinstance(layout_item, (QHBoxLayout, QVBoxLayout)):
                    for j in reversed(range(layout_item.count())): 
                        layout_item.itemAt(j).widget().setParent(None)
                    self.selector.removeItem(layout_item)
        if text == "Variable":
            self.variable = QComboBox()
            self.vardict = {}
            self.variable.addItem("", 0)
            for var, path in self.db.getcase(0).vars.items():
                try:
                    if (len(var)>0) and ("detected" not in path):
                        shape = [x for x in self.db.getcase(0).get(var).shape if x>1]
                        if (len(shape) > 0) and (var not in self.omitvar): 
                            self.vardict[var] = tuple(shape)
                except:
                    pass

            varlist = list(self.vardict.keys())
            varlist.sort()
            for var in varlist:
                self.variable.addItem(var, 1+len(self.vardict))
            self.variable.activated[str].connect(self.pick_location)
            self.selector.addWidget(self.variable)
        elif text == "Custom":
            self.selector.addWidget(QLabel("Command to evaluate for each Case"))
            self.customvar = MyLineEdit()
            self.selector.addWidget(self.customvar)
            self.complete=True
            message = QLabel("  ixpt1: {}".format(self.db.ixpt1))
            self.selector.addWidget(message)
            message = QLabel("  ixpt2: {}".format(self.db.ixpt2))
            self.selector.addWidget(message)
            message = QLabel("  ixmp: {}".format(self.db.ixmp))
            self.selector.addWidget(message)
            message = QLabel("  iysptrx: {}".format(self.db.iysptrx))
            self.selector.addWidget(message)
        elif text == "Predefined":
            self.selector.addWidget(QLabel("Location"))
            self.location = QComboBox()
            for label, ind in self.predef_loc.items():
                self.location.addItem(label, ind[0])
            self.selector.addWidget(self.location)
            self.location.activated[str].connect(self.predef_update)
            self.variable = QComboBox()
            i = 0
            for label in self.predef_var.keys():
                self.variable.addItem(label, i)
                i += 1
            self.variable.activated[str].connect(self.predef_update)
            self.selector.addWidget(self.variable)

            self.predef_update("")
            self.complete = True
    
            
    def predef_update(self, text):
        for i in reversed(range(3,self.selector.count())): 
            try:
                self.selector.itemAt(i).widget().setParent(None)
            except:
                layout_item = self.selector.itemAt(i)
                if isinstance(layout_item, (QHBoxLayout, QVBoxLayout)):
                    for j in reversed(range(layout_item.count())): 
                        layout_item.itemAt(j).widget().setParent(None)
                    self.selector.removeItem(layout_item)
        self.label.setText(",".join([
            self.predef_loc[self.location.currentText()][1],
            self.predef_var[self.variable.currentText()][0],
        ]))
        self.scaling.setText("{:.3g}".format(self.predef_var[self.variable.currentText()][1]))
        var = self.db.get(self.predef_var[self.variable.currentText()][0])[0]
        if len(var.shape) == 3:
            self.species = QComboBox()
            for i in range(var.shape[2]):
                self.species.addItem(str(i+1))
            self.species.activated[str].connect(self.update_species_label)
            self.selector.addWidget(self.species)
            self.label.setText(self.label.text() +  "[{}]".format(
                self.species.currentText()
            ))

    def update_species_label(self, text):
        text = self.label.text().split('[')[0]
        self.label.setText(text +  "[{}]".format(
            self.species.currentText()
        ))

    def pick_location(self, text):
        self.setFocus()
        self.var = self.variable.currentText()
        for i in reversed(range(1,self.selector.count())): 
            try:
                self.selector.itemAt(i).widget().setParent(None)
            except:
                layout_item = self.selector.itemAt(i)
                if isinstance(layout_item, (QHBoxLayout, QVBoxLayout)):
                    for j in reversed(range(layout_item.count())): 
                        layout_item.itemAt(j).widget().setParent(None)
                    self.selector.removeItem(layout_item)
        self.indices = QHBoxLayout()
        if len(self.vardict[text])==1:
            self.indices.addWidget(QLabel("Select index: "))
        else:
            self.indices.addWidget(QLabel("Select indices: "))
        self.indarr = []
        for i in self.vardict[text]:
            self.indarr.append(QComboBox())
            for j in range(1,i+1):
                self.indarr[-1].addItem(f"{j}", j)
            self.indices.addWidget(self.indarr[-1])
            self.indarr[-1].activated[str].connect(self.set_index)
        self.selector.addLayout(self.indices)
        try:
            self.selector.addWidget(QLabel("Array shape: "+str(self.vardict[text])))
        except:
            return
        message = QLabel("  ixpt1: {}".format(self.db.ixpt1))
        self.selector.addWidget(message)
        message = QLabel("  ixpt2: {}".format(self.db.ixpt2))
        self.selector.addWidget(message)
        message = QLabel("  ixmp: {}".format(self.db.ixmp))
        self.selector.addWidget(message)
        message = QLabel("  iysptrx: {}".format(self.db.iysptrx))
        self.selector.addWidget(message)
        self.set_index("")

    def set_index(self, text):
        var = self.db.getcase(0).get(self.var)
        i = 0
        self.indices = []
        for j in range(len(var.shape)):
            if var.shape[j] == 1:
                self.indices.append(0)
            else:
                self.indices.append(int(self.indarr[i].currentText()))
                i += 1
        self.complete = True

    def sortdb(self):
        from numpy import sum, mod
        if self.complete:
            sortarr = []
            if self.type.currentText() == "Variable":
                for _, case in self.db.cases.items():
                    var = case.get(self.var)
                    for i in self.indices:
                        var = var[i]
                    sortarr.append(var)
            elif self.type.currentText() == "Custom":
                for _, case in self.db.cases.items():
                    ixpt1 = case.get("ixpt1")[0]
                    ixpt2 = case.get("ixpt2")[0]
                    ixmp = case.get("ixmp")
                    iysptrx = case.get("iysptrx")
                    get = case.get
                    try:
                        exec("sortarr.append(" + self.customvar.text() + ")")
                    except:
                        self.raise_message("Expression could not be parsed")
                        return
                    if not isinstance(sortarr[-1], (int, float)):
                        self.raise_message("Expression did not yield int/float")
                        return
            if self.type.currentText() == "Predefined":
                # Determine location
                locind = self.predef_loc[self.location.currentText()][0]
                x = mod(locind, 10)
                y = int(locind/10)
                for _, case in self.db.cases.items():
                    indices = []
                    ixpt1 = case.get("ixpt1")[0]
                    ixpt2 = case.get("ixpt2")[0]
                    ixmp = case.get("ixmp")
                    iysptrx = case.get("iysptrx")
                    [varname, _] = self.predef_var[self.variable.currentText()]
                    var = case.get(varname)
                    if len(var.shape) == 3:
                        species = int(self.species.currentText())
                        var = var[:,:,species-1]

                    if x == 0:
                        var = var[ixmp]
                    elif x == 1:
                        var = var[1]
                    elif x == 2:
                        var = var[-2]
                    if y in [0, 1]:
                        var = var[iysptrx+1]
                    elif y == 2:
                        var = var.max()
                    elif y == 3:
                        var = var.min()

                    sortarr.append(var)
            try:
                exec("locals()['sort'] = [x*float({}) for x in sortarr]".format(
                    self.scaling.text()
                ))
            except:
                self.raise_message("Invalid scaling factor")
                return
            self.db.sort(locals()['sort'])
            self.db.sortlabel = self.label.text()
            self.accept()
        else:
            self.raise_message("Sorting setup is incomplete!")
            return


    def raise_message(self, message, time=5000):
        self.message.setText(message)
        self.timer = QTimer()
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(self.hide_message)
        self.timer.start(time)

    def hide_message(self):
        self.message.setText("")


class CaseDashboard(QWidget):
    """Main Window."""
    def __init__(self, case, parent=None):
        """Initializer."""
        super().__init__(parent)
        self.file = case.info['filename']

        self.centralWidget = QWidget(self)
        self.caseplot = HeatmapInteractiveFigure(case)
        # Cannibalize Case functions
        self.inplace = case.info['inplace']
        self.casevars = case.variables['stored']
        self.get = case.get
        self.nx = case.get('nx')
        self.ny = case.get('ny')
        self.ionarray = case.about.ionarray
        self.gasarray = case.about.gasarray
        self.casename = case.info['casename']
        self.ionspecies = 0
        self.gasspecies = 0
        self.species = ''


        self._createRadios()
        self._createVarButtons()
        self._createSwitches()

        self.plot_te()
        self.caseplot.setTitle("Electron temperature [eV]")

        self.settings = QVBoxLayout()
        self.settings.addWidget(self.switches['frame'])

        self.menu = QHBoxLayout()
        self.menu.addLayout(self.buttons['layout'])
        self.menu.addWidget(self.radios['frame'])
        self.fullmenu = QVBoxLayout()
        self.fullmenu.addLayout(self.menu)
        self.fullmenu.addStretch()
        self.fullmenu.addLayout(self.settings)

        self.layout = QHBoxLayout(self)
        self.layout.addLayout(self.fullmenu)
        self.layout.addWidget(self.caseplot)
    
        self.setLayout(self.layout)

        
    """ ===========================
               WIDGET CREATION
        ==========================="""
    def _createVarButtons(self):
        from numpy import sum, ndarray
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
        self.buttongroup = QButtonGroup()
        for key, setup in self.buttons['items'].items():
            self.buttons['items'][key] = QPushButton(setup)
            self.buttongroup.addButton(self.buttons['items'][key], i)
            self.buttons['items'][key].setCheckable(True)
            self.buttons['items'][key].clicked.connect(
                self.__getattribute__(f"plot_{key}"))
            cols[int(i/maxbuttons)].addWidget(
                self.buttons['items'][key], 
            )
            # TODO: Implement smarter check to deactivat buttons
            buttonvar = self.get(key)
            if buttonvar is not None:
                vis = True
                if len(buttonvar.shape) == 3:
                    varsum = 0
                    for species in range(buttonvar.shape[-1]):
                        varsum += abs(sum(buttonvar[:,:,species]))
                    if varsum == 0:
                        self.buttons['items'][key].setEnabled(False)
                else:
                    if sum(self.get(key)) == 0:
                        self.buttons['items'][key].setEnabled(False)
                i += 1
        self.checkedButton = self.buttongroup.button(0)
        self.checkedButton.setChecked(True)

        # Add drop-down
        self.buttons['items']['dropdown'] = QComboBox()
        varlist = list(self.buttons['items'].keys())
        self.buttons['items']['dropdown'].addItem("", 0)
        # TODO: Check that can be plotted --> len(shape)
        for vartype in ['centered', 'staggered']:
            for var, path in self.casevars.items():
                inpath = True
                if self.inplace:
                    inpath = (vartype in path)
                if inpath and (var not in varlist):
                    getvar = self.get(var)
                    if isinstance(getvar, ndarray):
                        if (len(getvar.shape) <= 3) and \
                            (len(getvar.shape) > 1) and \
                            (abs(getvar).max() > 0) and \
                            (getvar.shape[-1] != 5) and \
                            (getvar.shape[1] not in (1, 2)):
                            self.buttons['items']['dropdown'].addItem(var, i)

        cols[-1].addWidget(self.buttons['items']['dropdown'])
        self.buttons['items']['dropdown'].activated[str].connect(self.plot_dropdown)

        custom = QVBoxLayout()
        custom.addWidget(QLabel("Custom formula"))
        self.buttons['items']['custom'] = MyLineEdit()
        self.buttons['items']['custom'].mousePressEvent = \
            lambda _ : self.buttons['items']['custom'].selectAll()
        self.buttons['items']['custom'].editingFinished.connect(self.plot_custom)
        self.buttons['items']['custom'].setToolTip("Plots custom formula. "+
            "Access var by get('var'). Python syntax applies. Expects output"+
            "of shape ({},{}). Hit return to apply. ".format(self.nx, self.ny))
        custom.addWidget(self.buttons['items']['custom'])

        self.buttons['layout'].addWidget(self.buttons['title'],0,0,1, len(cols))
        for i in range(len(cols)):
            self.buttons['layout'].addLayout(cols[i], 1, i, 1, 1)
        self.buttons['layout'].addLayout(custom, 2,0,1,len(cols))




    def _createRadios(self):
        self.radios = {
            'frame': QFrame(),
            'layout': QVBoxLayout(),
        }

        self.ion_radio = {
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
        ind = 0
        for label in self.ion_radio['labels']:
            self.ion_radio['items'][label] = QRadioButton(label)
            self.ion_radio['layout'].addWidget(\
                    self.ion_radio['items'][label])
            self.ion_radio['group'].addButton(\
                    self.ion_radio['items'][label], ind)
            self.ion_radio['items'][label].clicked.connect(self.ion_radio_clicked)
            ind += 1
        self.ion_radio['items'][self.ionarray[0]].setChecked(True)
        self.ion_radio['layout'].addStretch()
        self.activeIonRadio = self.ion_radio['items'][self.ionarray[0]]

        self.gas_radio = {
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
        
        self.gas_radio['items'][self.gasarray[0]].setChecked(True)
        self.activeGasRadio = self.gas_radio['items'][self.gasarray[0]]

        radios = QHBoxLayout()
        radios.addLayout(self.ion_radio['layout'])
        radios.addLayout(self.gas_radio['layout'])
        label = QLabel(self.title("Species"))
        label.setAlignment(Qt.AlignHCenter)
        self.radios['layout'].addWidget(label)
        self.radios['layout'].addLayout(radios)
        self.radios['frame'].setFrameStyle(QFrame.Panel | QFrame.Raised)
        self.radios['frame'].setLayout(self.radios['layout'])
        self.radios['frame'].setMaximumWidth(150)
        self.radios['frame'].setMaximumHeight(250)


    """ ===========================
                 ACTIONS
        ==========================="""

    def set_cmap(self, cmap):
        self.cmap = cmap
        self.caseplot.setCmap(cmap)

    def toggle_abs(self):
        self.caseplot.toggleAbs()

    def change_varscale(self):
        self.caseplot.setVarscale(float(self.switches['items']['varscale'].text()))
        self.switches['items']['varscale'].clearFocus()

    def ion_radio_clicked(self):
        self.ionspecies = self.ion_radio['group'].checkedId()
        title = self.caseplot.getTitle()
        if not self.is_nonzero():
            self.activeIonRadio.setChecked(True)
            return
        if self.multispecies == 'ion':
            if " for " in title:
                title = title.split("for")[0] + "for {}".format(
                    self.ionarray[self.ionspecies])
            self.caseplot.setTitle(title)
        self.caseplot.setVar(self.get_speciesvar())
        self.activeIonRadio = self.ion_radio["group"].checkedButton()
        
    def gas_radio_clicked(self):
        self.gasspecies = self.gas_radio['group'].checkedId()
        if not self.is_nonzero():
            self.activeGasRadio.setChecked(True)
            return
        if self.multispecies == 'gas':
            title = self.caseplot.getTitle()
            if " for " in title:
                title = title.split("for")[0] + "for {}".format(
                    self.gasarray[self.gasspecies])
            self.caseplot.setTitle(title)
        self.caseplot.setVar(self.get_speciesvar())
        self.activeGasRadio = self.gas_radio["group"].checkedButton()


    def _createSwitches(self):
        from matplotlib.pyplot import colormaps
        self.switches = {
            'frame': QFrame(),
            'layout': QGridLayout(),
            'title': QLabel(self.title("Plot options")),
            'items': {}
        }
        
        grid = QHBoxLayout()
        col = QHBoxLayout()
        grid.addLayout(col)

        self.switches['layout'].addWidget(self.switches['title'], 0, 0, 1, 2)
        self.switches['layout'].setSpacing(2)
        self.switches['frame'].setLayout(self.switches['layout'])
        self.switches['frame'].setMaximumWidth(500)
        self.switches['frame'].setFrameStyle(QFrame.Panel | QFrame.Raised)

        self.switches['items']['Vessel'] = QPushButton("Vessel")
        self.switches['items']['Vessel'].setCheckable(True)
        self.switches['items']['Vessel'].setChecked(True)
        self.switches['items']['Vessel'].clicked.connect(self.vessel_switch_clicked)

        self.switches['items']['Plates'] = QPushButton("Plates")
        self.switches['items']['Plates'].setCheckable(True)
        self.switches['items']['Plates'].setChecked(True)
        self.switches['items']['Plates'].clicked.connect(self.plates_switch_clicked)

        self.switches['items']['Separatrix'] = QPushButton("Separatrix")
        self.switches['items']['Separatrix'].setCheckable(True)
        self.switches['items']['Separatrix'].setChecked(True)
        self.switches['items']['Separatrix'].clicked.connect(self.separatrix_switch_clicked)

        self.switches['items']['Grid'] = QPushButton("Grid")
        self.switches['items']['Grid'].setCheckable(True)
        self.switches['items']['Grid'].clicked.connect(self.grid_switch_clicked)

        self.switches['items']['Logarithmic'] = QPushButton("Logarithmic")
        self.switches['items']['Logarithmic'].setCheckable(True)
        self.switches['items']['Logarithmic'].clicked.connect(self.log_switch_clicked)

        self.switches['items']['Abs'] = QPushButton("Absolute Value")
        self.switches['items']['Abs'].setCheckable(True)
        self.switches['items']['Abs'].clicked.connect(self.toggle_abs)

        self.switches['items']['varscale'] =  MyLineEdit("1")
        self.switches['items']['varscale'].setValidator(QDoubleValidator())
        self.switches['items']['varscale'].editingFinished.connect(self.change_varscale)
        self.switches['items']['varscale'].mousePressEvent = lambda _ : self.switches['items']['varscale'].selectAll()
        self.switches['items']['varscale'].setToolTip("Hit return to apply.")
        
        self.switches['items']['cmap'] = QComboBox()
        i = 0
        cmaps = colormaps()
        cmaps.sort()
        for  colormap in cmaps:
            if colormap == self.caseplot.cmap:
                default = i
            self.switches['items']['cmap'].addItem(colormap, i)
            i+=1
        self.switches['items']["cmap"].setCurrentIndex(default) 
        self.switches['items']['cmap'].activated[str].connect(self.set_cmap)


        labels = [QLabel("Multiplier"), QLabel("Colormap")]
        for label in labels:
            label.setMaximumWidth(120)
            label.setAlignment( Qt.AlignRight | Qt.AlignVCenter )
        self.switches['layout'].addWidget(self.switches['items']['Vessel'],
                1, 0, 1, 1)
        self.switches['layout'].addWidget(self.switches['items']['Plates'],
                2, 0, 1, 1)
        self.switches['layout'].addWidget(self.switches['items']['Separatrix'],
                1, 1, 1, 1)
        self.switches['layout'].addWidget(self.switches['items']['Grid'],
                2, 1, 1, 1)
        self.switches['layout'].addWidget(self.switches['items']['Logarithmic'],
                1, 2, 1, 1)
        self.switches['layout'].addWidget(self.switches['items']['Abs'],
                2, 2, 1, 1)
        self.switches['layout'].addWidget(labels[0],
                1, 3, 1, 1)
        self.switches['layout'].addWidget(labels[1],
                2, 3, 1, 1)
        self.switches['layout'].addWidget(self.switches['items']['varscale'],
                1, 4, 1, 1)
        self.switches['layout'].addWidget(self.switches['items']['cmap'],
                2, 4, 1, 1)

        for switch in ['Vessel', 'Plates', 'Separatrix', 'Grid', 'Logarithmic', 
            'Abs', 'varscale', 'cmap'
        ]:
            self.switches['items'][switch].setMaximumWidth(120)

    def plates_switch_clicked(self):
        self.caseplot.togglePlates()

    def grid_switch_clicked(self):
        self.caseplot.toggleGrid()

    def separatrix_switch_clicked(self):
        self.caseplot.toggleSeparatrix()

    def vessel_switch_clicked(self):
        self.caseplot.toggleVessel()

    def log_switch_clicked(self):
        from matplotlib.colors import LogNorm, Normalize
        from numpy import log
        lims = self.caseplot.get_lims()
        if (lims[0] <= 0) and (self.caseplot.log is False):
            self.raise_message("Negative values in var: masking array!")
        self.caseplot.toggleLog()
        # TODO: detect log from figure?

    """ ===========================
                PLOT ACTIONS
        ==========================="""

    def get_speciesvar(self):
        if self.multispecies == 'gas':
            self.enable_ionradio(False)
            self.enable_gasradio(True)
            self.species = self.gasarray[self.gasspecies]
            return self.var[:, :, self.gasspecies]
        elif self.multispecies == 'ion':
            self.enable_ionradio(True)
            self.enable_gasradio(False)
            self.species = self.ionarray[self.ionspecies]
            return self.var[:, :, self.ionspecies]
        else:
            self.enable_ionradio(False)
            self.enable_gasradio(False)
            return self.var

    def is_nonzero(self):
        if abs(self.get_speciesvar()).max() == 0:
            self.raise_message("Requested variable unpopulated! "+
                "Selecting populated species index.")
            # Find acceptable species index
            if self.multispecies == 'gas':
                self.gasspecies = 0
                while not abs(self.get_speciesvar()).max():
                    self.gasspecies += 1
                    if self.gasspecies > 30:
                        raise ValueError("Out of loop!")
                self.gas_radio['group'].button(self.gasspecies).setChecked(True)
            elif self.multispecies == 'ion':
                self.ionspecies = 0
                while not abs(self.get_speciesvar()).max():
                    self.ionspecies += 1
                    if self.ionspecies > 30:
                        raise ValueError("Out of loop!")
                self.ion_radio['group'].button(self.ionspecies).setChecked(True)
            else:
                raise ValueError("SHOULD NOT BE HERE")
        else:
            self.currentButton = self.buttongroup.checkedButton()
        return True

    def enable_ionradio(self, status):
        for button in self.ion_radio['group'].buttons():
            button.setEnabled(status)
        if status is True:  
            for i in range(len(self.ionarray)):
                if abs(self.var[:,:,i]).max() == 0:
                    self.ion_radio['group'].button(i).setEnabled(False)

    def enable_gasradio(self, status):
        for button in self.gas_radio['group'].buttons():
            button.setEnabled(status)
        if status is True:  
            for i in range(len(self.gasarray)):
                if abs(self.var[:,:,i]).max() == 0:
                    self.gas_radio['group'].button(i).setEnabled(False)
        

    def uncheck_buttons(self):
        self.buttongroup.setExclusive(False)
        for button in self.buttongroup.buttons():
            button.setChecked(False)
        self.buttongroup.setExclusive(True)

    def plot_driver(self, var, title, multispecies=False):
        from numpy.ma import masked_array
        from numpy import sum
        # Save variable as masked_array to mask for log
        self.var = masked_array(var, mask=False)
        self.multispecies = multispecies
        if not self.is_nonzero():
            return
        if (self.get_speciesvar().min() <= 0) and self.caseplot.log:
            self.raise_message("Non-positive values in logarithmic plot "+
                "have been masked!")
        self.caseplot.setTitle(title + \
                f" for {self.species}"*(multispecies is not False))
        self.caseplot.setVar(self.get_speciesvar())
        self.buttons['items']['dropdown'].setCurrentIndex(0)
        command = self.buttons['items']['custom'].clear()
 
    def updateVar(self):
        if self.buttons['items']['dropdown'].currentIndex() == 0:
            self.lastfunc()
        else:
            self.lastfunc(self.buttons['items']['dropdown'].itemText(
                    self.buttons['items']['dropdown'].currentIndex()
            ))

    def plot_dropdown(self, text):
        # TODO: attempt auto-detecting shape
        if text == "":
            return
        self.raise_message(text)
        var = self.get(text)
        if len(var.shape) == 3:
            if var.shape[-1] == len(self.ionarray):
                self.multispecies = 'ion'
                self.var = var
                if not self.is_nonzero():
                    return
                self.caseplot.setTitle(text + " for {}".format(\
                    self.ionarray[self.ionspecies]))
                self.caseplot.setVar(self.get_speciesvar())
                self.enable_ionradio(True)
                self.enable_gasradio(False)
            elif var.shape[-1] == len(self.gasarray):
                self.multispecies = 'gas'
                self.var = var
                if not self.is_nonzero():
                    return
                self.caseplot.setTitle(text + " for {}".format(\
                    self.gasarray[self.gasspecies]))
                self.caseplot.setVar(self.get_speciesvar())
                self.enable_ionradio(False)
                self.enable_gasradio(True)
            else:
                self.raise_message("Could not determine species"+
                    f" index of {text}")
                return
        elif len(var.shape) == 2:
            self.multispecies = False
            self.var = var
            self.caseplot.setVar(self.var)
            if not self.is_nonzero():
                return
            self.caseplot.setTitle(text)
            self.enable_ionradio(False)
            self.enable_gasradio(False)
        else:
            self.raise_message(f"Could not determine shape of {text}")
            return
        self.buttons['items']['custom'].clear()
        self.lastfunc = self.plot_dropdown
        self.uncheck_buttons()

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
            self.buttons['items']['custom'].setFocus()
            self.buttons['items']['custom'].selectAll()
            return
        if len(locals()['var'].shape) != 2:
            self.raise_message("Command variable has wrong shape: "+
                "{}, expected ({},{})".format(\
                    locals()['var'].shape,
                    self.nx+2,
                    self.ny+2
                )
            )
            self.buttons['items']['custom'].setFocus()
            self.buttons['items']['custom'].selectAll()
            return
        else:
            self.var = locals()['var']
            self.multispecies = False
            if not self.is_nonzero():
                self.buttons['items']['custom'].setFocus()
                self.buttons['items']['custom'].selectAll()
                return
            self.caseplot.setVar(self.var)
            for substr in ['get', '("', '")', "('", "')",
                '"', "'", ".", "bbb", "com",
                "grd", "flx", "aph", "api", "self"]:
                command = command.replace(substr, '')
            self.caseplot.setTitle(command)
            self.buttons['items']['custom'].clearFocus()
        self.lastfunc = self.plot_custom
        self.enable_gasradio(False)
        self.enable_ionradio(False)
        self.uncheck_buttons()

    def plot_te(self):
        self.plot_driver(
            self.get('te')/1.602e-19,
            'Electron temperature [eV]',
        )
        self.lastfunc = self.plot_te

    def plot_ti(self):
        self.plot_driver(
            self.get('ti')/1.602e-19,
            'Ion temperature [eV]',
        )
        self.lastfunc = self.plot_ti

    def plot_tg(self):
        self.plot_driver(
            self.get('tg')/1.602e-19,
            'Gas temperature [eV]',
            'gas'
        )
        self.lastfunc = self.plot_tg

    def plot_ne(self):
        self.plot_driver(
            self.get('ne'),
            r'Electron density [m$\mathrm{{}^{-3}}$]',
        )
        self.lastfunc = self.plot_ne

    def plot_ni(self):
        self.plot_driver(
            self.get('ni'),
            r'Ion density [m$\mathrm{{}^{-3}}$]',
            'ion'
        )
        self.lastfunc = self.plot_ni

    def plot_ng(self):
        self.plot_driver(
            self.get('ng'),
            r'Gas density [m$\mathrm{{}^{-3}}$]',
            'gas'
        )
        self.lastfunc = self.plot_ng

    def plot_phi(self):
        self.plot_driver(
            self.get('phi'),
            'Electrical potential [V]',
        )
        self.lastfunc = self.plot_phi

    def plot_prad(self):
        self.plot_driver(
            self.get('prad')+self.get('pradhyd'),
            'Total radiated power [W/m$\mathrm{{}^{-3}}$]',
        )
        self.lastfunc = self.plot_prad

    def plot_pradhyd(self):
        self.plot_driver(
            self.get('pradhyd'),
            r'Hydrogenic radiated power [W/m$\mathrm{{}^{-3}}$]',
        )
        self.lastfunc = self.plot_pradhyd

    def plot_pradimp(self):
        self.plot_driver(
            self.get('prad'),
            r'Impurity radiated power [W/m$\mathrm{{}^{-3}}$]',
        )
        self.lastfunc = self.plot_pradimp

    def plot_psorc(self):
        self.plot_driver(
            self.get('psorc'),
            'Ion ionization source [parts/s]',
            'ion'
        )
        self.lastfunc = self.plot_psorc

    def plot_psorgc(self):
        self.plot_driver(
            self.get('psorgc'),
            'Gas ionization sink [parts/s]',
            'gas'
        )
        self.lastfunc = self.plot_psorgc

    def plot_psorxrc(self):
        self.plot_driver(
            self.get('psorxrc'),
            'Ion recombination + CX sink [parts/s]',
            'ion'
        )
        self.lastfunc = self.plot_psorxrc

    def plot_psorrgc(self):
        self.plot_driver(
            self.get('psorrgc'),
            'Gas recombination source [parts/s]',
            'gas'
        )
        self.lastfunc = self.plot_psorrgc

    def plot_psorcxg(self):
        self.plot_driver(
            self.get('psorcxg'),
            'Gas CX source [parts/s]',
            'gas'
        )
        self.lastfunc = self.plot_psorcxg

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
        self.canvas = WidgetPlot(self, width=3, height=4, dpi=100)
        self.canvas.setMinimumWidth(550)
        self.canvas.setMinimumHeight(700)
        self.case = case        
        self.case.plot.te2D(ax=self.canvas.axes)
        self.case.about.set_speciesarrays()
        self.verts = self.case.plot.Qvertices
        
        self._createSlider()
        # TODO: moce to MplCanvas??
        self.sliderControl = QGridLayout()
        self.sliderControl.addWidget(self.slider['items']['ulim'],
            0,0,1,3)
        self.sliderControl.addWidget(self.slider['items']['slider'],
            1,1,1,1)
        self.sliderControl.addWidget(self.slider['items']['llim'],
            2,0,1,3)
        self.canvaslayout = QHBoxLayout(self)
        self.canvaslayout.addWidget(self.canvas)
        self.canvaslayout.addLayout(self.sliderControl)

        (ax, cax) = self.canvas.fig.get_axes()
        old_cax = cax.get_position()
        old_ax = ax.get_position()
        dx = 0.03
        cax.set_position([old_cax.x0+dx, old_cax.y0, old_cax.x1+dx, 0.85])
        ax.set_anchor("C")
        ax.set_position([0.1, old_ax.y0, 0.75, 0.85])

        self.casename = self.case.info['casename']
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



        self.centralWidget = QWidget(self)
        self.canvas.draw()

    def _createSlider(self):
        lims = self.get_lims()
        self.slider = {
            'layout': QVBoxLayout(),
            'items': {
                'ulim': QLineEdit(),
                'slider': RangeSlider(Qt.Vertical),
                'llim': QLineEdit(),
            }
        }

        font = self.slider['items']['ulim'].font()
        font.setPointSize(18)               # change it's size
        font.setBold(True)            # change it's size

        self.slider['items']['ulim'].setFont(font) 
        self.slider['items']['ulim'].setText(
            "{:.3g}".format(lims[1])
        )


        self.slider['items']['ulim'].setValidator(QDoubleValidator())
        self.slider['items']['ulim'].editingFinished.connect(
            self.changeUpperSlider
        )
        

        self.slider['items']['llim'].setValidator(QDoubleValidator())
        self.slider['items']['llim'].editingFinished.connect(
            self.changeLowerSlider
        )

        self.slider['items']['llim'].setFont(font) 
        self.slider['items']['llim'].setText(
            "{:.3g}".format(lims[0])
        )

        self.slider['items']['llim'].setAlignment(\
                        Qt.AlignHCenter | \
                        Qt.AlignVCenter
        )
        self.slider['items']['ulim'].setAlignment(\
                        Qt.AlignHCenter | \
                        Qt.AlignVCenter
        )
        self.slider['items']['ulim'].setFixedWidth(90)
        self.slider['items']['llim'].setFixedWidth(90)

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
        self.slider['items']['llim'].setText("{:.3g}".format(nl))
        self.slider['items']['ulim'].setText("{:.3g}".format(nu))
        self.canvas.draw()



    def setTitle(self, title):
        self.suptitle = title
        self.canvas.fig.suptitle(self.suptitle)
        self.canvas.draw()

    def getTitle(self):
        return self.suptitle

    def setVar(self, var):
        from numpy.ma import masked_array
        self.var = masked_array(var)
        self.updatePlot()

    def getVar(self):
        return self.var

    def setCmap(self, cmap):
        self.cmap = cmap
        self.verts.set_cmap(cmap)
        self.canvas.draw()

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
        self.updatePlot()

    def toggleVessel(self):
        for line in self.canvas.axes.lines:
            if line.get_label() == 'vessel':
                line.set_visible(not line.get_visible())
        self.canvas.draw()

    def toggleLog(self):
        from matplotlib.colors import LogNorm, Normalize
        from numpy.ma import is_masked
        from numpy import log
        lims = self.get_lims()

        if (lims[0] <= 0) and (self.log is False):
            var = self.verts.get_array()
            if var is not None:
                var.mask=(var<=0)
                self.verts.set_array(var)
                lims = list(self.get_lims())
                if is_masked(lims[0]):
                    lims[0]=1e-100
                if is_masked(lims[1]):
                    lims[1]=1e100
        else:
            var = self.verts.get_array()
            if var is not None:
                var.mask=False
                self.verts.set_array(var)
                lims = self.get_lims()
        if self.log:
            self.verts.set_norm(Normalize(*lims))#vmin=lims[0], vmax=lims[1])
        else:
            self.verts.set_norm(LogNorm(*lims))
        vals = self.get_slider_values()
        self.log = not self.log 
        newpos = [int(self.value2position(x)) for x in vals]
        if newpos[0] == newpos[1]:
            newpos = [0, 1000]
        self.slider['items']['slider'].setLow(newpos[0])
        self.slider['items']['slider'].setHigh(newpos[1])
        clim = list(self.verts.get_clim())
        if self.log:
            if clim[0] <= 0:
                clim[0]=1e-100
            if (clim[1] <= 0) or (clim[1]<clim[0]):
                clim[1]=1e100
        self.updatePlot()
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


        if self.log:
            if self.abs:
                self.var.mask = False
            else:
                self.var.mask=(self.var<=0)
        else:
            self.var.mask=False
        self.suptitle = self.suptitle.split("(")[0] \
            + (self.varscale != 1)*" (x{:.3g})".format(self.varscale)
        var = self.varscale*self.var
        if sum(var) == 0:
            self.raise_message("Requested variable unpopulated! "+
                "Select another variable or species to plot.")
            return
        if self.abs:
            var = abs(var)
        # Record old values slider values
        self.verts.set_array(var[1:-1,1:-1].reshape(self.xy))
        [llim, ulim] = deepcopy(self.verts.get_clim())
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
        self.verts.set_clim((llim, ulim))
        self.slider['items']['ulim'].setText(\
                "{:.3g}".format(ulim)
        )
        self.slider['items']['llim'].setText(\
                "{:.3g}".format(llim)
        )
        self.slider['items']['slider'].setHigh(
            int( self.value2position(ulim)
        ))
        self.slider['items']['slider'].setLow(
            int( self.value2position(llim)
        ))
        self.verts.set_cmap(self.verts.get_cmap())
        self.canvas.fig.suptitle(self.suptitle)
        self.canvas.draw()

    def changeUpperSlider(self):
        if self.setUpperSlider(float(self.slider['items']['ulim'].text())):
            self.slider['items']['ulim'].clearFocus()
        else:
            self.slider['items']['ulim'].setFocus()
            self.slider['items']['ulim'].selectAll()
        
    def changeLowerSlider(self):
        if self.setLowerSlider(float(self.slider['items']['llim'].text())):
            self.slider['items']['llim'].clearFocus()
        else:
            self.slider['items']['llim'].setFocus()
            self.slider['items']['llim'].selectAll()

    def setUpperSlider(self, ulim):
        clim = self.verts.get_clim()
        lims = self.get_lims()
        if ulim <= clim[0]:
            self.raise_message("Requested upper range {:.3g} is".format(ulim)+
                " below the current lower limit {:.3g}".format(clim[0]))
            return False
        ulim = min(ulim, lims[1])
        self.verts.set_clim((clim[0], ulim))
        self.slider['items']['ulim'].setText(\
                "{:.3g}".format(ulim)
        )
        self.slider['items']['slider'].setHigh(
           int( self.value2position(ulim)
        ))
        self.canvas.draw()
        return True

    def setLowerSlider(self, llim):
        clim = self.verts.get_clim()
        lims = self.get_lims()
        if llim >= clim[1]:
            self.raise_message("Requested lower range {:.3g} is".format(llim)+
                " above the current upper limit {:.3g}".format(clim[1]))
            return False
        llim = max(llim, lims[0])
        self.verts.set_clim((llim, clim[1]))
        self.slider['items']['llim'].setText(\
                "{:.3g}".format(llim)
        )
        self.slider['items']['slider'].setLow(
           int( self.value2position(llim)
        ))
        self.canvas.draw()
        return True


    def raise_message(self, message, time=5000):
        QApplication.sendEvent(self, QStatusTipEvent(message))
        self.timer = QTimer()
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(self.hide_message)
        self.timer.start(time)

    def hide_message(self):
        QApplication.sendEvent(self, QStatusTipEvent(''))

class WidgetPlot(QWidget):
    def __init__(self, *args, width=2, height=4, dpi=300, parent=None, **kwargs):
        from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
        super(WidgetPlot, self).__init__(parent)
        
        self.setLayout(QVBoxLayout())
        self.canvas = MplCanvas(self, width=width, height=height, dpi=dpi)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        self.layout().addWidget(self.toolbar)
        self.layout().addWidget(self.canvas)
        self.axes = self.canvas.axes
        self.fig = self.canvas.fig
        self.draw = self.canvas.draw

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=2, height=4, dpi=300):
        from matplotlib.figure import Figure
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MplCanvas, self).__init__(self.fig)

class MyLineEdit(QLineEdit):
        pass

class MyClearLineEdit(QLineEdit):
        pass
 
class PopUp(QMainWindow):
    def __init__(self, parent=None):
        super(PopUp, self).__init__(parent)
 
 
class StandaloneDatabaseDashboard(QMainWindow):
    def __init__(self, case, parent = None):
        super().__init__(parent)
        
        self.case = case
        self.setFocus()
        self.resize(1250, 850)
    
        self.statusbar = self.statusBar()
        self.centralWidget =  DatabaseDashboard(self.case)
        self.centralWidget.setMinimumWidth(1250)
        self.centralWidget.setMinimumHeight(850)
        self.setCentralWidget(self.centralWidget)


class StandaloneDashboard(QMainWindow):
    def __init__(self, case, parent = None):
        super().__init__(parent)
        
        self.case = case
        self.setFocus()
        self.resize(1250, 850)
    
        self.statusbar = self.statusBar()
        self.centralWidget =  CaseDashboard(self.case)
        self.centralWidget.setMinimumWidth(1250)
        self.centralWidget.setMinimumHeight(850)
        self.setCentralWidget(self.centralWidget)

class MainMenu(QMainWindow):
    """Main Window."""

    # TODO: pop-out tab option

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
        self.popups = []

        self._createActions()
        self._createMenuBar()
        self._connectActions()
        self._createStatusBar()
        

        self.centralWidget = QTabWidget(self)#QLabel("Hello, World")
        self.centralWidget.setMinimumWidth(1300)
        self.centralWidget.setMinimumHeight(800)
        self.centralWidget.setTabsClosable(True)
        self.centralWidget.tabCloseRequested.connect(self.closeTab)
        self.setCentralWidget(self.centralWidget)
        self.tabMenu.setDisabled(True)



    def closeTab(self, index):
        self.centralWidget.removeTab(index)
        if self.centralWidget.count() == 0:
            self.tabMenu.setDisabled(True)


    def openHDF5(self):#, caseobj=None):
        from uetools import Case
        # Logic for opening an existing file goes here...
        if 1==0:
            file = QFileDialog.getOpenFileName(self, 'Open UETOOLS save', 
            self.lastpath, "All files (*.*)")[0]
            self.lastpath = "/".join(file.split("/")[:-1])
        else:
            file = "/Users/holm10/Documents/fusion/uedge/src/"+\
                    "UETOOLS/dashboard_test/testcase_hires/nc20.hdf5"
            print("USING TUTORIAL CASE")
        if len(file.strip()) == 0:
            return
        case =  CaseDashboard(Case(file, inplace=True))
        try:
            case =  CaseDashboard(Case(file, inplace=True))
        except:
            self.raise_message(f"File {file} is not a valid UETOOLS save file!")
            return
        case.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        # TODO: Figure out how to resize Widget with Window?!
        self.tablist.append(self.centralWidget.addTab(case,
            f"{case.casename} heatmap")
        )
        self.centralWidget.setTabToolTip(self.tablist[-1], file)
        self.centralWidget.setCurrentIndex(self.tablist[-1])
        self.raise_message(f"File > Opened HDF5 file {file}.")
        self.tabMenu.setDisabled(False)

    def openDatabase(self):
        from uetools import Database
        if 1==0:
            path = "/Users/holm10/Documents/fusion/uedge/src/UETOOLS/dashboard_test/testcase_hires"
            self.tablist.append(self.centralWidget.addTab(
                DatabaseDashboard(Database(path)), f"test DB")
            )
            print("USING TEST DATABASE")
        else:
            path = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
            if len(path) == 0:
                return
            name = path.split('/')[-1]
            try:
                self.tablist.append(self.centralWidget.addTab(
                    DatabaseDashboard(Database(path)), f"{name} DB heatmap")
                )
            except:
                self.raise_message(f"Database could not be created from {path}!")
                return

        self.centralWidget.setTabToolTip(self.tablist[-1], path)
        self.centralWidget.setCurrentIndex(self.tablist[-1])
        self.raise_message(f"File > Opened Database at {path}.")
        self.tabMenu.setDisabled(False)
#        self.centralWidget.widget(self.tablist[-1]).sort()
            


    def _createMenuBar(self):
        menuBar = self.menuBar()
        menuBar.setNativeMenuBar(False)
        # File menu
        fileMenu = QMenu("&File", self)
        menuBar.addMenu(fileMenu)
        fileMenu.addAction(self.openHDF5Action)
        fileMenu.addAction(self.openDBAction)
        fileMenu.addSeparator()
        fileMenu.addAction(self.exitAction)
        self.tabMenu = menuBar.addMenu("&Tab")
        self.tabMenu.addAction(self.renameTabAction)
        self.tabMenu.addAction(self.popOutTabAction)
#        helpMenu = menuBar.addMenu("&Help")
#        helpMenu.addAction(self.helpContentAction)
#        helpMenu.addAction(self.aboutAction)



    def _createActions(self):
        # Creating actions using the second constructor
        self.openHDF5Action = QAction("Open heatmap from HDF5 save...", self)
        self.openDBAction = QAction("Create database from folder...", self)
        self.renameTabAction = QAction("Rename tab", self)
        self.popOutTabAction = QAction("Pop out figure", self)
        self.exitAction = QAction("&Exit", self)
        self.helpContentAction = QAction("&Help Content", self)
        self.aboutAction = QAction("&About", self)
        
    def _createStatusBar(self):
        self.statusbar = self.statusBar()

    def _connectActions(self):
        # Connect File actions
        self.openHDF5Action.triggered.connect(self.openHDF5)
        self.openDBAction.triggered.connect(self.openDatabase)
        self.renameTabAction.triggered.connect(self.renameTab)
        self.popOutTabAction.triggered.connect(self.popTab)
        self.exitAction.triggered.connect(self.close)

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

    def renameTab(self):
        # Logic for opening an existing file goes here...

        newname, ok = QInputDialog().getText(self, "Edit tab name",
                "New tab name:", QLineEdit.Normal,
                self.centralWidget.tabText(self.activeTab()))
        if ok and newname:
            self.centralWidget.setTabText(
                self.centralWidget.currentIndex(), 
                newname
            )
    
    def activeTab(self):
        return self.centralWidget.currentIndex()

    def transfer_HeatmapInteractiveFigure(self, plot):
        newplot=HeatmapInteractiveFigure(plot.case)
        lines = {}
        for line in plot.canvas.axes.lines:
            lines[line.get_label()] = line.get_visible()
        for line in newplot.canvas.axes.lines:
            line.set_visible(lines[line.get_label()])
        if plot.grid:
            newplot.toggleGrid()
        if plot.log:
            newplot.toggleLog()
            var = newplot.verts.get_array()
            var.masked = (var <= 0)
            newplot.verts.set_array(var)
        clim = plot.verts.get_clim()
        newplot.setTitle(plot.suptitle)
        newplot.verts.set_clim(clim)
        newplot.canvas.axes.set_xlim(plot.canvas.axes.get_xlim())
        newplot.canvas.axes.set_ylim(plot.canvas.axes.get_ylim())
        newplot.verts.set_array(plot.verts.get_array())
        newplot.setUpperSlider(clim[1])
        newplot.setLowerSlider(clim[0])
        newplot.verts.set_cmap(plot.verts.get_cmap())
        newplot.canvas.draw()
        return newplot


    def popTab(self):
        widget = self.centralWidget.widget(self.activeTab())
        plot = widget.caseplot
        newplot = self.__getattribute__(
            "transfer_{}".format(type(plot).__name__))(plot)
        self.popups.append(PopUp(self))
        self.popups[-1].show()
        title = plot.suptitle
        for substr in ["$", r"\mathrm", "{","}"]:
            title = title.replace(substr, "")
        title = self.centralWidget.tabText(self.activeTab())
        self.popups[-1].setWindowTitle(title)
        self.popups[-1].setCentralWidget(newplot)



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

    def raise_message(self, message, time=5000):
        self.statusbar.showMessage(message, time)

def uedashboard():
    import sys
    import matplotlib
    matplotlib.use('Qt5Agg')

    app = QApplication(sys.argv)
    win = MainMenu()
    win.show()
    sys.exit(app.exec_())
    


if __name__ == "__main__":
    import sys
    import matplotlib
    matplotlib.use('Qt5Agg')
    app = QApplication(sys.argv)
    win = MainMenu()
    win.openHDF5()
    win.show()
    sys.exit(app.exec_())
#else:
#    app = QApplication(sys.argv)
#    win = CaseDashboard(self)    
#    win.show()
#    sys.exit(app.exec_())
