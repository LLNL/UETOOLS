# Package setting up solvers, time-stepping, etc
# Holm10 Dec 10 2022, based on rundt.py
from Forthon import packageobject
from uedge.rundt import UeRun

class Solver(UeRun):

    def exmain(self):
#        self.exmain_evals = self.get('exmain_evals') + 1
        packageobject('bbb').getpyobject('exmain')()

