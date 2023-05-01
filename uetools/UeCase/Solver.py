# Package setting up solvers, time-stepping, etc
# Holm10 Dec 10 2022, based on rundt.py
from uedge.rundt import UeRun
from Forthon import packageobject

class Solver():
#    def __init__():
#        super(Solver, self).__init__()
#        return
    def converge(self, **kwargs):
        from uedge.rundt import rundt
        rundt(**kwargs) 

    def exmain(self):
        packageobject('bbb').getpyobject('exmain')()

