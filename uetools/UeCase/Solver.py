# Package setting up solvers, time-stepping, etc
# Holm10 Dec 10 2022, based on rundt.py
from Forthon import packageobject
from uedge.rundt import UeRun
import os
from copy import deepcopy


class Solver(UeRun):
    def exmain(self):
        #        self.exmain_evals = self.get('exmain_evals') + 1
        original_wd = os.getcwd()
        try:
            # Run in case directory
            os.chdir(self.location)
            packageobject("bbb").getpyobject("exmain")()
        finally:
            os.chdir(original_wd)

    def takestep(self, **kwargs):
        # Modify settings, storing their original values
        original = {}
        for name, value in kwargs.items():
            original[name] = deepcopy(self.getue(name))
            self.setue(name, value)
        try:
            self.exmain()
        finally:
            # Restore settings
            for name, value in original.items():
                self.setue(name, value)
            # Following steps should reuse state
            self.setue("restart", 1)

    def converge(self, *args, **kwargs):
        original_wd = os.getcwd()
        try:
            # Run in case directory
            os.chdir(self.location)
            UeRun.converge(self, *args, **kwargs)
        finally:
            # Restore original directory
            os.chdir(original_wd)

    # TODO: Fix this so that methods can be inherited directly!
    def convergenceanalysis(self, filename, **kwargs):
        return UeRun.convergenceanalysis(filename)

    def failureanalysis(self, filename, **kwargs):
        return UeRun.failureanalysis(filename)
