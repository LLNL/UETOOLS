# Package setting up solvers, time-stepping, etc
# Holm10 Dec 10 2022, based on rundt.py
import os
from copy import deepcopy
try:
    from Forthon import packageobject
except:
    pass
try:
    from uedge.rundt import UeRun
except:
    # Dummy class for standalone evaluation
    class UeRun:
        pass


class Solver(UeRun):
    def __init__(self, case, *args, **kwargs):
        self.restored_from_hdf5 = case.info['restored_from_hdf5']
        self.location = case.info['location']

    def exmain(self):
        #        self.exmain_evals = self.get('exmain_evals') + 1
        if self.restored_from_hdf5 is False:
            original_wd = os.getcwd()
            try:
                # Run in case directory
                os.chdir(self.location)
                packageobject("bbb").getpyobject("exmain")()
            finally:
                os.chdir(original_wd)
        else:
            packageobject("bbb").getpyobject("exmain")()
            

    def takestep(self, **kwargs):
        """
        Take a single timestep by calling `bbb.exmain()`.

        Returns
        -------

        A dictionary containing
        - success: bool
          True if the timestep succeeded
        - failed: bool
          True if the timestep failed.
          Note that a step can neither fail nor succeed if iterm == 2.
        - numvar: int
            Number of equations solved per cell
        - neq : int
            Total number of equations
        - dtreal: float
            The time step attempted, in seconds
        - ftol: float
            Tolerance
        - exmain_aborted: bool
            True if exmain aborted
        - iterm : int
            Termination status
        - fnrm : float
            Function norm (RMS)
        - fn_maxabs : float
            Function (yldot * sfscal) maximum absolute value
        - yldot_maxabs : float
            Unscaled time derivative maximum absolute value
        - nfe : int
            Number of function evaluations
        - nni : int
            Number of nonlinear iterations
        - itermx : int
            Maximum number of nonlinear iterations allowed
        """
        from uedge import bbb

        # Modify settings, storing their original values
        original = {}
        for name, value in kwargs.items():
            original[name] = deepcopy(self.getue(name))
            self.setue(name, value)
        try:
            bbb.exmain_aborted = 0  # Reset flag
            self.exmain()
        finally:
            # Gather data on the call

            yldot = bbb.yldot[: bbb.neq - 1]
            sfscal = bbb.sfscal[: bbb.neq - 1]

            result = {
                "success": bbb.iterm == 1,
                "failed": bbb.iterm not in [1, 2],
                "exmain_aborted": bbb.exmain_aborted != 0,
                "fnrm": sum((yldot * sfscal) ** 2) ** 0.5,
                "yldot_maxabs": max(abs(yldot)),
                "fn_maxabs": max(abs(yldot * sfscal)),
                "nfe": int(bbb.nfe[0]),
                "nni": int(bbb.nni[0]),
            }
            for name in ["iterm", "numvar", "neq", "dtreal", "ftol", "rlx", "itermx"]:
                result[name] = getattr(bbb, name)

            # Restore settings
            for name, value in original.items():
                self.setue(name, value)
            # Following steps should reuse state
            self.setue("restart", 1)

        return result

    def conv_step(self, increment, name, var, ivar=None, stop=None, inversestep=False, **kwargs):
        from copy import deepcopy
        from numpy import ndarray
        from uedge import bbb

        if isinstance(self.get(var), ndarray) and (ivar is None):
            print("No index for var specified. Aborting!")
            return
        if stop is None:
            try:
                stop = self.get(var)[ivar]
            except:
                stop = deepcopy(bbb.__getattribute__(var))
            if increment < 0:
                stop /= 5
            else:
                stop *= 5

        try:
            self.get(var)[ivar]
            varstr = "{}[{}]".format(var, ivar)
        except:
            varstr = var
        while True:
            _var = self.getue(var, cp=False)
            if inversestep is True:
                currval = 1 / deepcopy(_var)
                currval += increment
                currval = 1 / currval
                self.setue(var, currval)
            else:
                try:
                    _var[ivar] += increment
                except:
                    _var += increment

            try:
                self.getue(var, copy=False)[ivar] = _var
            except:
                self.setue(var, _var)
            try:
                currval = self.getue(var)[ivar]
            except:
                currval = self.getue(var)
 
            casename = "{}_{}={:.3e}".format(name, var, currval)
            print("======================================")
            print("Solving for {}={:.2E}".format(varstr, currval))
            print("======================================")
            try:
                self.setue("dtreal", kwargs["dtreal"])
            except:
                self.setue("dtreal", 1e-9)
            self.setue("issfon", 1)
            self.setue("ftol", 1e-5)
            self.exmain()
            self.converge(savefname=casename, **kwargs)
            if self.getue("iterm") != 1:
                break
            if (increment > 0) and (currval > stop):
                break
            elif (increment < 0) and (currval < stop):
                break

    def conv_step_igaso(self, increment, name, inwsor, stop=None, **kwargs):
        self.conv_step(increment, name, "igaso", ivar=inwsor, stop=stop, **kwargs)

    def conv_step_ncore(self, increment, name, iisp=0, stop=None, **kwargs):
        self.conv_step(increment, name, "ncore", ivar=iisp, stop=stop, **kwargs)

    def conv_step_b0(self, name, increment=0.03, stop=1, **kwargs):
        self.conv_step(increment, name, "b0", stop=stop, b0=True)
