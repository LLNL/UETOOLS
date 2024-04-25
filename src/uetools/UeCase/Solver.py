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
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

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

