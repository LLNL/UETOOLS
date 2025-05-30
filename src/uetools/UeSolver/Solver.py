# Package setting up solvers, time-stepping, etc
# Holm10 Dec 10 2022, based on rundt.py
from uetools.UeUtils import Tools
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
    """Class providing UETOOLS solver routines

    Inherits the uedge.UeRun class

    Methods
    -------
    exmain()
        takes a time step by calling bbb.exmain
    takestep(**kwargs)
        takes a time step by calling self.exmain and returns info
    conv_step(increment, name, var, ivar=None, stop=None,
            inversestep=False, **kwargs)
        converges variable in steps sized by increment time-dependetly
    conv_step_igaso(increment, name, inwsro, **kwargs)
        converges igaso in steps time-dependently
    conv_step_ncore(increment, name, iisp=0, **kwargs)
        converges ncore in steps time-dependently
    conv_step_b0(self, name, increment=0.04, stop=1, **kwargs)
        converges b0 in uniform 1/b0 steps time-dependently
    store_oldgrid()
        stores the old grid locally
    gridmorph(newgrid, savedir, var={}, gridtarget=1, **kwargs)
        morphs the current grid to the new one using continutaion solver
    morphed_mesh(newgrid, fraction, standalone=True)
        sets the current grid to a linear interp. from old to new grid
    populate(silent=True, verbose=None, **kwargs)
        populates all UEDGE arrays without preconditioning


    """

    def __init__(self, case, *args, **kwargs):
        """Initializes the solver and links Case methods required"""
        self.info = case.info
        self.get = case.get
        self.setue = case.setue
        self.getue = case.getue
        self.tools = Tools()
        self.save = case.save
        self.mutex = case.mutex
        self.update = case.update
        # Intialize parent class
        super().__init__(*args, **kwargs)

    def exmain(self):
        """Take a single time step by calling `bbb.exmain()`

        Ensures the time-step is taken at the correct location
        to map paths.

        NOTE: This might be obsolete with recent UETOOLS developments
        """
        if self.info["restored_from_hdf5"] is False:
            original_wd = os.getcwd()
            try:
                # Run in case directory
                os.chdir(self.info["location"])
                packageobject("bbb").getpyobject("exmain")()
            finally:
                os.chdir(original_wd)
        else:
            packageobject("bbb").getpyobject("exmain")()

    def takestep(self, **kwargs):
        """Take a single timestep by calling `self.exmain()`.

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

    def conv_step(
        self, increment, name, var, ivar=None, stop=None, inversestep=False, **kwargs
    ):
        """Converges var in steps of increment size using dt-solver

        **kwargs passed UeRun.converge

        Arguments
        ---------
        increment - float of increments to vary var with (var+increment)
        name - name of the case to be written, appended with var value
        var - str of var to be edited

        Keyword arguments
        -----------------
        ivar : int (default = None)
            index of var to be varied. Applies to whole array if None
        stop : float (default = None)
            value to stop scanning var at. if none, scans to var/5 or
            var*5, depending on direction
        inversestep : bool (default = False)
            if true, takes increments 1/var instead of var

        Returns
        -------
        None
        """
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

    def conv_step_igaso(self, increment, name, inwsor, **kwargs):
        """Increments igaso in steps of increment using dt solver

        Calls and passes **kwargs to Solver.conv_step

        Arguments
        ---------
        increment - float of increments to igaso
        name - name of cases writte, value of var appended
        inwsor - int of index of source to be incremented

        Returns
        -------
        None
        """
        self.conv_step(increment, name, "igaso", ivar=inwsor, **kwargs)

    def conv_step_ncore(self, increment, name, iisp=0, **kwargs):
        """Increments ncore in steps of increment using dt solver

        Calls and passes **kwargs to Solver.conv_step

        Arguments
        ---------
        increment - float of increments to igaso
        name - name of cases writte, value of var appended

        Keyword arguments
        -----------------
        iisp : int (default = 0)
           int of index of ncore species to be incremented

        Returns
        -------
        None
        """
        self.conv_step(increment, name, "ncore", ivar=iisp, **kwargs)

    def conv_step_b0(self, name, increment=0.03, stop=1, **kwargs):
        """Increments 1/b0 in steps of increment using dt solver

        Calls and passes **kwargs to Solver.conv_step

        Arguments
        ---------
        name - name of cases writte, value of var appended

        Keyword arguments
        -----------------
        increment : float (default = 0.03)
            increments to b0
        stop : float (default = 1.)
            b0 value to stop at

        Returns
        -------
        None
        """
        self.conv_step(increment, name, "b0", stop=stop, inversestep=True)

    def store_oldgrid(self):
        """Stores current UEDGE grid parameters to Solver.oldgrid dict"""
        from copy import deepcopy

        self.oldgrid = {}
        for grdvar in ["rm", "zm", "psi", "br", "bz", "bpol", "bphi", 
            'sibdrys', 'simagxs', "b"]:
            self.oldgrid[grdvar] = deepcopy(self.getue(grdvar))

    def gridmorph(self, newgrid, savedir, var={}, gridtarget=1, **kwargs):
        """Uses the continuation solver to morph the UEDGE grid

        Calls morphed_mesh at every increment of continuation_solve,
        passes **kwargs to continuation_solve

        Arguments
        ---------
        newgrid - path to ne ASCII or HDF5 grid (depending on
            com.isgriduehdf5 value)
        savedir - str of folder where to write intermediate saves

        Keyword arguments
        -----------------
        var : dict (default = {})
            dictionary containing additional variables to be evolved
            simulataneously as grid, follows continuation_solve defs.
        gridtarget : float (default = 1)
            how far to evaluate the grid, 1=100%
        """
        from copy import deepcopy

        manualgrid = deepcopy(self.get("manualgrid"))
        self.store_oldgrid()

        self.setue("gridmorph", 0)
        self.setue("manualgrid", 1)
        var["gridmorph"] = {"target": gridtarget}

        self.continuation_solve(
            var,
            savedir=savedir,
            commands=[
                f"self.morphed_mesh('{newgrid}',"
                + "self.getue('gridmorph'), standalone=False)",
            ],
            newgeo=True,
            **kwargs,
        )
        self.setue("manualgrid", manualgrid)

    def morphed_mesh(self, newgrid, fraction, standalone=True):
        """Morphs the physical and magnetic mesh between old and new grids

        Grids must have same number of poloidal and radial grid cells.
        The number of cells in each macro-region should be identical.
        First use interpolation routines to interpolate solution
        to new grid. Requires evaluating Solver.store_oldgrid before
        execution.

        Arguments
        ---------
        newgrid - path to ne ASCII or HDF5 grid (depending on
            com.isgriduehdf5 value)
        fraction - linear interpolation between Solution.oldgrid and
            newgrid

        Keyword arguments
        -----------------
        standalone : bool (default = True)
            restores the default in memory if True, necessary to restore
            grid outside of gridmorph

        """
        from h5py import File
        from Forthon import packageobject
        from copy import deepcopy

        # List of variables to morph
        variables = ["rm", "zm", "psi", "br", "bz", "bpol", "bphi", "b", 
                "sibdrys", "simagxs"]
        # "nlim", "xlim", "ylim", "nplate1", "nplate2", "rplate1", "rplate2",
        # "zplate1", "zplate2"
        if "griddata" not in dir(self):
            self.griddata = {"old": {}, "new": {}, "delta": 0}
            for var in variables:
                self.griddata["old"][var] = deepcopy(self.getue(var))
                self.griddata["new"][var] = self.tools.hdf5search(newgrid, var)
        self.griddata["old"] = self.oldgrid
        for var in variables:
            self.setue(
                var,
                (1 - fraction) * self.griddata["old"][var]
                + fraction * self.griddata["new"][var],
            )
        self.griddata["delta"] = fraction
        packageobject("bbb").__getattribute__("guardc")()
        if standalone:
            iprint = deepcopy(self.get("iprint"))
            self.setue("iprint", 0)
            self.populate(silent=True, verbose=False)
            self.setue("iprint", iprint)
        del self.griddata

    def populate(self, silent=True, verbose=None, **kwargs):
        """Populates all UEDGE arrays by evaluating static 'time-step'

        Outputs prompt assessing whether case is converged or not.

        Keyword arguments
        -----------------
        silent : bool (default = True)
            Tells Case object to silence UEDGE exmain writes.
        verbose : bool (default = None)
            Tells Case to output UETOOLS prompts if True. If
            verbose = None, uses Case.verbose default.

        Modifies
        --------
        UEDGE memory : updates all UEDGE array to correspond to
            restored state

        Returns
        -------
        None
        """
        from copy import deepcopy
        import os

        if self.mutex() is False:
            raise Exception("Case doesn't own UEDGE memory")

        if verbose is None:
            verbose = self.info["verbose"]

        if silent is True:
            old_iprint = self.getue("iprint")
            self.setue("iprint", 0)

        issfon = deepcopy(self.getue("issfon"))
        ftol = deepcopy(self.getue("ftol"))
        try:
            self.setue("issfon", 0)
            self.setue("ftol", 1e20)
            self.exmain()
            self.setue("gengrid", 0)  # Ensure that grids aren't generated again
        finally:
            # Ensure that original settings and working directory are restored
            self.setue("issfon", issfon)
            self.setue("ftol", ftol)

        self.update()  # Reloads variables from UEDGE

        if verbose is True:
            fnrm = sum(self.getue("yldot") ** 2) ** 0.5
            prtstr = "\n*** UEDGE arrays populated: {} ***"
            if fnrm < 10:
                print(prtstr.format("Case appears converged"))
                print("fnrm without preconditioning: {:.2e}\n".format(fnrm))
            elif fnrm < 100:
                print(prtstr.format("Warning, case may noy be fully " "converged"))
                print("fnrm without preconditioning: {:.1f}\n".format(fnrm))
            else:
                print(prtstr.format("WARNING, case NOT converged"))
                print("fnrm without preconditioning: {:.2e}\n".format(fnrm))
        if silent is True:
            self.setue("iprint", old_iprint)  # Restore
