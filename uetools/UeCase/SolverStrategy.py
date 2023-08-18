from math import exp


class SolverStrategy:
    """
    Implements a solver strategy, using diagnostics from Solver.takestep
    to decide how to vary the next timestep to take.
    """

    # A list of potentially good settings. This is likely problem-dependent.
    good_settings = [
        [
            0.5,  # On failure, divide timestep by 3.5
            0.5,  # On success, multiply timestep by 3.5
            0.5,  # Sensitivity of nonlinear iteration dt factor
            0.5,  # Offset, setting expected number of iterations
            0.0,  # Importance of yldot dt factor on success
            0.0,  # Importance of nni dt factor on success
            1.0,  # Max -> Min weighting on success
            0.0,  # Importance of yldot dt factor on failure
            0.0,  # Importance of nni dt factor on failure
            1.0,  # Max -> Min weighting on failure
        ],
        [
            0.0786295,
            0.03415136,
            0.85525621,
            0.13271134,
            0.08776637,
            0.33095625,
            0.44373128,
            0.93586442,
            0.03066646,
            0.73827909,
        ],
        [
            0.23123474,
            0.86637406,
            0.69858623,
            0.09671805,
            0.11583457,
            0.20482517,
            0.5334955,
            0.84187559,
            0.3483742,
            0.05839172,
        ],
        [
            0.41424802,
            0.7417847,
            0.52739859,
            0.14051955,
            0.18287689,
            0.44381368,
            0.3586923,
            0.92812855,
            0.17109219,
            0.09163184,
        ],
    ]

    def __init__(self, settings):
        """
        Initialize with settings that decide the solver strategy

        Arguments
        ---------

        settings: list or array
            A 1D array of 10 floating point numbers between 0 and 1

        The optimum settings are likely to be problem-dependent.
        Some settings that may be good starting points are in
        SolverStrategy.good_settings

        Example
        -------

        import uetools
        strategy = uetools.UeCase.SolverStrategy(
            uetools.UeCase.SolverStrategy.good_settings[0]
        )


        """

        assert len(settings) == 10
        assert min(settings) >= 0.0
        assert max(settings) <= 1.0

        self.settings = settings

        # Settings are expected to be between 0 and 1
        self.dtdiv_failed = 1.0 + settings[0] * 5.0
        self.dtmul_succeeded = 1.0 + settings[1] * 5.0
        self.nni_sensitivity = settings[2]
        self.nni_offset = settings[3] * 20

        self.yldot_factor_succeeded = settings[4]
        self.nni_factor_succeeded = settings[5]
        self.dtmaxmin_succeeded = settings[6]

        self.yldot_factor_failed = settings[7]
        self.nni_factor_failed = settings[8]
        self.dtmaxmin_failed = settings[9]

    def start(self, case) -> dict:
        """
        Returns
        -------
        Dictionary returned by Case.takestep
        """
        self._case = case

        # First step
        self._last = self._case.takestep(restart=1, ftol=1e-8, dtreal=1e-14, itermx=30)
        return self._last

    def next(self) -> dict:
        """
        Take another step

        Returns
        -------
        Dictionary returned by Case.takestep
        """
        # Start with previous values
        dtreal = self._last["dtreal"]
        ftol = self._last["ftol"]
        rlx = self._last["rlx"]

        # Diagnostics from last step that can be used
        # to inform the next step settings

        # Timestep multiplier estimate
        dtmul_yldot = rlx / self._last["yldot_maxabs"]
        # nonlinear steps taken
        dtmul_nni = exp(self.nni_sensitivity * (self.nni_offset - self._last["nni"]))

        if self._last["failed"]:
            # Candidates for timestep change that can be combined
            dtdiv_1 = (1.0 / dtmul_yldot) * self.yldot_factor_failed
            dtdiv_2 = (1.0 / dtmul_nni) * self.nni_factor_failed
            dtdiv_3 = self.dtdiv_failed

            # Weight between maximum and minimum
            factor = self.dtmaxmin_failed * max([dtdiv_1, dtdiv_2, dtdiv_3]) + (
                1 - self.dtmaxmin_failed
            ) * min([dtdiv_1, dtdiv_2, dtdiv_3])

            dtreal /= factor
        else:
            # Candidates for timestep change that can be combined
            dtmul_1 = dtmul_yldot * self.yldot_factor_succeeded
            dtmul_2 = dtmul_nni * self.nni_factor_succeeded
            dtmul_3 = self.dtmul_succeeded

            # Weight between maximum and minimum
            factor = self.dtmaxmin_succeeded * max([dtmul_1, dtmul_2, dtmul_3]) + (
                1 - self.dtmaxmin_succeeded
            ) * min([dtmul_1, dtmul_2, dtmul_3])

            dtreal *= factor

        self._last = self._case.takestep(restart=1, ftol=1e-8, dtreal=dtreal, itermx=30)
        return self._last

    @staticmethod
    def optimize(start_settings, processes: int, *case_args, **case_kwargs):
        """Use SciPy.optimize.differential_evolution to optimize the
        strategy settings. Needs to be given the settings to start
        with, and the arguments needed to construct a Case to be
        solved.

        Arguments
        ---------

        start_settings: List[float]
            The strategy settings to start
        processes : int
            The number of processes to use

        All following arguments and keywords are passed to the Case constructor.

        Example
        -------

        In a directory with an "input.yaml" file:

        import uetools
        uetools.UeCase.SolverStrategy.optimize(
            uetools.UeCase.SolverStrategy.good_settings[0],
            32,
            "input.yaml"
        )

        """
        from .Parallel import Pool
        import scipy.optimize as optimize

        pool = Pool(processes=processes, rate_limit=1.0)

        runner = RunSolver(*case_args, **case_kwargs)

        result = optimize.differential_evolution(
            runner,
            [(0, 1.0)] * 10,
            maxiter=1000,
            popsize=15,
            x0=start_settings,
            workers=lambda f, i: pool.map(f, i, timeout=300, default=1e6),
        )

        with open("result.pkl", "wb") as f:
            pickle.dump(result)


class RunSolver:
    def __init__(self, *args, **kwargs):
        self.case_args = args
        self.case_kwargs = kwargs

    def __call__(self, settings):
        from . import Case
        from uedge import bbb
        import uuid
        import pickle
        import time

        score = 1e6
        try:
            bbb.iprint = 0  # Turn off printing from UEDGE

            solver = SolverStrategy(settings)
            case = Case(*self.case_args, **self.case_kwargs, verbose=False)

            # Take first step
            result = solver.start(case)

            # Make list of results for each step
            steps = [result]

            iters = 0
            start_time = time.time()
            total_nfe = 0  # Total number of function evaluations
            while iters < 100:
                result = solver.next()
                total_nfe += result["nfe"]
                steps.append(result)
                iters += 1
                if result["success"] and result["dtreal"] > 1.0:
                    # Successfully took a timestep over 1s
                    score = time.time() - start_time
                    break

            total_time = time.time() - start_time

            # Pickle for later analysis
            filename = str(uuid.uuid4()) + ".pkl"
            with open(filename, "wb") as f:
                pickle.dump(
                    {
                        "settings": settings,
                        "steps": steps,
                        "time": total_time,
                        "nfe": total_nfe,
                        "score": score,
                    },
                    f,
                )
        except Exception as e:
            print(f"Exception {e}")

        return score
