from typing import Tuple, Callable, TypeVar, List, Dict, Optional
from pathlib import Path
import random
import uuid
import time
import datetime
import numpy as np
import pickle

T = TypeVar("T")


def distance(values: dict, targets: dict) -> float:
    """
    Calculate a metric of distance between given values and
    target values

    Both input dictionaries should have the same keys
    """
    reldiff = 0.0
    for var, target in targets.items():
        value = values[var]
        if np.isscalar(target) and not np.isscalar(value):
            value = value[0]
        reldiff += np.sum(np.abs(target - value)) / (
            np.sum(np.abs(target)) + np.sum(np.abs(value)) + 1e-10
        )
    return reldiff


def converge_case(case) -> Tuple[bool, dict]:
    """
    Returns True as first argument if converged successfully.
    Second return argument describes the convergence state
    """
    import uetools

    s = uetools.UeCase.SolverStrategy(uetools.UeCase.SolverStrategy.good_settings[0])

    s.start(case)
    for i in range(100):
        result = s.next()
        if result["success"] and result["dtreal"] > 1.0:
            return (True, result)
    return (False, result)


def run_case(input_file, changes, output_file):
    """
    - Read state from an input file
    - Change settings
    - Run UEDGE and attempt to converge
    - Save state to an output file

    """
    import uetools

    # Load the state
    # Note: Case can't handle Path or PosixPath objects because HDF5 tries to hash them.
    case = uetools.Case(str(input_file))

    # Take a small step to ensure that we're syncronised
    case.takestep(restart=1, ftol=1e-8, dtreal=1e-7, itermx=30)

    # Modify the settings
    for var, value in changes.items():
        if var in ["pcoree", "pcorei"]:
            # Setting a core power
            case.setue_memory("iflcore", 1)
        elif var in ["tcoree", "tcorei"]:
            # Setting a core temperature
            case.setue_memory("iflcore", 0)
        case.setue_memory(var, value)

    case.takestep(restart=1, ftol=1e-8, dtreal=1e-7, itermx=30)
    converged, solver_info = converge_case(case)

    # Save state to the  output file
    case.save(output_file)

    return {
        "filename": output_file,
        "previous": input_file,
        "values": changes,
        "converged": converged,
        "solver_info": solver_info,
    }


class Campaign:
    """
    Tries to find a converged UEDGE state with a specified target
    set of parameters.

    Uses a random walk process, optionally using a pool of processes
    to explore paths to the target parameters in parallel.

    Members
    -------

    path: Path
        The directory where datasets are stored (absolute path)

    target: dict
        Keys are the parameter to vary, values are the target values.

    start_values: dict
        Same keys as target, containing values from initial state.
        Can be used to calculate a measure of progress.

    states: dict
        Known state files. Key is File name.
        Each entry contains:
        - "filename"   The file, same as the key.
        - "previous"   Name of the previous state file, or None if initial state.
        - "values"     A dictionary of values, with the same keys as the target.
        - "converged"  Boolean, whether converged or not.
        - "distance"   A measure of distance between values and target

    in_progress: dict
        One entry per currently running calculation. Key is the resulting file name.
        Each entry contains:
         - "filename"        The file that will be written when finished
         - "previous"        Existing file containing the starting state
         - "values"          Parameter values being attempted
         - "async_result"    An AsyncResult object representing future result
         - "submitted_time"  Time when the job was submitted

    Example
    -------

    Starting from a Case object, case, begin a campaign
    to get to core electron power of 1MW and core density 6e19m^-3.
    Store state files in the current directory ".":

        import uetools
        campaign = uetools.Campaign(case, {'pcoree': 1e6, 'ncore': 6e19}, ".")

    Running in serial
    ~~~~~~~~~~~~~~~~~

    To run the Campaign in serial until it reaches the target values:

        campaign.run()

    To perform one step at a time:

        campaign.run_next()

    Will modify the settings and try to converge UEDGE.
    It will return a dictionary describing the run, and add the new
    state to the campaign.states dictionary.

    Progress towards the target values can be monitored using
    campaign.fraction_progress(), which returns 0 at the start and 1
    when target is reached.

    Running in parallel
    ~~~~~~~~~~~~~~~~~~~

    Campaigns perform a kind of random walk search to find a way
    to get to the target values. This search can be run in parallel
    with a pool of workers. The primary process coordinates workers,
    and decides what to try next based on the results of all workers.

    To run asyncronously with 2 processors:

        campaign.run_async(processes = 2)

    While this is running it will display a status line with
    the elapsed time, number of states calculated, and progress e.g:

        00:03:51 | 003 states |###########>         |55.0%

    When the target is reached the file containing the result is printed:

        Final state in file: c238b767-ca33-46d6-8c6c-746dbfe5de12.hdf5

    The state can also be retrieved by calling campaign.closest_state()

    To run asyncronously while doing other things, you can set up a
    multiprocessing pool:

        pool = uetools.Parallel.QuietPool(processes=2)

    then launch one job at a time:

        campaign.run_next_async(pool)

    which returns information about the submitted job.
    It's important to check for finished jobs to update the states:

        campaign.check_async()

    That will return a list of finished jobs.

    Output logging
    ~~~~~~~~~~~~~~

    The default run_async method will write to output every 1/2 second
    as it updates the status display. This quickly adds up if the output
    is being logged. To only print when something interesting happens,
    set output='log':

        campaign.run_async(processes=2, output='log')

    This will print something like:

        2023-09-07 23:31:28.074046: Started cb26 (87.6%) from init (0.0%)
        2023-09-07 23:31:28.074760: Started fae0 (21.2%) from init (0.0%)
        2023-09-07 23:35:48.201753: Finished cb26 (87.6%)
        2023-09-07 23:35:48.202120: Started 869a (69.8%) from init (0.0%)
        2023-09-07 23:35:54.705231: Finished fae0 (21.2%)
        2023-09-07 23:35:54.705568: Started 48ec (100.0%) from cb26 (87.6%)
        2023-09-07 23:40:07.828991: Finished 869a (69.8%)
        2023-09-07 23:40:07.829406: Started 2482 (100.0%) from cb26 (87.6%)
        2023-09-07 23:40:22.336841: Finished 48ec (100.0%)
        Final state in file: 48ec63ab-f48b-4435-a7b4-b68ba1024fcd.hdf5

    Where 4 letters of the start of the file name are printed, along with
    the percentage progress towards the target.

    Saving and restoring
    ~~~~~~~~~~~~~~~~~~~~

    Campaign objects can be pickled, or serialised
    and stored in another format. That file will contain a record of all
    the datasets, and can be unpickled to continue a campaign.

        # To save:
        import pickle
        with open("campaign.pickle", "wb") as f:
            pickle.dump(campaign, f)

        # To restore:
        import pickle
        with open("campaign.pickle", "rb") as f:
            campaign = pickle.load(f)

    """

    def __init__(self, case, target: dict, path: Path):
        """
        Start a new UEDGE campaign, from a given starting case
        and a dictionary of target values.

        case: Case
            A UEDGE simulation case
        target: dict
            Dictionary of UEDGE settings to aim for
        path: Path
            A directory to be used to store UEDGE states

        """

        path = Path(path).absolute()
        assert path.is_dir()
        self.path = path
        self.target = target

        # Store the initial state
        filename = path / "init.hdf5"
        case.save(filename)

        # Get the value of the parameters being varied
        self.start_values = {var: case.getue_memory(var) for var in target}
        self.states = {
            filename: {
                "filename": filename,
                "previous": None,  # No previous state => initial
                "values": self.start_values,
                "converged": True,
                "distance": distance(self.start_values, target),
            }
        }
        self.total_distance = distance(self.start_values, self.target)

        # No simulations currently running or queued
        self.in_progress = {}

    def closest_state(self) -> Optional[Dict]:
        """ Return information about the converged state that is
        closest to the target

        May return None if no states have converged
        """
        closest = None
        for state in self.states.values():
            if not state["converged"]:
                continue  # Skip non-converged cases
            if (closest is None) or (state["distance"] < closest["distance"]):
                closest = state
        return closest

    def closest_distance(self) -> float:
        """
        Return a measure of the closest distance to the target
        """
        state = self.closest_state()
        if state is None:
            return float('inf')
        return state["distance"]

    def fraction_progress(self) -> float:
        """
        Returns a measure of progress made. 0 at the start of the campaign,
        1 when the target values are reached
        """
        return 1.0 - self.closest_distance() / self.total_distance

    def __str__(self):
        """
        Print a human-readable summary
        """
        # First a table of targets, start and current best values

        closest = self.closest_state()
        closest_values = closest["values"]
        s = ""
        for var, target in self.target.items():
            s += f"{var}\t{self.start_values[var]}\t{closest_values[var]}\t{target}\n"

        return s

    def __getstate__(self):
        """
        Pickle this class, saving only the datasets
        """
        return {
            "path": self.path,
            "target": self.target,
            "states": self.states,
            "start_values": self.start_values,
        }

    def __setstate__(self, state):
        """
        Restore state from e.g. pickle file
        """
        self.path = state["path"]
        self.target = state["target"]
        self.states = state["states"]
        self.start_values = state["start_values"]
        self.total_distance = distance(self.start_values, self.target)
        self.in_progress = {}

    def next_action(self) -> Tuple[Path, dict]:
        """
        Using the dictionary of datasets and in_progress simulations,
        work out what simulation to run next.

        Returns the path to the starting dataset, and a dictionary
        of values to try and reach.

        Parameters
        ----------

        prefer_near : float, > 0
            Higher values increase the chance of states near the target being chosen
            i.e. more effort is put into advancing existing paths, rather than
            starting new ones

        prefer_no_jump : float, > 0
            Higher values reduce the chance of jumping to the target parameters

        prefer_small_step : float, > 0
            Higher values reduce the average step size
        """

        prefer_near = 2
        prefer_no_jump = 2
        prefer_small_step = 3

        # Choose a starting dataset, preferring those that are close
        # to the target

        def choose(states: List[T], weight_fn: Callable[[T], float]) -> T:
            """
            Choose a random case using a function to calculate the weight
            """
            states_list = list(states)
            weights = [weight_fn(state) for state in states_list]
            sum_weight = sum(weights)
            r = random.uniform(0.0, sum_weight)
            for state, weight in zip(states_list, weights):
                r -= weight
                if r <= 0.0:
                    return state
            return states_list[0]

        state = choose(
            list(self.states.values()),
            lambda state: 1.0 / (state["distance"] + 1e-10) ** prefer_near,
        )  # Avoid divide-by-zero

        if not state["converged"]:
            # Failed to converge. Check if the previous case converged
            previous = self.states[state["previous"]]
            if previous["converged"]:
                # Previous case did converge.
                # Try again with the same settings
                # and set distance to a large value so it's not chosen again
                state["distance"] = 1e20
                return (state, state["values"])

            # Previous also failed. Go back and try a different way
            # Question: Should we delete the failed cases?
            #           Or try going in the same direction with a smaller step?
            # Set distances to large values, so these states aren't used again
            state["distance"] = 1e20
            previous["distance"] = 1e20
            state = self.states[previous["previous"]]

        # Have a converged starting case.
        # Decide whether to try the target values. The closer we are the higher the chance.
        d = distance(state["values"], self.target) / self.total_distance
        if random.uniform(0, 1) ** prefer_no_jump > d:
            # Try jumping to the target
            return (state, self.target)

        # Choose a direction to go in towards the target.
        values = {}
        for var, target in self.target.items():
            # Note: Higher powers bias towards smaller steps
            weight = random.uniform(0, 1) ** prefer_small_step
            values[var] = (1.0 - weight) * state["values"][var] + weight * target
        return (state, values)

    def run_next(self) -> dict:
        """
        Run the next case syncronously,
        returning information about the resulting state.
        """

        # Get the next action to perform
        input_state, changes = self.next_action()
        input_file = input_state["filename"]

        # Choose a file name for the result
        output_file = self.path / (str(uuid.uuid4()) + ".hdf5")

        # Run the case syncronously
        result = run_case(input_file, changes, output_file)

        # Sanity check round-trip
        assert result["previous"] == input_file
        assert result["filename"] == output_file

        # Add the "distance" metric and timestamp
        result["distance"] = distance(result["values"], self.target)
        result["finished_time"] = time.time()

        # Add to dictionary of states
        self.states[output_file] = result
        return result

    def run(self):
        """
        Run this campaign in serial until it finishes
        """
        while True:
            self.run_next()
            progress = self.fraction_progress()
            print(f"Campaign progress made: {progress * 100}%")
            if progress > 0.999:
                return

    def run_next_async(self, pool) -> dict:
        """
        Asyncronously start a new case using the given state, values and Pool.

        pool: Multiprocessing Pool
           Must have an apply_async(f, args) method.

        Returns a dictionary describing the simulation
        """

        # Get the next action to perform
        input_state, values = self.next_action()
        input_file = input_state["filename"]

        # Choose a file name for the result
        output_file = self.path / (str(uuid.uuid4()) + ".hdf5")

        # Submit job to the pool
        async_result = pool.apply_async(run_case, (input_file, values, output_file))

        # Store the handle and some information about the run
        new_run = {
            "filename": output_file,
            "previous": input_file,
            "values": values,
            "async_result": async_result,
            "submitted_time": time.time(),
        }
        self.in_progress[output_file] = new_run
        return new_run

    def check_async(self) -> List[Dict]:
        """
        Check if any asyncronous jobs have finished. If so,
        move them from `in_progress` to the `states` dictionary,
        and return information about the finished runs.
        """
        finished = []
        for run in self.in_progress.values():
            if run["async_result"].ready():
                try:
                    result = run["async_result"].get()
                    result["distance"] = distance(result["values"], self.target)
                    result["submitted_time"] = run["submitted_time"]
                    result["finished_time"] = time.time()
                    finished.append(result)

                    # Insert into dictionary of states
                    self.states[result["filename"]] = result
                except:
                    # This run raised an exception
                    finished.append(
                        {
                            "filename": run["filename"],
                            "submitted_time": run["submitted_time"],
                            "finished_time": time.time(),
                            "converged": False,
                            "distance": 1e20,
                        }
                    )
        # Remove finished runs from in_progress
        for run in finished:
            del self.in_progress[run["filename"]]
        return finished

    def run_async(
        self,
        processes: int = 2,
        delay: float = 0.5,
        output: str = "term",
        bar_length: int = 20,
        save_file: Optional[str] = None,
    ):
        """
        Run UEDGE simulations asyncronously until the target is reached

        processes: int
            Number of processors to use

        delay: float
            Time between updates, in seconds

        output: str
            One of 'term', 'log'

        bar_length: int
            For 'term' output, the length of the progress bar in chars

        save_file: str
            Pickle the campaign as new results come in.
            This helps if the campaign has to be stopped and restarted.

        """
        # Create a pool that suppresses output from subprocesses
        # Otherwise with many tasks the output is very cluttered
        import uetools

        pool = uetools.Parallel.QuietPool(processes)

        assert output in ["term", "log"]

        def log(msg: str):
            print(f"{datetime.datetime.now()}: " + msg)

        start_time = time.time()
        while True:
            while len(self.in_progress) < processes:
                # Start another job
                started = self.run_next_async(pool)

                if output == "log":
                    from_progress = 1.0 - (
                        distance(
                            self.states[started["previous"]]["values"], self.target
                        )
                        / self.total_distance
                    )
                    to_progress = 1.0 - (
                        distance(started["values"], self.target) / self.total_distance
                    )
                    log(
                        "Started "
                        + started["filename"].name[:4]
                        + f" ({to_progress * 100.:.1f}%) from "
                        + started["previous"].name[:4]
                        + f" ({from_progress * 100:.1f}%)"
                    )

            time.sleep(delay)

            # Check for finished jobs
            finished = self.check_async()
            progress = self.fraction_progress()

            if output == "term":
                # Update a simple one-line status with progress bar
                nruns = len(self.states) - 1  # Number of completed runs
                elapsed_time = datetime.timedelta(seconds=time.time() - start_time)

                nbars = int(bar_length * progress)
                bar = "|" + "#" * nbars + ">" + " " * (bar_length - nbars) + "|"

                print(
                    f"\r{elapsed_time} | {nruns:03} states "
                    + bar
                    + f"{progress*100:.1f}%",
                    end="",
                )

            elif output == "log":
                # Print logs of finished jobs
                for run in finished:
                    to_progress = 1.0 - (run["distance"] / self.total_distance)
                    log(
                        "Finished "
                        + run["filename"].name[:4]
                        + f" ({to_progress * 100.:.1f}%)"
                    )

            if len(finished) > 0 and save_file is not None:
                # Pickle the campaign
                with open(save_file, "wb") as f:
                    pickle.dump(self, f)

            if progress > 0.999:
                print("\nFinal state in file: " + str(self.closest_state()["filename"]))
                return
