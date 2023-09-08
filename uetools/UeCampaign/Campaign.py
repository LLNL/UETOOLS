from typing import Tuple, Callable, TypeVar
from pathlib import Path
import random
import uuid
import time
import numpy as np

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
        - "filename": The file, same as the key.
        - "previous": Name of the previous state file, or None if initial state.
        - "values": A dictionary of values, with the same keys as the target.
        - "converged": Boolean, whether converged or not.
        - "distance": A measure of distance between values and target

    in_progress: dict
        One entry per currently running calculation. Key is the resulting file name.
        Each entry contains:
         - 

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

        # No simulations currently running or queued
        self.in_progress = {}

    def closest_state(self) -> dict:
        """
        Return information about the state that is closest to the target
        """
        closest = None
        for state in self.states.values():
            if (closest is None) or (state["distance"] < closest["distance"]):
                closest = state
        return closest

    def closest_distance(self) -> float:
        """
        Return a measure of the closest distance to the target 
        """
        return self.closest_state()["distance"]

    def fraction_progress(self) -> float:
        """
        Returns a measure of progress made. 0 at the start of the campaign,
        1 when the target values are reached
        """
        return 1.0 - self.closest_distance() / distance(self.start_values, self.target)

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
        self.in_progress = {}

    def next_action(self) -> Tuple[Path, dict]:
        """
        Using the dictionary of datasets and in_progress simulations,
        work out what simulation to run next.

        Returns the path to the starting dataset, and a dictionary
        of values to try and reach.
        """

        # Choose a starting dataset, preferring those that are close
        # to the target

        def choose(states: list[T], weight_fn: Callable[[T], float]) -> T:
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

        state = choose(
            list(self.states.values()), lambda state: 1.0 / (state["distance"] + 1e-10)
        )  # Avoid divide-by-zero

        if not state["converged"]:
            # Failed to converge. Check if the previous case converged
            previous = self.states[state["previous"]]
            if previous["converged"]:
                # Previous case did converge.
                # Try again with the same settings
                return (state, state["values"])

            # Previous also failed. Go back and try a different way
            # Question: Should we delete the failed cases?
            #           Or try going in the same direction with a smaller step?
            state = self.states[previous["previous"]]

        # Have a converged starting case.
        # Decide whether to try the target values. The closer we are the higher the chance.
        d = distance(state["values"], self.target) / distance(
            self.start_values, self.target
        )
        if random.uniform(0, 1) > d:
            # Try jumping to the target
            return (state, self.target)

        # Choose a direction to go in towards the target.
        values = {}
        for var, target in self.target.items():
            weight = random.uniform(0, 1)
            values[var] = weight * state["values"][var] + (1.0 - weight) * target
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
        input_state, changes = self.next_action()

        # Choose a file name for the result
        output_file = self.path / (str(uuid.uuid4()) + ".hdf5")

        # Submit job to the pool
        async_result = pool.apply_async(
            run_case, (input_state["filename"], changes, output_file)
        )

        # Store the handle and some information about the run
        new_run = {
            "filename": filename,
            "previous": case["filename"],
            "values": values,
            "async_result": async_result,
            "submitted_time": time.time(),
        }
        self.in_progress[filename] = new_run
        return new_run

    def check_async(self) -> list[dict]:
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
                    finished.append({"filename": run["filename"], "converged": False})
        # Remove finished runs from in_progress
        for run in finished:
            del self.in_progress[run["filename"]]
        return finished
