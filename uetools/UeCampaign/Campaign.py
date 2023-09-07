from typing import Tuple, Callable
from pathlib import Path
import random
import uuid
import time


def distance(values: dict, targets: dict) -> float:
    """
    Calculate a metric of distance between given values and
    target values

    Both input dictionaries should have the same keys
    """
    reldiff = 0.0
    for var, target in targets.items():
        value = values[var]
        reldiff += np.abs(target - value) / (np.abs(target) + np.abs(value) + 1e-10)
    return reldiff


class Campaign:
    """

    Members
    -------

    path: Path
        The directory where datasets are stored (absolute path)

    datasets: dict
       - key: File name
       - dictionary of values
       - Whether converged or not
       - The last convergence report, containing e.g. norms
       - distance measure from values to target
    
    """

    def __init__(self, case_args, case, target: dict, path: Path):
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
        self.case_args = case_args
        self.path = path
        self.target = target

        # Create the initial state
        filename = path / "init.hdf5"
        case.dump(filename)

        values = {var: case.getue_memory(var) for var in target}
        self.datasets = {
            filename: {
                "filename": filename,
                "previous": None,
                "values": values,
                "converged": True,
                "distance": distance(values, target),
            }
        }

        # Simulations currently running or queued
        self.in_progress = {}

    def __getstate__(self):
        """
        Pickle this class, saving only the datasets
        """
        return {
            "case_args": self.case_args,
            "path": self.path,
            "target": self.target,
            "datasets": self.datasets,
        }

    def __setstate__(self, state):
        """
        Restore state from e.g. pickle file
        """
        self.case_args = state["case_args"]
        self.path = state["path"]
        self.target = state["target"]
        self.datasets = state["datasets"]
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

        def choose(cases, weight_fn):
            """
            Choose a random case using a function to calculate the weight
            """
            cases_list = list(cases)
            weights = [weight_fn(case) for case in cases_list]
            sum_weight = sum(weights)
            r = random.uniform(0.0, sum_weight)
            for case, weight in zip(cases_list, weights):
                r -= weight
                if r <= 0.0:
                    return case

        case = choose(self.datasets, lambda case: 1.0 / case["distance"])

        if not case["converged"]:
            # Failed to converge. Check if the previous case converged
            previous = self.datasets[case["previous"]]
            if previous["converged"]:
                # Previous case did converge.
                # Try again with the same settings
                return (case, case["values"])

            # Previous also failed. Go back and try a different way
            # Question: Should we delete the failed cases?
            #           Or try going in the same direction with a smaller step?
            case = self.datasets[previous["previous"]]

        # Have a converged starting case. Choose a direction
        # to go in.

        values = {}
        for var, target in self.target.items():
            weight = random.uniform(0, 1)
            values[var] = weight * case["values"][var] + (1.0 - weight) * target

        return (case, values)

    def start_new_case(self, pool, case, values: dict) -> dict:
        """
        Start a new case using the given Pool. Inserts a new case
        into self.in_progress.

        pool: Multiprocessing Pool
           .apply_async()

        Returns the case dictionary
        """
        # Choose a file name for the new result
        filename = self.path / (str(uuid.uuid4()) + ".hdf5")

        async_result = pool.apply_async()

        new_case = {
            "filename": filename,
            "previous": case["filename"],
            "values": values,
            "async_result": async_result,
            "submitted_ts": time.time(),
        }
        self.in_progress[filename] = new_case
        return new_case


def console(campaign: Campaign):
    """
    Display a console view of the Campaign
    """
