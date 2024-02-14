"""
Routines for executing multiple UEDGE sessions in parallel

UeCase.ParallelCase is used in the same way as UeCase.Case, but
uses multiprocessing to launch a subprocess, and in the subprocess
create a UeCase.Case object.

Method calls on the ParallelCase object are forwarded through a
pipe to the subprocess.

Example
-------

In this process:
>>> case = uetools.UeCase.Case("186841/nc11/nc11_UeDB.hdf5", inplace=True)

In a new process:
>>> pcase = uetools.UeCase.ParallelCase("186841/nc11/nc11_UeDB.hdf5", inplace=True)

pcase and case should now produce the same results:
>>> case.get('te')
array([[3.75810267e-19, 3.76900342e-19, 4.47382024e-19, ...,

>>> pcase.get('te')
array([[3.75810267e-19, 3.76900342e-19, 4.47382024e-19, ...,

Docstrings are also copied, so
>>> help(case)
and
>>> help(pcase)
should contain the same docstring. Similarly for most methods.

"""

from .Case import Case

import multiprocessing as mp
import inspect
import time
import uuid
import pickle

from typing import Optional



def _case_handler(pipe, args, kwargs):
    """
    This function is launched on the remote process, and marshals
    data through the pipe
    """
    import sys
    from io import StringIO

    # Capture stdout into a string buffer.
    # - If instead we used a pipe to send output to the main process
    #   then we might block when the pipe is full
    # - When polled we will send the contents of this buffer
    # - This may not redirect Forthon code output
    sys.stdout = StringIO()

    case = Case(*args, **kwargs)
    while True:
        # Wait for a method call
        name, args, kwargs = pipe.recv()

        # Handle special cases
        if name == "_get_stdout":
            # Get the contents of the stdout buffer
            result = sys.stdout.getvalue()
            # Reset buffer
            sys.stdout = StringIO()
        elif name == "_exit":
            # Exit this process
            return
        else:
            # A method call for the Case object
            try:
                result = getattr(case, name)(*args, **kwargs)
            except Exception as e:
                result = e  # Can be pickled?

        # Send the result back
        pipe.send(result)


class ParallelCase:
    """
    A wrapper around Case that uses multiprocessing to start a subprocess
    and forwards method calls through a pipe.
    """

    def __init__(self, **kwargs):
        """
        Initializes a UeCase object in a separate process
        """

        # Get a list of methods that UeCase.Case supports
        # Extract the docstring for each method
        self._case_methods = {
            name: {"doc": func.__doc__}
            for name, func in inspect.getmembers(Case, inspect.isfunction)
        }

        # Special cases. These methods are added in the Case constructor.
        self._case_methods["get"] = {"doc": ""}
        self._case_methods["set"] = {"doc": ""}
        self._case_methods["getue"] = {"doc": ""}
        self._case_methods["setue"] = {"doc": ""}

        # Copy the class docstring from Case
        self.__doc__ = Case.__doc__

        # Create a pipe for communication with the subprocess
        parent_conn, child_conn = mp.Pipe()

        # Start a new process, forwarding the constructor arguments
        self.process = mp.Process(target=_case_handler, args=(child_conn, args, kwargs))
        self.process.start()
        self.pipe = parent_conn

    def getstdout(self):
        """
        Gets any text written to stdout by the subprocess since the last call.
        This is buffered so that processes don't block when printing.
        """
        self.pipe.send(("_get_stdout", [], {}))
        return self.pipe.recv()

    def __getattr__(self, name):
        if name not in self._case_methods:
            raise AttributeError("Attribute " + name + " not in ParallelCase")

        properties = self._case_methods[name]

        class CaseMethod:
            """Represents a method on a remote Case object that can be
            called through a pipe.

            The method call is blocking, waiting for the reply from
            the subprocess.
            """

            def __init__(self, pipe, name):
                self.pipe = pipe  # The bidirectional pipe for communication
                self.name = name  # The method to be called

            def __call__(self, *args, **kwargs):
                # Send the method and arguments
                self.pipe.send((self.name, args, kwargs))
                # Wait for the reply
                return self.pipe.recv()

        method = CaseMethod(self.pipe, name)
        method.__doc__ = properties["doc"]
        return method

    def __del__(self):
        """
        Clean up, close process
        """
        # First ask nicely
        self.pipe.send(("_exit", [], {}))
        time.sleep(0.1)
        # Then force process to close
        self.process.terminate()


def _run_func_and_send(pipe, func, args):
    """
    Internal function used by Pool on subprocesses

    A standalone function is used rather than a method, so that the Pool
    object (self) does not need to be pickled.
    """
    result = func(args)
    pipe.send(result)


class Pool:
    """
    A replacement for multiprocessing.Pool
    https://docs.python.org/3/library/multiprocessing.html#multiprocessing.pool.Pool

    Includes ability to limit the rate that new processes are spawned,
    and a timeout. Always starts every task in a new process, like
    setting maxtasksperchild=1 and chunksize=1 in multiprocessing.Pool.

    """

    def __init__(self, processes: int = 1, rate_limit: float = 0.0):
        """Initialize Pool object

        Arguments
        ---------

        processes: int
            The maximum number of processors to use. Must be >= 1

        rate_limit: float
            The shortest time (in seconds) to wait between launching
            processes. This is to avoid starting too many processes
            simultaneously that all require the same resources.

        """
        assert processes >= 1

        self.max_processes = processes
        self.rate_limit = rate_limit
        # create a process context
        self.ctx = mp.get_context("spawn")

    def map(self, func, iterable, timeout: Optional[float] = None, default=None):
        """Maps a function `func` over each item in the iterable.

        Arguments
        ---------

        func : callable
            A function that will be called
        iterable
            A list or other iterable of items to pass to `func`

        Keyword arguments
        -----------------

        timeout : float or None
            Time (in seconds) that processes will be allowed to run
            before being killed.

        default
            The value returned when a subprocess failed or timed out.

        Returns
        -------
        A list of the same length as `iterable`

        """
        active_pipes = []  # Connections to receive the results
        active_index = []  # Index in the input
        active_process = []  # Process handles
        active_started = []  # Start times

        last_started_time = 0.0

        results = [default for i in range(len(iterable))]

        def wait_for_results():
            """Waits for one or more processes to finish.
            If `timeout` is not None, then it kills any processes
            that have been running for longer than `timeout` seconds.

            Modifies
            --------

            active_pipes: list of Connections
                One for each running process
                https://docs.python.org/3/library/multiprocessing.html#multiprocessing.connection.Connection
            active_index: list of ints
                Index into the `result` list for each running process
            active_process: list of Processes
                Handle for each Process that is still running
            active_started: list of floats
                The time in seconds when each running process was started

            The lengths of these lists should always be the same, and
            the lists syncronised so e.g. active_pipes[i] is the
            connection to active_process[i] that started at time
            active_started[i], whose result will be put into
            result[active_index[i]]

            """
            # Wait for a process to finish
            print(f"Pool: Waiting on {len(active_pipes)} processes")
            ready_pipes = mp.connection.wait(active_pipes, timeout=timeout)
            current_time = time.time()
            print(f"Pool: {len(ready_pipes)} processes ready")

            for pipe in ready_pipes:
                try:
                    result = pipe.recv()
                except Exception as e:
                    print(f"Pool: Exception in pipe.recv: {e}")
                    result = default

                task_index = active_pipes.index(pipe)
                pipe = active_pipes.pop(task_index)
                try:
                    pipe.close()
                except Exception as e:
                    print(f"Pool: Exception in pipe.close: {e}")

                result_index = active_index.pop(task_index)
                results[result_index] = result
                process = active_process.pop(task_index)
                try:
                    process.kill()
                except Exception as e:
                    print(f"Pool: Exception in process.kill: {e}")
                start_time = active_started.pop(task_index)

            if timeout is not None:
                # Check for processes that have been running too long
                inds = [
                    ind
                    for ind, time in enumerate(active_started)
                    if time < (current_time - timeout)
                ]
                print(f"Pool: {len(inds)} processes timed out")
                # Remove in reverse order so that indices are valid
                for ind in reversed(inds):
                    pipe = active_pipes.pop(ind)
                    pipe.close()

                    result_index = active_index.pop(ind)
                    # Result is already default

                    process = active_process.pop(ind)
                    process.kill()

                    start_time = active_started.pop(ind)

        for index, item in enumerate(iterable):
            if len(active_process) >= self.max_processes:
                wait_for_results()

            # At least one process has finished
            # => Can now start a new process
            current_time = time.time()
            if current_time - last_started_time < self.rate_limit:
                # Wait until the minimum time has elapsed
                time.sleep(self.rate_limit - (current_time - last_started_time))

            # Create a unidirectional pipe (duplex=False)
            receiver, sender = mp.Pipe(False)

            # Start new process
            proc = self.ctx.Process(
                target=_run_func_and_send, args=(sender, func, item)
            )
            proc.start()
            last_started_time = time.time()

            active_pipes.append(receiver)
            active_index.append(index)
            active_process.append(proc)
            active_started.append(last_started_time)

        # All tasks now allocated. Wait for them to finish
        while active_process:
            print("Pool: Draining")
            wait_for_results()

        return results

