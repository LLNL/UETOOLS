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

    def __init__(self, *args, **kwargs):
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
