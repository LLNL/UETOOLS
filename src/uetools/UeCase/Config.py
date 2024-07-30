# Reads/creates a config environment for UEDGE
from uetools.UeUtils.Lookup import Lookup


class Config:
    """Class for creating and reading .uedgerc config file

    Methods
    -------
    case(verbose=True, **kwargs)
        Looks for ~/.uedgerc. Reads it if found, prints instructions
        for creating it if not
    create()
        Creates ~/.uedgerc interactively from Python prompt

    """

    def __init__(self):
        """Creates Config object and links functions"""
        self.configured = False
        self.search = Lookup()
        self.configs = {}

    def case(self, verbose=True, **kwargs):
        """Looks for and reads ~/.uedgerc

        Prints instructions to creating it if file not found

        Keyword arguments
        -----------------
        verbose : bool (default = False)
            Switch whether to write progress to stdout (True) or not

        Returns
        -------
        None
        """
        from os import path
        from yaml import safe_load
        from pathlib import Path

        # True if succeeds
        self.configured = False

        searchpath = path.expanduser("~")
        try:
            config = safe_load(Path("{}/.uetoolsrc".format(searchpath)).read_text())
            if verbose is True:
                print("UEDGE configuration file {}/.uetoolsrc read.".format(searchpath))
        except:
            if verbose is True:
                print(
                    "No UETOOLS config file found: Configure file by calling Case().config.create() or uetools.config.create()"
                )
                print(
                    "Alternatively, manually create the .uetoolsrc configuration YAML in your home directory."
                )
            config = None
        if config is not None:
            for key, variable in config.items():
                self.configs[key] = variable
            # NOTE: what other information to write/store?

    def create(self):
        """Creates ~/.uetoolsrc based on queried input"""
        from os import path
        from yaml import dump

        paths = {}
        searchpath = path.expanduser("~")
        yes = ["yes", "y"]
        no = ["no", "n"]
        print(
            "Do you want to create a configuration file at "
            + "{}? [y/n]".format(searchpath)
        )
        try:
            # Note: This can fail in subprocess
            create = input()
        except:
            create = "n"
        if (create.lower() in yes) or (len(create) == 0):
            for dirpath in ["aphdir", "apidir"]:
                defpath = "x"
                while not path.exists(defpath):
                    print('Define path to "{}":'.format(dirpath))
                    defpath = input()
                    defpath = defpath.replace("~", searchpath)
                    if path.exists(defpath) is False:
                        print("Directory does not exist, please try again.")
                    else:
                        paths[dirpath] = path.abspath(defpath)
                        print("    Path defined successfully!")
            print("Do you want to define a path to a YAML variable file?")
            print("(Used to deifne the variables to be saved by UETOOLS)")
            create = input()
            defpath = "x"
            if (create.lower() in yes) or (len(create) == 0):
                while not path.exists(defpath):
                    print("Define path to variableyaml:")
                    defpath = input()
                    defpath = defpath.replace("~", searchpath)
                    if path.exists(defpath) is False:
                        print("Directory does not exist, please try again.")
                    else:
                        paths["variableyamlfile"] = defpath
                        print("    Path defined successfully!")
            else:
                print("Using standard YAML variable file.")
        else:
            print("Please create .uetoolsrc manually in your home directory")
            print("Aborting!")
            self.configured = False
            return False
        with open("{}/.uetoolsrc".format(searchpath), "w") as file:
            dump(paths, file)
        print("UEDGE config file {}.uetoolsrc successfully created!".format(searchpath))
