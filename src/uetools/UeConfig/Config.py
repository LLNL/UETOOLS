# Reads/creates a config environment for UEDGE
from uetools.UeUtils.Lookup import Lookup


class Config(Lookup):
    def configcase(self, verbose=True, new=True, **kwargs):
        from os import path
        from yaml import safe_load
        from pathlib import Path

        # True if succeeds
        self.configured = False
        try:
            self.verbose
            verbose = self.verbose
        except:
            pass

        searchpath = path.expanduser("~")
        super().__init__()
        try:
            config = safe_load(Path("{}/.uetoolsrc".format(searchpath)).read_text())
            if verbose is True:
                print("UEDGE configuration file {}/.uetoolsrc read.".format(searchpath))
        except:
            print("No UETOOLS config file found: Configure file by calling Case.CreateConfig()")
            print("Alternatively, manually create the .uetoolsrc configuration YAML in your home directory.")
            config = None
        if config is not None:
            for key, variable in config.items():
                setattr(self, key, variable)
            # NOTE: what other information to write/store?

    def CreateConfig(self):
        from os import path
        from yaml import dump

        paths = {}
        searchpath = path.expanduser("~")
        yes = ["yes", "y"]
        no = ["no", "n"]
        print("Do you want to create a configuration file at "+\
                "{}? [y/n]".format(searchpath))
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
            print('Do you want to define a path to a YAML variable file?')
            print('(Used to deifne the variables to be saved by UETOOLS)')
            create = input()
            defpath='x'
            if (create.lower() in yes) or (len(create) == 0):
                while not path.exists(defpath):
                    print("Define path to variableyaml:")
                    defpath = input()
                    defpath = defpath.replace("~", searchpath)
                    if path.exists(defpath) is False:
                        print("Directory does not exist, please try again.")
                    else:
                        paths['variableyamlfile'] = defpath
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

    def configured(self):
        return self.configured
