# Reads/creates a config environment for UEDGE
from uetools.UeUtils.Lookup import Lookup


class Config(Lookup):
    def configcase(self, verbose=True, **kwargs):
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
            config = safe_load(Path("{}/.uedgerc".format(searchpath)).read_text())
            if verbose is True:
                print("UEDGE configuration file {}/.uedgerc found.".format(searchpath))
        except:
            if self.createuedgerc() is False:
                return False
            else:
                config = safe_load(Path("{}/.uedgerc".format(searchpath)).read_text())

        for dirpath in ["aphdir", "apidir"]:
            packobj = self.getpackobj(dirpath)
            try:
                strlen = len(packobj.getpyobject(dirpath)[0])
                packobj.getpyobject(dirpath)[0] = config[dirpath].ljust(strlen)
            except:
                print(
                    'Required path "{}" not found in .uedgerc. Aborting!'.format(
                        dirpath
                    )
                )
                return False
        # NOTE: what other information to write/store?
        return True

    def createuedgerc(self):
        from os import path
        from yaml import dump

        paths = {}
        searchpath = path.expanduser("~")
        yes = ["yes", "y"]
        no = ["no", "n"]
        print("UEDGE config file not found!")
        print("Create it here? [y/n]")
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
        else:
            print("Please create .uedgerc manually in your home directory")
            print("Aborting!")
            self.configured = False
            return False
        with open("{}/.uedgerc".format(searchpath), "w") as file:
            dump(paths, file)
        print("UEDGE config file .uedgerc successfully created!")

    def configured(self):
        return self.configured
