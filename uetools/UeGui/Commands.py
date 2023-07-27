import tkinter as tk
import tkinter.ttk as ttk
from functools import partial
import uetools
import os
from os.path import exists
import sys
import matplotlib.pyplot as plt

try:
    from ruamel import yaml
except:
    import ruamel_yaml as yaml
import numpy as np
from uedge import *
from uedge.uedgeplots import *
from . import ScrollableNotebook as s
import os

root = None


class ToolTip(object):
    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0

    def showtip(self, text):
        "Display text in tooltip window"
        self.text = text
        if self.tipwindow or not self.text:
            return
        x, y, cx, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 57
        y = y + cy + self.widget.winfo_rooty() + 27
        self.tipwindow = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        tw.wm_geometry("+%d+%d" % (x, y))
        label = ttk.Label(
            tw, text=self.text, justify=tk.LEFT, relief=tk.SOLID, borderwidth=1
        )
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()


def CreateToolTip(widget, text):
    toolTip = ToolTip(widget)

    def enter(event):
        toolTip.showtip(text)

    def leave(event):
        toolTip.hidetip()

    widget.bind("<Enter>", enter)
    widget.bind("<Leave>", leave)


def cmdCallBack(cmd):
    exec(cmd)
    plt.pause(2)


class CmdGui(object):
    def __init__(self, yamlfilename):
        """
        This routine builds a graphics user interface as
        specified in the yamlfile.

        gui = uetools.UeGui.CmdGui(yamlfilename,nonblock=True)
        yamlfilename - either a yaml filename from the user space or
                       one that is installed with the module
        nonblock - True or False. If False then the Tkinter mainloop
                   is called and this class blocks until all windows
                   are closed. Setting to True results can vary but
                   the Python command interface should remain active.

        """
        global root
        root = tk.Tk()
        self.root = root
        self.root.title("Uedge Case Input")
        yamlfile = "{}/{}/{}".format(uetools.__path__[0], "yamls", yamlfilename)
        print(yamlfile)
        if exists(yamlfilename):
            self.filename = yamlfilename
        elif exists(yamlfile):
            self.filename = yamlfile
        else:
            print("File not found: ", yamlfilename)
            return
        with open(self.filename, "r") as f:
            self.config = yaml.load(f, Loader=yaml.Loader)

        notebook = s.ScrollableNotebook(self.root, wheelscroll=True, tabmenu=True)
        for tab in self.config.keys():
            frame = tk.Frame(notebook)
            notebook.add(frame, text=tab)
            notebook.pack(fill="both", expand=True)
            for i, en in enumerate(self.config[tab].keys()):
                state = "normal"
                cmd = self.config[tab][en]["cmd"]
                b1 = tk.Button(frame, text=en, command=partial(cmdCallBack, cmd))
                b1.grid(column=0, row=i, sticky="ew")
                if "help" in self.config[tab][en]:
                    l1 = tk.Label(frame, text="?", width=1)
                    l1.grid(column=1, row=i, sticky="e")
                    CreateToolTip(l1, text=self.config[tab][en]["help"])
        if "idlelib" not in sys.modules:
            print("")
            print("Blocking in mainloop()")
            print("Close all windows to return to commandline")
            print("Run IDLE/IDLE3 to avoid blocking")
            print("")
            self.root.mainloop()
