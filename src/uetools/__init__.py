# Try to import UEDGE packages into the namespace
try:
    from uedge import bbb, com, grd, flx, aph, api
# In case no UEDGE is installed, only standalone mode is allowed
except:
    print('No UEDGE install found. Importing UEDGE Toolbox in standalone mode')
from .UeCase import Case
from .UeDatabase import Database
from . import UeCase
from . import UeConfig
#from . import Dashboard
from .UeDashboard import uedashboard
from .UeDashboard import StandaloneDashboard
from .UeDashboard import StandaloneDatabaseDashboard
from os import path

with open(path.join(__path__[0],"VERSION")) as f:
    __version__ = f.read().replace('\n', '').strip()
