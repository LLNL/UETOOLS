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
