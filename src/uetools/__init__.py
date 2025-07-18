# Try to import UEDGE packages into the namespace
try:
    from uedge import bbb, com, grd, flx, aph, api
# In case no UEDGE is installed, only standalone mode is allowed
except:
    print('No UEDGE install found. Importing UEDGE Toolbox in standalone mode')
try:
    from uedge import ppp
except:
    # Using version of UEDGE not contining the ppp package
    pass

from .UeCase import Case, Config
from .UeDatabase import Database
from . import UeCase
from . import UeCase
from .UeCase import Parallel
from .UeCampaign import Campaign
try:
    from .UeDashboard import uedashboard
except:
    pass
try:
    from .UeDashboard import StandaloneDashboard
except:
    pass
try:
    from .UeDashboard import StandaloneDatabaseDashboard
except:
    pass
from os import path

with open(path.join(__path__[0],"VERSION")) as f:
    __version__ = f.read().replace('\n', '').strip()

config = Config()
config.case()
