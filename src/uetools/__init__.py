"""
uetools package init.

Banner behavior:
- If UEDGE_BANNER=0: print nothing.
- If banner enabled:
  - Print ASCII logo unless UEDGE_LOGO=0.
  - Print versions unless UEDGE_VERSION=0.

Defaults (if not set by user/environment):
- UEDGE_LOGO=0
- UEDGE_VERSION=0

Additionally:
- Banner output is printed at most once per Python session.
"""

from __future__ import annotations
from os import path, environ
import os
with open(path.join(__path__[0],"VERSION")) as f:
    __version__ = f.read().replace('\n', '').strip()




# ---------- Printed-once guard (per Python session) ----------
# This flag lives in the module namespace, so repeated imports won't re-print.
# It will reset only when the interpreter restarts.
_BANNER_ALREADY_PRINTED = False

ASCII_LOGO = r"""
        __  _________  _________
       / / / / __/ _ \/ ___/ __/
      / /_/ / _// // / (_ / _/ 
      \____/___/____/\___/___/
              POWERED BY 
   __  ________________  ____  __   ____
  / / / / __/_  __/ __ \/ __ \/ /  / __/
 / /_/ / _/  / / / /_/ / /_/ / /___\ \  
 \____/___/ /_/  \____/\____/____/___/  
""".strip("\n")


def _check_newer_uetools_ver():
    pkg = "uedgetools"

    try:
        import json
        import urllib.request

        try:
            import importlib.metadata

            thisver = importlib.metadata.version(pkg)
        except:
            import pkg_resources

            thisver = pkg_resources.get_distribution(pkg).version

        contents = urllib.request.urlopen("https://pypi.org/pypi/" + pkg + "/json").read()
        data = json.loads(contents.decode())
        thatver = data["info"]["version"]

        if thisver < thatver:
            return f"\nAn update to UETOOLS v{thatver} is available via PyPi (pip)"+ \
                    "\nTo update: pip install uedgetools"
    except Exception as err:
        return "\nError checking pypi version: {}".format(err)


def _env_flag(name: str, default: str = "1") -> bool:
    """
    Interpret env var as boolean-ish flag.
    Returns True if enabled, False if explicitly disabled.
    We treat "0", "false", "no", "off" (case-insensitive) as disabled.
    """
    val = environ.get(name, default)
    return str(val).strip().lower() not in {"0", "false", "no", "off"}


def _get_dist_version(dist_name: str) -> str | None:
    """
    Return installed distribution version, or None if not installed / not found.
    """
    try:
        from importlib.metadata import version as dist_version  # py3.8+
    except Exception:
        try:
            from importlib_metadata import version as dist_version  # backport
        except Exception:
            return None

    try:
        return dist_version(dist_name)
    except Exception:
        return None


def _print_banner_once() -> None:
    global _BANNER_ALREADY_PRINTED

    # Hard stop if banner explicitly disabled.
    if not _env_flag("UEDGE_BANNER", default="1"):
        return

    # Only print once per Python session.
    if _BANNER_ALREADY_PRINTED:
        return
    _BANNER_ALREADY_PRINTED = True

    show_logo = _env_flag("UEDGE_LOGO", default="1")
    show_versions = _env_flag("UEDGE_VERSION", default="1")

    lines: list[str] = []

    if show_logo and ASCII_LOGO:
        lines.append(ASCII_LOGO)

    if show_versions:
        # Adjust these distribution names if your install uses different names.
        uedge_ver = _get_dist_version("uedge") or "unknown"
        uetools_ver = __version__


        lines.append("Versions:")
        lines.append(f"  uedge  : {uedge_ver}")
        lines.append(f"  uetools: {uetools_ver}")
        ver = _check_newer_uetools_ver()
        if ver is not None:
            lines.append(ver)
    if lines:
        print("\n".join(lines))


# Print banner at import time (but only once per session)
_print_banner_once()

# ---------- Defaults to match UEDGE-style behavior ----------
environ.setdefault("UEDGE_LOGO", "0")
environ.setdefault("UEDGE_VERSION", "0")
## Optional: expose version programmatically
#__version__ = _get_dist_version("uetools") or "unknown"

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
config = Config()
config.case()
