"""The u14e V&V family's library: the trusted reference, the windows,
the routes, the pricing, the table.  The legs (``u14e_vv_*.py``) are
its readers -- each leg a case list and its checks; the probes
(``probe_u14e_*.py``) print only.  Importing the package puts the
repository's ``src`` on the path, so the legs run as plain scripts."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))
