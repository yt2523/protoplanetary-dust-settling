"""Add ../src to sys.path so the test scripts can import the solver modules.

Import this at the top of any test script:

    import _path_setup  # noqa: F401
"""

import os
import sys

_SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src")
_SRC = os.path.normpath(_SRC)

if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
