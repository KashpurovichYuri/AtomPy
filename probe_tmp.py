import sys
sys.path.insert(0, "src"); sys.path.insert(0, "tests")
import numpy as np
from conftest import _frozen_mf  # no: conftest import by name fails; build inline
