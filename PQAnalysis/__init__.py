"""
PQAnalysis is a Python package for the analysis of molecular topologies and trajectories.
"""

from pathlib import Path
from beartype.claw import beartype_this_package

beartype_this_package()

__base_path__ = Path(__file__).parent
