"""
PQAnalysis is a Python package for the analysis of molecular topologies and trajectories.
"""

import logging

from pathlib import Path
from beartype.claw import beartype_this_package

beartype_this_package()

__base_path__ = Path(__file__).parent

package_name = __name__

package_logger = logging.getLogger(__name__)
package_logger.setLevel(logging.INFO)
