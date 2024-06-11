"""
PQAnalysis is a Python package for the analysis of molecular topologies and trajectories.
"""
# pylint: disable=invalid-name

import logging
import os
import time

from pathlib import Path
from beartype.claw import beartype_this_package

import PQAnalysis.config as config  # pylint: disable=consider-using-from-import

from PQAnalysis.utils.custom_logging import CustomLogger

__base_path__ = Path(__file__).parent
__package_name__ = __name__

##################
# BEARTYPE SETUP #
##################

__beartype_default_level__ = "RELEASE"
__beartype_level__ = os.getenv(
    "PQANALYSIS_BEARTYPE_LEVEL", __beartype_default_level__
)

if __beartype_level__.upper() == "DEBUG":
    beartype_this_package()

#################
# LOGGING SETUP #
#################

logging_env_var = os.getenv("PQANALYSIS_LOGGING_LEVEL")

if logging_env_var and logging_env_var not in logging.getLevelNamesMapping():
    raise ValueError(
        f"Invalid logging level: {logging_env_var}. Valid logging "
        f"levels are: {logging.getLevelNamesMapping()}"
    )

if (
    'execution_start_time' not in vars() and
    'execution_start_time' not in globals()
):
    execution_start_time = time.strftime('%Y-%m-%d_%H-%M-%S', time.localtime())

logging.setLoggerClass(CustomLogger)
logging.basicConfig(level=os.getenv("PQANALYSIS_LOGGING_LEVEL", "INFO"))
package_logger = logging.getLogger(__name__)

log_file_env_var = os.getenv("PQANALYSIS_LOG_FILE")

if log_file_env_var and logging_env_var.lower() != "off":
    config.use_log_file = True

    if log_file_env_var.lower() != "on" and len(log_file_env_var) > 0:
        config.log_file_name = log_file_env_var

if config.log_file_name is None:
    config.log_file_name = f"PQAnalysis_{execution_start_time}.log"
