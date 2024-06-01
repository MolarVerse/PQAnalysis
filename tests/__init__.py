"""
Unit tests for the PQAnalysis package.
"""

import os

__beartype_level__ = os.environ.get("PQANALYSIS_BEARTYPE_LEVEL", "RELEASE")
__beartype_level__ = __beartype_level__.upper()

os.environ['CONSTANT_SEED'] = '42'
