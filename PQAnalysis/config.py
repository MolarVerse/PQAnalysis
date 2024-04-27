"""
A module defining configurations of the PQAnalysis package.
"""
# pylint: disable=invalid-name

base_url = "https://molarverse.github.io/PQAnalysis/"
code_base_url = base_url + "code/"

PQ_docs_url = "https://molarverse.github.io/PQ/"

with_progress_bar = True
use_log_file = False  # defined via environment variable and __init__.py of PQAnalysis
log_file_name = None  # defined via environment variable and __init__.py of PQAnalysis
