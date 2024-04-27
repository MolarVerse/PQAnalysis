"""
A subpackage to handle .gen files.
"""

from .gen_file_reader import GenFileReader
from .gen_file_writer import GenFileWriter

from .api import read_gen_file, write_gen_file
