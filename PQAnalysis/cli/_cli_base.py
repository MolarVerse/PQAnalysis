"""
This module contains the abstract base class for all the CLI classes.
"""

from abc import ABCMeta, abstractmethod



class CLIBase(metaclass=ABCMeta):

    """
    Abstract base class for all the CLI classes.
    """

    @classmethod
    @abstractmethod
    def program_name(cls) -> str:
        """
        Return the name of the program.

        Returns
        -------
        str
            The name of the program.
        """

    @classmethod
    @abstractmethod
    def add_arguments(cls, parser):
        """
        Add the arguments to the parser.

        Parameters
        ----------
        parser : _type_
            _description_
        """

    @classmethod
    @abstractmethod
    def run(cls, args):
        """
        Run the CLI.
        """
