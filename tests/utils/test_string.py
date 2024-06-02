"""
Test cases for the string utility functions.
"""

from PQAnalysis.utils.string import is_comment_line

from . import pytestmark  # pylint: disable=unused-import



def test_is_comment_line():
    """
    Test the is_comment_line function.
    """

    assert not is_comment_line("This is not a comment line.")

    assert is_comment_line("# This is a comment line.")
    assert not is_comment_line(" ")
    assert is_comment_line(" ", empty_line=True)

    assert not is_comment_line("# This is a comment line.", comment_token="!")
    assert is_comment_line("! This is a comment line.", comment_token="!")
    assert is_comment_line("     !", comment_token="!", empty_line=False)
