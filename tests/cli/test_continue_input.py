import pytest
import argparse
import os

from filecmp import cmp as filecmp
from unittest.mock import patch

from PQAnalysis.cli.continue_input import main
from PQAnalysis.exceptions import PQNotImplementedError



@pytest.mark.parametrize(
    "example_dir",
    ["inputFileReader/PQ_input"],
    indirect=False
)
def test_continue_input(test_with_data_dir, capsys):

    with pytest.raises(PQNotImplementedError) as exception:
        with patch('argparse._sys.argv',
            ['continue_input.py',
            'input.in',
            '-n',
            '1',
            '--input-format',
            'qmcfc']):
            main()
    assert str(
        exception.value
    ) == "Format InputFileFormat.QMCFC not implemented yet for continuing input file."

    with patch('argparse._sys.argv',
        ['continue_input.py',
        'run-08.in',
        '-n',
        '2']):
        main()

    assert filecmp("run-09.in", "run-09.in.ref")
    assert filecmp("run-10.in", "run-10.in.ref")

    with patch('argparse._sys.argv',
        ['continue_input.py',
        'run-08.rpmd.in',
        '-n',
        '3',
        '--input-format',
        'PQ']):
        main()

    assert filecmp("run-09.rpmd.in", "run-09.rpmd.in.ref")
    assert filecmp("run-10.rpmd.in", "run-10.rpmd.in.ref")
    assert os.path.exists("run-11.rpmd.in")
