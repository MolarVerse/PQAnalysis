import pytest

from unittest import mock

from PQAnalysis.cli.xyz2gen import main, XYZ2GENCLI
from PQAnalysis.cli._argument_parser import _ArgumentParser



def test_name():
    assert XYZ2GENCLI.program_name() == "xyz2gen"



@pytest.mark.parametrize(
    "value, expected",
    [
    ("True", True),
    ("False", False),
    ("None", None),
    ("true", True),
    ("false", False),
    ("none", None),
    ]
)
def test_periodic_argument_parsing(value, expected):
    args = _parse(["md-01.xyz", "--periodic", value])

    assert args.periodic is expected



def test_periodic_argument_default_and_invalid(capsys):
    args = _parse(["md-01.xyz"])
    assert args.periodic is None

    with pytest.raises(SystemExit):
        _parse(["md-01.xyz", "--periodic", "maybe"])

    captured = capsys.readouterr()
    assert "invalid choice: 'maybe'" in captured.err



def _parse(argv):
    parser = _ArgumentParser(description="test")
    XYZ2GENCLI.add_arguments(parser)
    return parser.parse_args(argv)



@pytest.mark.parametrize("example_dir", ["xyz2rst"], indirect=False)
def test_main_periodic(test_with_data_dir):
    with mock.patch(
        "sys.argv",
        ["xyz2gen", "md-01.xyz", "--periodic", "True", "-o", "box.gen"],
    ):
        main()

    with open("box.gen", "r", encoding="utf-8") as file:
        box_lines = file.read().splitlines()

    assert box_lines[0].split() == ["4", "S"]
    assert len(box_lines) == 10  # header + elements + 4 atoms + 4 cell lines

    with mock.patch(
        "sys.argv",
        ["xyz2gen", "md-01.xyz", "--periodic", "False", "-o", "nobox.gen"],
    ):
        main()

    with open("nobox.gen", "r", encoding="utf-8") as file:
        nobox_lines = file.read().splitlines()

    assert nobox_lines[0].split() == ["4", "C"]
    assert len(nobox_lines) == 6  # header + elements + 4 atoms, no cell
