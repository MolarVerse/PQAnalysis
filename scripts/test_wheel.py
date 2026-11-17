"""Test an installed PQAnalysis wheel against the repository test suite.

Usage: python scripts/test_wheel.py PROJECT_DIR

The tests run from a copy outside PROJECT_DIR, so ``import PQAnalysis``
loads the installed wheel and not the source tree.
"""

import importlib
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

EXTENSIONS = (
    "PQAnalysis.io.traj_file.process_lines",
    "PQAnalysis.io.traj_file._slab_parser",
    "PQAnalysis.analysis.momentum._momentum_kernel",
    "PQAnalysis.analysis.msd._msd_kernel",
    "PQAnalysis.analysis.rdf._rdf_kernel",
    "PQAnalysis.analysis.vacf._vacf_kernel",
)

COMMANDS = ("pqanalysis", "rdf", "msd", "vacf", "traj2com")



def check_installed(project):
    """Fail if a compiled module is missing or comes from the source tree."""
    for name in EXTENSIONS:
        path = Path(importlib.import_module(name).__file__).resolve()
        if project in path.parents or path.suffix not in (".so", ".pyd"):
            sys.exit(f"{name} is not a compiled module of the wheel: {path}")
        print(f"ok {name}: {path.name}")



def check_commands():
    """Fail if a console script of the wheel does not start."""
    for command in COMMANDS:
        executable = shutil.which(command)
        if executable is None:
            sys.exit(f"{command} is not on PATH")
        subprocess.run(
            [executable, "--help"], check=True, stdout=subprocess.DEVNULL
        )
        print(f"ok {command} --help")



def run_tests(project):
    """Run the test suite from a copy; type-check levels are tested in ci."""
    with tempfile.TemporaryDirectory() as work:
        for name in ("tests", "examples"):
            shutil.copytree(project / name, Path(work) / name)
        shutil.copy(project / "pytest.ini", work)

        env = dict(os.environ, PQANALYSIS_BEARTYPE_LEVEL="RELEASE")
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "pytest",
                "tests",
                "-q",
                "-o",
                "addopts=",
                "-p",
                "no:cacheprovider"
            ],
            cwd=work,
            env=env,
            check=False,
        )
        if result.returncode != 0:
            sys.exit(result.returncode)



def main():
    """Check the compiled modules and commands, then run the tests."""
    project = Path(sys.argv[1]).resolve()
    check_installed(project)
    check_commands()
    run_tests(project)



if __name__ == "__main__":
    main()
