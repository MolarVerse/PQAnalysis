"""Regression tests for deferred package and CLI imports."""

import importlib
import json
import subprocess
import sys


PACKAGE_NAMES = (
    "PQAnalysis.analysis",
    "PQAnalysis.core",
    "PQAnalysis.io",
    "PQAnalysis.io.input_file_reader",
    "PQAnalysis.io.traj_file",
    "PQAnalysis.topology",
    "PQAnalysis.traj",
    "PQAnalysis.utils",
)



def _run_python(source: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-c", source],
        check=True,
        capture_output=True,
        text=True,
    )



def _run_cli_and_list_loaded(arguments: list[str]) -> dict:
    completed = _run_python(
        f"""
import contextlib
import io
import json
import sys
import PQAnalysis.cli.main as cli_main

sys.argv = ["pqanalysis", *{arguments!r}]
try:
    with contextlib.redirect_stdout(io.StringIO()):
        with contextlib.redirect_stderr(io.StringIO()):
            cli_main.main()
except SystemExit as exception:
    exit_code = exception.code
else:
    exit_code = 0

command_modules = [
    f"PQAnalysis.cli{{module_name}}"
    for module_name, _, _ in cli_main._COMMANDS.values()
]
print(json.dumps({{
    "exit_code": exit_code,
    "loaded": sorted(name for name in command_modules if name in sys.modules),
}}))
"""
    )

    return json.loads(completed.stdout)



def test_cli_main_import_is_silent_and_does_not_load_commands():
    completed = _run_python(
        """
import json
import sys
import PQAnalysis.cli.main as cli_main

command_modules = [
    f"PQAnalysis.cli{module_name}"
    for module_name, _, _ in cli_main._COMMANDS.values()
]
print(json.dumps(sorted(name for name in command_modules if name in sys.modules)))
"""
    )

    assert completed.stdout == "[]\n"
    assert completed.stderr == ""



def test_root_help_does_not_load_deferred_commands():
    assert _run_cli_and_list_loaded(["--help"]) == {
        "exit_code": 0,
        "loaded": [],
    }



def test_subcommand_help_loads_only_selected_cli_module():
    assert _run_cli_and_list_loaded(["rdf", "--help"]) == {
        "exit_code": 0,
        "loaded": ["PQAnalysis.cli.rdf"],
    }



def test_analysis_star_import_preserves_function_exports():
    completed = _run_python(
        """
import json

namespace = {}
exec("from PQAnalysis.analysis import *", namespace)
print(json.dumps({name: callable(namespace[name]) for name in ("msd", "rdf", "vacf")}))
"""
    )

    assert json.loads(completed.stdout) == {
        "msd": True,
        "rdf": True,
        "vacf": True,
    }



def test_direct_analysis_submodule_imports_do_not_shadow_api_functions():
    completed = _run_python(
        """
import importlib
import json

analysis = importlib.import_module("PQAnalysis.analysis")
for name in ("msd", "rdf", "vacf"):
    importlib.import_module(f"PQAnalysis.analysis.{name}")

print(json.dumps({
    name: callable(getattr(analysis, name))
    for name in ("msd", "rdf", "vacf")
}))
"""
    )

    assert json.loads(completed.stdout) == {
        "msd": True,
        "rdf": True,
        "vacf": True,
    }



def test_analysis_import_does_not_initialize_progress_backend():
    completed = _run_python(
        """
import json
import sys
from PQAnalysis.analysis.momentum import Momentum

print(json.dumps({
    "momentum": Momentum.__name__,
    "tqdm": "tqdm" in sys.modules,
    "tqdm_auto": "tqdm.auto" in sys.modules,
}))
"""
    )

    assert json.loads(completed.stdout) == {
        "momentum": "Momentum",
        "tqdm": False,
        "tqdm_auto": False,
    }



def test_all_declared_lazy_exports_resolve():
    completed = _run_python(
        f"""
import importlib
import json

package_names = {PACKAGE_NAMES!r}
failures = []
for package_name in package_names:
    package = importlib.import_module(package_name)
    for name in package.__all__:
        try:
            getattr(package, name)
        except (AttributeError, ImportError) as exception:
            failures.append((package_name, name, type(exception).__name__))

print(json.dumps(failures))
"""
    )

    assert json.loads(completed.stdout) == []



def test_declared_lazy_exports_are_visible_to_dir():
    for package_name in PACKAGE_NAMES:
        package = importlib.import_module(package_name)

        assert set(package.__all__).issubset(dir(package))
