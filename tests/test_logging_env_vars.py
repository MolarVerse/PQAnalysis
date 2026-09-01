import os
import subprocess
import sys

from pathlib import Path

import PQAnalysis

PACKAGE_ROOT = Path(PQAnalysis.__file__).resolve().parents[1]


def run_import(log_file_value):
    env = os.environ.copy()
    env.pop("PQANALYSIS_LOGGING_LEVEL", None)
    env["PQANALYSIS_LOG_FILE"] = log_file_value
    env["PYTHONPATH"] = str(PACKAGE_ROOT)

    return subprocess.run(
        [
            sys.executable,
            "-c",
            "import PQAnalysis.config as config; import PQAnalysis; "
            "print(config.use_log_file, config.log_file_name)",
        ],
        cwd=PACKAGE_ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )


def test_log_file_env_var_without_logging_level():
    result = run_import("on")

    assert result.returncode == 0, result.stderr

    use_log_file, log_file_name = result.stdout.split()
    assert use_log_file == "True"
    assert log_file_name.startswith("PQAnalysis_")


def test_log_file_env_var_off():
    result = run_import("off")

    assert result.returncode == 0, result.stderr

    use_log_file, log_file_name = result.stdout.split()
    assert use_log_file == "False"
    assert log_file_name != "off"


def test_log_file_env_var_custom_name():
    result = run_import("my_custom.log")

    assert result.returncode == 0, result.stderr
    assert result.stdout.split() == ["True", "my_custom.log"]
