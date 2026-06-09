import os
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_release_type_checking_error_has_no_stacktrace():
    env = os.environ.copy()
    env["PQANALYSIS_BEARTYPE_LEVEL"] = "RELEASE"
    env["PYTHONPATH"] = os.pathsep.join(
        [str(PROJECT_ROOT), env.get("PYTHONPATH", "")]
    )

    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "from PQAnalysis.analysis import RDF\nRDF(1, 1, 1)\n",
        ],
        cwd=PROJECT_ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 1
    assert "PQAnalysis.TypeChecking - PQTypeError" in result.stderr
    assert "Argument 'traj'" in result.stderr
    assert "Traceback" not in result.stderr
    assert "BeartypeCallHintParamViolation" not in result.stderr
