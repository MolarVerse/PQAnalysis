"""
Tests for the partial output of add_molecule.
"""

import logging
import os
import signal
import threading
import time

import pytest

from PQAnalysis import __package_name__
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.atomic_system.exceptions import AtomicSystemError
from PQAnalysis.tools.add_molecule import (
    AddMolecule,
    add_molecule,
    _sigterm_as_keyboard_interrupt,
)

from ..conftest import caplog_for_logger

FRAMEWORK = """Box 12 12 12 90 90 90
C 0 0 0.0 0.0 0.0 0 0 0 0 0 0
C 0 0 3.0 0.0 0.0 0 0 0 0 0 0
"""

WATER = """3

O 0 0 0
H 0.7572 0.5865 0
H -0.7572 0.5865 0
"""

FRAMEWORK_TOP = """SHAKE 1 1 0
1 2 1.0
END
"""

WATER_TOP = """SHAKE 1 1 0
1 2 1.0
END
"""



@pytest.fixture
def inputs(tmpdir, monkeypatch):  # pylint: disable=unused-argument
    monkeypatch.setenv("CONSTANT_SEED", "5")
    with open("framework.rst", "w", encoding="utf-8") as file:
        file.write(FRAMEWORK)
    with open("water.xyz", "w", encoding="utf-8") as file:
        file.write(WATER)
    with open("framework.top", "w", encoding="utf-8") as file:
        file.write(FRAMEWORK_TOP)
    with open("water.top", "w", encoding="utf-8") as file:
        file.write(WATER_TOP)



def stop_after(monkeypatch, n_placed, exception):
    original = AtomicSystem._fit_atomic_system
    calls = {"n": 0}

    def fit(self, *args, **kwargs):
        if calls["n"] == n_placed:
            raise exception
        calls["n"] += 1
        return original(self, *args, **kwargs)

    monkeypatch.setattr(AtomicSystem, "_fit_atomic_system", fit)



def atom_lines(path):
    with open(path, "r", encoding="utf-8") as file:
        return [line for line in file if not line.startswith("Box")]



class TestPartialOutput:

    def test_complete_run(self, inputs):  # pylint: disable=unused-argument
        add_molecule(
            "framework.rst",
            "water.xyz",
            output_file="out.rst",
            number_of_additions=4,
        )

        assert len(atom_lines("out.rst")) == 2 + 4 * 3

    @pytest.mark.parametrize(
        "exception",
        [AtomicSystemError("max iterations"), KeyboardInterrupt()],
    )
    def test_writes_placed_molecules(
        self, inputs, monkeypatch, caplog, exception
    ):  # pylint: disable=unused-argument
        stop_after(monkeypatch, 3, exception)

        with caplog_for_logger(
            caplog,
            __package_name__ + "." + AddMolecule.__qualname__,
            logging.WARNING
        ):
            with pytest.raises(type(exception)):
                add_molecule(
                    "framework.rst",
                    "water.xyz",
                    output_file="out.rst",
                    number_of_additions=10,
                    topology_file="framework.top",
                    topology_file_to_add="water.top",
                    topology_file_output="out.top",
                )

        lines = atom_lines("out.rst")
        assert len(lines) == 2 + 3 * 3
        assert [line.split()[1] for line in lines[2:]
                ] == [str(i) for i in (1, 1, 1, 2, 2, 2, 3, 3, 3)]

        messages = " ".join(record.getMessage() for record in caplog.records)
        assert "3/10" in messages

        with open("out.top", "r", encoding="utf-8") as file:
            assert file.read().split()[1] == "4"

    def test_nothing_placed_writes_nothing(self, inputs, monkeypatch):  # pylint: disable=unused-argument
        stop_after(monkeypatch, 0, AtomicSystemError("max iterations"))

        with pytest.raises(AtomicSystemError):
            add_molecule(
                "framework.rst",
                "water.xyz",
                output_file="out.rst",
                number_of_additions=2,
            )

        assert not os.path.exists("out.rst")



class TestSigterm:

    def test_sigterm_raises_keyboard_interrupt(self):
        previous = signal.getsignal(signal.SIGTERM)

        with pytest.raises(KeyboardInterrupt):
            with _sigterm_as_keyboard_interrupt():
                os.kill(os.getpid(), signal.SIGTERM)
                time.sleep(5)

        assert signal.getsignal(signal.SIGTERM) is previous

    def test_other_threads_are_ignored(self):
        previous = signal.getsignal(signal.SIGTERM)
        seen = []

        def run():
            with _sigterm_as_keyboard_interrupt():
                seen.append(signal.getsignal(signal.SIGTERM))

        thread = threading.Thread(target=run)
        thread.start()
        thread.join()

        assert seen == [previous]
