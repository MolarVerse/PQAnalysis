"""Tests for deferred progress-bar backend selection."""

import sys

from types import ModuleType

from PQAnalysis.utils import progress



def test_terminal_backend_is_selected_outside_notebooks(monkeypatch):
    calls = []
    expected = object()
    terminal_backend = ModuleType("tqdm")

    def terminal_tqdm(*args, **kwargs):
        calls.append((args, kwargs))
        return expected

    terminal_backend.tqdm = terminal_tqdm
    monkeypatch.delitem(sys.modules, "ipykernel", raising=False)
    monkeypatch.setitem(sys.modules, "tqdm", terminal_backend)

    actual = progress.tqdm("frames", total=3)

    assert actual is expected
    assert calls == [(("frames",), {"total": 3})]



def test_notebook_backend_is_selected_when_ipykernel_is_loaded(monkeypatch):
    calls = []
    expected = object()
    notebook_backend = ModuleType("tqdm.auto")

    def notebook_tqdm(*args, **kwargs):
        calls.append((args, kwargs))
        return expected

    notebook_backend.tqdm = notebook_tqdm
    monkeypatch.setitem(sys.modules, "ipykernel", ModuleType("ipykernel"))
    monkeypatch.setitem(sys.modules, "tqdm.auto", notebook_backend)

    actual = progress.tqdm("frames", disable=True)

    assert actual is expected
    assert calls == [(("frames",), {"disable": True})]
