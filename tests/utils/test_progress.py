"""Tests for deferred progress-bar backend selection."""

import sys

from types import ModuleType

from PQAnalysis.utils import progress



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
    assert calls == [(('frames',), {"disable": True})]
