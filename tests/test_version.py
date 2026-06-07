import builtins
import importlib



def test_version_fallback(monkeypatch):
    import PQAnalysis.version as version

    real_import = builtins.__import__

    def fail_version_import(
        name, globals=None, locals=None, fromlist=(), level=0
    ):
        if level == 1 and name == "_version":
            raise ModuleNotFoundError("No module named 'PQAnalysis._version'")

        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", fail_version_import)

    reloaded_version = importlib.reload(version)

    assert reloaded_version.__version__ == "0+unknown"
    assert reloaded_version.__version_tuple__ == (0, "unknown")

    monkeypatch.undo()
    importlib.reload(version)
