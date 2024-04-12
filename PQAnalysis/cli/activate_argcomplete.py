import os

from pathlib import Path


def main():
    __dir__ = Path(__file__).parent

    files = []
    for f in __dir__.iterdir():
        if f.is_file() and f.suffix == ".py" and not str(f).split("/")[-1].startswith("_") and "activate_argcomplete" not in str(f):
            files.append(f.stem)

    for file in files:
        print(file)
        command = f"eval $(register-python-argcomplete {file})"
        print(command)
        os.system(command)
