import pytest
import os
import shutil


@pytest.fixture(scope="function")
def tmpdir():

    tmpdir = "tmpdir"

    if os.path.exists(tmpdir) and os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    os.mkdir(tmpdir)

    os.chdir(tmpdir)

    yield tmpdir

    os.chdir("..")
    shutil.rmtree(tmpdir)
