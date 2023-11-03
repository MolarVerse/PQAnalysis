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


@pytest.fixture(scope="function")
def test_with_data_dir(example_dir):

    tmpdir = "tmpdir"

    if os.path.exists(tmpdir) and os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)

    shutil.copytree(os.path.join("tests/data", example_dir), tmpdir)

    os.chdir(tmpdir)

    yield tmpdir

    os.chdir("..")
    shutil.rmtree(tmpdir)
