export PQANALYSIS_BEARTYPE_LEVEL=DEBUG
python setup.py pytest --addopts $@
if [ $? -ne 0 ]; then
    exit 1
fi

export PQANALYSIS_BEARTYPE_LEVEL=RELEASE
python setup.py pytest --addopts $@
