export PQANALYSIS_BEARTYPE_LEVEL=DEBUG
#check if $@ is empty
if [ -z "$@" ]; then
    python setup.py pytest
else
    python setup.py pytest --addopts $@
fi
if [ $? -ne 0 ]; then
    exit 1
fi

export PQANALYSIS_BEARTYPE_LEVEL=RELEASE
#check if $@ is empty
if [ -z "$@" ]; then
    python setup.py pytest
else
    python setup.py pytest --addopts $@
fi
