export PQANALYSIS_BEARTYPE_LEVEL=DEBUG
python -m pytest $@
if [ $? -ne 0 ]; then
    exit 1
fi

export PQANALYSIS_BEARTYPE_LEVEL=RELEASE
python -m pytest $@
