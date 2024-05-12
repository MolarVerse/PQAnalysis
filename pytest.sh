export PQANALYSIS_BEARTYPE_LEVEL=DEBUG
python -m pytest $@
export PQANALYSIS_BEARTYPE_LEVEL=RELEASE
python -m pytest $@
