# conda-forge recipe

Draft [conda-forge](https://conda-forge.org/) recipe for distributing PQAnalysis
through the `conda-forge` channel. It is kept here for reference and
maintenance; the recipe conda-forge actually builds lives in the `pqanalysis`
*feedstock* created from
[`conda-forge/staged-recipes`](https://github.com/conda-forge/staged-recipes).

PQAnalysis ships a compiled Cython extension
(`PQAnalysis/io/traj_file/process_lines.pyx`), so this is a per-platform
(non-`noarch`) recipe: it needs a C compiler (`{{ compiler('c') }}` /
`{{ stdlib('c') }}`) and `build.sh` / `bld.bat`. It requires Python `>=3.12`
(`skip: true  # [py<312]`). All run dependencies (`numpy`, `scipy`, `beartype`,
`multimethod`, `lark`, `tqdm`, `decorator`, `argcomplete`, `rich-argparse`) are
already on conda-forge.

## Before submitting

- **Maintainer(s):** `extra.recipe-maintainers` lists `galjos`. Add any other
  GitHub usernames who should co-maintain the feedstock.
- **Version & hash** are pinned to the current PyPI release (1.3.0). To refresh
  for a new release, bump `version` and replace `sha256` with the sdist hash:

  ```bash
  curl -sL https://pypi.org/pypi/PQAnalysis/json \
    | python -c "import json,sys; d=json.load(sys.stdin); \
        print(next(u['digests']['sha256'] for u in d['urls'] if u['packagetype']=='sdist'))"
  ```

  After the feedstock exists, conda-forge's `regro-cf-autotick-bot` opens
  version-bump PRs automatically.

## Local check (optional)

With `conda-build` / `conda-smithy` installed:

```bash
conda smithy recipe-lint conda-recipes/pqanalysis
conda build conda-recipes/pqanalysis -c conda-forge
```

See the conda-forge [contributing guide](https://conda-forge.org/docs/maintainer/adding_pkgs/).
