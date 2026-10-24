# Development checks and releases

Use `python -m invoke` (or `inv`) from the repository root. Install the development
extras in your working environment: `python -m pip install -e '.[dev]'`.

## Check development without preparing a release

```
python -m invoke run-tests
python -m invoke run-demos --timeout=300
```

The demo runner attempts every `.py` and `.ipynb` under `demos`, excluding notebook
checkpoints. It copies the demo tree to `tmp/demo-check-*`, runs each example from
its own directory, and writes logs, executed notebooks and `results.json` there.
A failure or timeout makes the command exit nonzero after completing the sweep.
Examples requiring Genesis, MPI, external data, or site-specific paths may fail
on a workstation: inspect their logs and record the missing prerequisites.
A timeout means incomplete validation, not a physics error. Increase `--timeout`
or rerun a subset with `--pattern=rk` (substring of the relative path).
Use `--exclude-scripts`, `--exclude-notebooks`, or `--critical-demos-only` to select
work. `--dry-run` lists selected examples without running them.

On macOS with conflicting OpenMP runtimes, set `NUMBA_THREADING_LAYER=workqueue`
as described in README.md. Matplotlib uses Agg by default. The runner uses your
current Python environment; optional demo dependencies must be installed there.

## Prepare a release later

Versioning is `YY.MM.patch`, for example `26.10.0`. Substitute the chosen version
in the commands below. Start on an up-to-date, clean `dev` branch.

1. Review changes since the previous tag, resolve relevant Dependabot alerts, and
   finish the unreleased changelog. Rename that section to the version/date and
   start a new empty Unreleased section.
2. Run the full tests and demos. Investigate failures; do not call a timed-out or
   externally blocked example a pass.
3. Update all version metadata and commit it:

   ```
   python -m invoke bump-version --version=26.10.0
   git diff
   python -m invoke commit-release --version=26.10.0
   python -m invoke build
   ```

   The version task also changes the Conda recipe's source tag. Test installation
   of the built wheel in a fresh environment, outside the source directory.
   Each build writes to a new `dist/build-*` directory and checks only its own
   artifacts; previous builds are kept.
4. Review/merge `dev` into `master` (a pull request is convenient), then tag the
   tested release commit. For a local merge:

   ```
   python -m invoke merge --version=26.10.0
   python -m invoke tag --version=26.10.0
   python -m invoke push --version=26.10.0
   ```

   Push includes only the specified release tag, not all local tags. The legacy
   `do-full-release` task performs local tests, demos, version update, commit,
   build, merge and tag; it now stops before any push or upload. Explicit steps
   are easier to review and resume.

## Publish explicitly

Publishing requires clean source at `v<version>` and matching Python/Conda
metadata. It builds into a new directory under `dist/release-<version>` and runs
Twine's package checks. No upload occurs by default:

```
python -m invoke publish --version=26.10.0
```

Once ready, explicitly select the destination:

```
python -m invoke publish --version=26.10.0 --upload-pypi --repository=testpypi
python -m invoke publish --version=26.10.0 --upload-pypi
```

Each invocation builds fresh artifacts from the same tagged source. TestPyPI and
PyPI use separate credentials. Configure Twine credentials using a keyring or
its documented environment variables; never commit tokens. Follow the
[Python packaging guide](https://packaging.python.org/en/latest/tutorials/packaging-projects/).

For Conda, install `conda-build` and `anaconda-client` in your release environment.
The recipe downloads its source from the release tag, so push the tag first.
Build without uploading, then explicitly publish using your configured account:

```
python -m invoke publish --version=26.10.0 --conda
anaconda login
python -m invoke publish --version=26.10.0 --upload-conda
```

Conda builds use conda-forge and disable automatic upload. All requested builds
and checks complete before either upload starts. PyPI and Anaconda are separate
services: a failure during the second upload cannot undo the first one.
Create a GitHub release from the same tag and use the changelog as release notes.

Keeping these local commands is the simplest next step. Later, GitHub Actions
with PyPI Trusted Publishing can replace local PyPI credentials, after a separate
review of the workflow and repository/PyPI permissions.
