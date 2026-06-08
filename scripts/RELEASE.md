# Release Checklist

Use this order for a normal OCELOT release. The examples below use `26.06.1`;
replace it with the intended version.

## 1. Decide the version

OCELOT uses versions like `26.06.0`. A bugfix after `26.06.0` should be
`26.06.1`.

Update `CHANGELOG.md` before running the release scripts. The scripts update
version metadata, but they do not write release notes.

## 2. Prepare from a clean `dev` branch

```bash
git checkout dev
git pull
git status
python scripts/prepare_release.py 26.06.1 --full --dry-run
python scripts/prepare_release.py 26.06.1 --full
```

`prepare_release.py` updates version metadata, runs tests and smoke demos,
builds package artifacts, commits the release metadata, merges `dev` into
`master`, and creates the `v26.06.1` tag. It does not push or publish unless
explicit flags are passed.

After checking the result:

```bash
git push origin dev master --tags
```

## 3. Publish artifacts explicitly

Run this first to rebuild and check artifacts without uploading:

```bash
python scripts/publish_release.py 26.06.1
```

After credentials are configured and the artifacts look correct:

```bash
python scripts/publish_release.py 26.06.1 --upload-pypi --upload-conda
```

PyPI upload uses Twine credentials. Anaconda upload uses `anaconda-client`;
run `anaconda login` first if needed.

## 4. Optional GitHub release

If a GitHub release is needed, either create it manually from tag `v26.06.1` or
run `prepare_release.py` with `--github-release` during the preparation step.

## Notes

- Start from a clean worktree. Avoid `--allow-dirty` for normal releases.
- Use `--skip-tests`, `--skip-demos`, or `--skip-build` only for an intentional
  emergency workflow.
- `publish_release.py` expects `HEAD` to be tagged as `v<version>` unless
  `--allow-tag-mismatch` is passed.
