#!/usr/bin/env python3
"""
Common tasks for developers only.
"""

from __future__ import annotations

from invoke import task
import os
import re
import subprocess
import sys
import shutil
import tempfile
import json
import signal
import time
import tomllib
from pathlib import Path

VERBOSE = True

SOURCE_BRANCH = "dev"
TARGET_BRANCH = "master"

CRITICAL_DEMOS = (
    "demos/ebeam/rk_track.py",
    "demos/ebeam/dba.py",
    "demos/ebeam/dba_tracking.py",
)

VERSION_FILES = (
    ("pyproject.toml", r"version\s*=\s*['\"][^'\"]+['\"]", r"version = '{version}'"),
    (
        "ocelot/__init__.py",
        r"__version__\s*=\s*['\"][^'\"]+['\"]",
        r"__version__ = '{version}'",
    ),
    (
        "conda-recipe/meta.yaml",
        r'version:\s*"[^"]+"',
        r'version: "{version}"',
    ),
)

BUILD_DIR = Path("dist")


def repo_root() -> Path:
    return Path(__file__).resolve().parent


def dev_env() -> dict[str, str]:
    env = os.environ.copy()
    # Set MPLBACKEND to Agg to prevent matplotlib from trying to open a window and blocking the execution of the script until manual intervention.
    env.setdefault("MPLBACKEND", "Agg")
    env["PYTHONPATH"] = str(repo_root()) + os.pathsep + env.get("PYTHONPATH", "")
    return env


def run(
    command: list[str],
    *,
    cwd: Path,
    env: dict[str, str] | None = None,
    dry_run: bool,
    supress_stdout: bool = False,
) -> subprocess.CompletedProcess:
    """
    Run a command in a subprocess.

    Return:
        Either a subprocess.CompletedProcess of the run command or a subprocess.CompletedProcess with returncode 0 if `dry_run` is True.

    Raises:
        subprocess.CalledProcessError if the command fails.
    """

    if VERBOSE:
        print("+ " + " ".join(command))

    if dry_run:
        return subprocess.CompletedProcess(command, 0)
    return subprocess.run(
        command,
        cwd=cwd,
        env=env,
        check=True,
        stdout=subprocess.DEVNULL if supress_stdout else None,
    )


def capture(command: list[str], *, cwd: Path) -> str:
    return subprocess.check_output(command, cwd=cwd, text=True).strip()


def is_clean_worktree() -> bool:
    root = repo_root()
    status = capture(["git", "status", "--porcelain"], cwd=root)
    return status == ""


@task
def clean(
    c,
    include_tmp: bool = False,
    dry_run: bool = False,
) -> None:
    """Removes some temporary files and directories."""
    root = repo_root()

    if include_tmp:
        tmp = root / "tmp"
        # rm all except .gitkeep files in tmp
        for path in tmp.rglob("*"):
            if path.name != ".gitkeep":
                if VERBOSE:
                    print(f"clean {path}")
                if not dry_run:
                    if path.is_dir():
                        shutil.rmtree(path)
                    else:
                        path.unlink()


def clean_directory(path: Path, *, dry_run: bool) -> None:
    if path.exists():
        if VERBOSE:
            print(f"clean {path}")
        if not dry_run:
            shutil.rmtree(path)
    if not dry_run:
        path.mkdir(parents=True, exist_ok=True)


@task()
def run_tests(c, dry_run: bool = False) -> bool:
    """Runs unit tests."""
    try:
        run(
            [sys.executable, "-m", "pytest", "unit_tests"],
            cwd=repo_root(),
            dry_run=dry_run,
            env=dev_env(),
        )
    except subprocess.CalledProcessError as e:
        raise SystemExit(e.returncode)
    return True


@task()
def run_demos(c, critical_demos_only=False, exclude_scripts=False,
              exclude_notebooks=False, dry_run=False, timeout=300, pattern=""):
    """Run demos in a copied tree; retain logs, executed notebooks and JSON results."""
    if timeout <= 0:
        raise SystemExit("Demo timeout must be positive")
    demos = list(CRITICAL_DEMOS) if critical_demos_only else get_all_demos()
    demos = [d for d in demos if pattern in d
             and not (exclude_scripts and d.endswith(".py"))
             and not (exclude_notebooks and d.endswith(".ipynb"))]
    if not demos:
        raise SystemExit("No demos selected")
    if dry_run:
        for demo in demos:
            print("Would run:", demo)
        return True
    root = repo_root()
    (root / "tmp").mkdir(exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix="demo-check-", dir=root / "tmp"))
    shutil.copytree(root / "demos", output / "demos",
                    ignore=shutil.ignore_patterns("__pycache__", ".ipynb_checkpoints"))
    results = []
    for demo in demos:
        path = output / demo
        log = path.with_suffix(path.suffix + ".log")
        command = [sys.executable, str(path)]
        if path.suffix == ".ipynb":
            command = [sys.executable, "-m", "jupyter", "nbconvert", "--execute",
                       "--to", "notebook", "--ExecutePreprocessor.timeout=" + str(timeout),
                       "--output", path.stem + ".executed.ipynb", str(path)]
        started = time.monotonic()
        with log.open("w") as stream:
            process = subprocess.Popen(command, cwd=path.parent, env=dev_env(),
                                       stdin=subprocess.DEVNULL, stdout=stream,
                                       stderr=subprocess.STDOUT, start_new_session=True)
            try:
                code = process.wait(timeout=timeout)
                status = "passed" if code == 0 else "failed"
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                status, code = "timeout", None
        results.append(dict(demo=demo, status=status, returncode=code,
                            seconds=round(time.monotonic()-started, 2), log=str(log)))
        (output / "results.json").write_text(json.dumps(results, indent=2))
        print(f"{len(results)}/{len(demos)} {status}: {demo}", flush=True)
    print("Report:", output / "results.json")
    if any(r["status"] != "passed" for r in results):
        raise SystemExit(1)
    return True


def get_all_demos() -> list[str]:
    """Return a list of all demo scripts and notebooks in the demos directory, relative to the root directory."""
    root = repo_root()
    demos_dir = root / "demos"
    if not demos_dir.exists():
        raise SystemExit(f"Demos directory not found: {demos_dir}")
    demo_files = []
    for ext in (".py", ".ipynb"):
        demo_files.extend(
            str(path.relative_to(root)) for path in demos_dir.rglob(f"*{ext}")
            if ".ipynb_checkpoints" not in path.parts and "__pycache__" not in path.parts
        )
    return sorted(demo_files)


@task
def bump_version(c, version: str, dry_run: bool = False) -> None:
    """Updates the version number in all relevant files."""
    if not re.fullmatch(r"\d{2}\.\d{2}\.\d+", version):
        raise SystemExit("Version must use YY.MM.patch")
    root = repo_root()
    for file, pattern, replacement in VERSION_FILES:
        if not (root / file).exists():
            raise SystemExit(f"Expected version file not found: {file}")
        replace_once(
            root / file,
            pattern,
            replacement.replace("{version}", version),
            dry_run=dry_run,
        )
    replace_once(root / "conda-recipe/meta.yaml", r"git_tag:\s*[^\n]+",
                 f"git_tag: v{version}", dry_run=dry_run)


def replace_once(path: Path, pattern: str, replacement: str, *, dry_run: bool) -> None:
    text = path.read_text(encoding="utf-8")
    new_text, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise SystemExit(f"Expected exactly one version match in {path}")
    if text == new_text:
        return
    if VERBOSE:
        print(f"update {path.relative_to(repo_root())}")
    if not dry_run:
        path.write_text(new_text, encoding="utf-8")


@task()
def commit_release(
    c,
    version: str,
    message: str = "",
    source_branch: str = SOURCE_BRANCH,
    dry_run: bool = False,
) -> None:
    """Commits to the current branch in order to prepare the release."""
    root = repo_root()
    env = dev_env()

    # check if current branch is source_branch
    current_branch = capture(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=root)
    if current_branch != source_branch:
        raise SystemExit(
            f"Current branch is '{current_branch}', but expected '{source_branch}'. Please checkout the correct branch first."
        )

    allowed = {file for file, _, _ in VERSION_FILES} | {"CHANGELOG.md"}
    staged = set(capture(["git", "diff", "--cached", "--name-only"], cwd=root).splitlines())
    if staged - allowed:
        raise SystemExit("Unrelated staged changes; commit them separately first")

    # Skip when there are no changes to commit.
    comittable_changes = capture(["git", "status", "--porcelain"], cwd=root)
    if not comittable_changes:
        print("skip commit: nothing to commit")
        return

    run(
        ["git", "add", *[file for file, _, _ in VERSION_FILES] + ["CHANGELOG.md"]],
        cwd=root,
        env=env,
        dry_run=dry_run,
    )
    if not dry_run and not capture(["git", "diff", "--cached", "--name-only"], cwd=root):
        print("skip commit: no release metadata changes staged")
        return

    run(
        ["git", "commit", "-m", message or f"Prepare OCELOT {version} release"],
        cwd=root,
        env=env,
        dry_run=dry_run,
    )


@task()
def build(
    c, build_dir: str = BUILD_DIR, allow_dirty: bool = False, dry_run: bool = False
) -> None:
    """Builds the package."""

    if not allow_dirty and not is_clean_worktree():
        raise SystemExit(
            "Release preparation should start from a clean worktree. "
            "Commit, stash, or remove local changes first, or pass --allow-dirty."
        )

    root = repo_root()
    env = dev_env()

    build_dir = Path(build_dir).expanduser()
    if not build_dir.is_absolute():
        build_dir = root / build_dir
    if not dry_run:
        build_dir.mkdir(parents=True, exist_ok=True)

    # Keep old artifacts on disk, but never mix them with this build's output.
    artifacts_dir = (build_dir / "build-DRY-RUN" if dry_run else
                     Path(tempfile.mkdtemp(prefix="build-", dir=build_dir)))
    run([sys.executable, "-m", "build", "--outdir", str(artifacts_dir)],
        cwd=root, dry_run=dry_run, env=env)
    artifacts = sorted(str(path) for pattern in ("*.whl", "*.tar.gz")
                       for path in artifacts_dir.glob(pattern))
    if not dry_run and not artifacts:
        raise SystemExit("No Python artifacts built")
    if artifacts:
        run([sys.executable, "-m", "twine", "check", *artifacts],
            cwd=root, dry_run=dry_run, env=env)
    print("Build artifacts:", artifacts_dir)


@task()
def merge(
    c,
    version: str,
    source_branch: str = SOURCE_BRANCH,
    target_branch: str = TARGET_BRANCH,
    dry_run: bool = False,
) -> None:
    """Merges the source branch into the target branch with a merge commit."""

    if not is_clean_worktree():
        raise SystemExit(
            "Release preparation should start from a clean worktree. "
            "Commit, stash, or remove local changes first."
        )

    root = repo_root()

    run(["git", "checkout", target_branch], cwd=root, dry_run=dry_run)
    run(
        [
            "git",
            "merge",
            "--no-ff",
            source_branch,
            "-m",
            f"Merge {source_branch} into {target_branch} for {version} release",
        ],
        cwd=root,
        dry_run=dry_run,
    )


@task()
def tag(c, version: str, dry_run: bool = False) -> None:
    """Adds a special version tag to the current commit."""
    root = repo_root()
    tag = f"v{version}"

    run(
        ["git", "tag", "-a", tag, "-m", f"OCELOT {version}"],
        cwd=root,
        dry_run=dry_run,
    )


@task()
def push(
    c,
    source_branch: str = SOURCE_BRANCH,
    target_branch: str = TARGET_BRANCH,
    dry_run: bool = False,
    version: str = "",
) -> None:
    """Push branches and, optionally, one specified release tag."""
    root = repo_root()

    run(
        ["git", "push", "origin", source_branch, target_branch]
        + ([f"refs/tags/v{version}"] if version else []),
        cwd=root,
        dry_run=dry_run,
    )


def check_release(version):
    """Require clean, tagged source and consistent Python/Conda metadata."""
    root = repo_root()
    if not is_clean_worktree():
        raise SystemExit("Commit changes before publishing")
    if capture(["git", "rev-parse", "HEAD"], cwd=root) != capture(
            ["git", "rev-parse", f"v{version}^{{}}"], cwd=root):
        raise SystemExit("HEAD must be the requested release tag")
    metadata = tomllib.loads((root / "pyproject.toml").read_text())
    recipe = (root / "conda-recipe/meta.yaml").read_text()
    init = (root / "ocelot/__init__.py").read_text()
    if (metadata["project"]["version"] != version
            or not re.search(r"__version__\s*=\s*['\"]" + re.escape(version) + r"['\"]", init)
            or f'version: "{version}"' not in recipe
            or not re.search(r"git_tag:\s*v" + re.escape(version) + r"\s*$", recipe, re.M)):
        raise SystemExit("Release versions and Conda source tag must agree")


@task()
def publish(c, version, upload_pypi=False, upload_conda=False,
            conda=False, repository="pypi", dry_run=False):
    """Build/check tagged release; upload only with explicit upload flags."""
    check_release(version)
    root = repo_root()
    # Each build gets its own directory, preventing stale artifacts being uploaded.
    output = root / "dist" / f"release-{version}"
    if dry_run:
        print("Would build/check", version, "in a fresh directory under", output)
        print("Uploads:", dict(pypi=upload_pypi, conda=upload_conda, repository=repository))
        return
    output.mkdir(parents=True, exist_ok=True)
    artifacts_dir = Path(tempfile.mkdtemp(prefix="artifacts-", dir=output))
    run([sys.executable, "-m", "build", "--outdir", str(artifacts_dir)], cwd=root, dry_run=False)
    artifacts = sorted(str(p) for p in artifacts_dir.iterdir() if p.suffix == ".whl" or p.name.endswith(".tar.gz"))
    if not artifacts:
        raise SystemExit("No Python artifacts built")
    run([sys.executable, "-m", "twine", "check", *artifacts], cwd=root, dry_run=False)
    packages = []
    if conda or upload_conda:
        command = ["conda", "build", "--no-anaconda-upload", "--croot",
                   str(artifacts_dir / "conda"), "-c", "conda-forge", "conda-recipe"]
        run(command, cwd=root, dry_run=False)
        packages = [str(p) for p in (artifacts_dir / "conda").rglob("*")
                    if p.name.endswith((".conda", ".tar.bz2")) and p.parent.name in ("noarch", "linux-64", "osx-64", "osx-arm64", "win-64")]
        if not packages:
            raise SystemExit("No Conda artifacts built")
    if upload_pypi:
        run([sys.executable, "-m", "twine", "upload", "--repository", repository, *artifacts], cwd=root, dry_run=False)
    if upload_conda:
        run(["anaconda", "upload", *packages], cwd=root, dry_run=False)
    print("Release artifacts:", artifacts_dir)


@task()
def validate_for_commit(
    c,
    dry_run: bool = False,
) -> None:
    """Runs all checks which ensure code quality (e.g. tests and demos)."""
    run_tests(c, dry_run=dry_run)
    run_demos(c, critical_demos_only=True, dry_run=dry_run)


@task()
def do_full_release(
    c,
    version: str,
    commit_message: str = "",
    source_branch: str = SOURCE_BRANCH,
    target_branch: str = TARGET_BRANCH,
    dist_dir: str = BUILD_DIR,
    dry_run: bool = False,
) -> None:
    """Prepare locally, including merge and tag; never push or upload."""
    print("Local release preparation only; update CHANGELOG.md first.")
    clean(c, dry_run=dry_run)
    if not run_tests(c, dry_run=dry_run):
        print("\nWARNING: Some tests failed.")
        raise SystemExit(1)
    if not run_demos(c, dry_run=dry_run):
        print("\nWARNING: Some demos failed.")
        raise SystemExit(1)
    bump_version(c, version=version, dry_run=dry_run)
    commit_release(
        c,
        version=version,
        message=commit_message,
        source_branch=source_branch,
        dry_run=dry_run,
    )
    build(c, build_dir=dist_dir, dry_run=dry_run)
    merge(
        c,
        version=version,
        source_branch=source_branch,
        target_branch=target_branch,
        dry_run=dry_run,
    )
    tag(c, version=version, dry_run=dry_run)
    print("Local release prepared. Push and publish explicitly; see RELEASE.md.")
