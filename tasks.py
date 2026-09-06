#!/usr/bin/env python3
"""
Common tasks for developers only.
"""

from __future__ import annotations

from invoke import task
import importlib.util
import os
import re
import subprocess
import sys
import shutil
import traceback
import tempfile
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

BUILD_DIR = Path("build")


def repo_root() -> Path:
    return Path(__file__).resolve().parent


def dev_env() -> dict[str, str]:
    env = os.environ.copy()
    # Set MPLBACKEND to Agg to prevent matplotlib from trying to open a window and blocking the execution of the script until manual intervention.
    env.setdefault("MPLBACKEND", "Agg")
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


def promt_continue() -> None:
    """Prompt the user to continue, abort if the answer is not yes."""
    input_str = input("Continue [Y/n]? ").strip().lower()
    if input_str not in ("y", "yes", ""):
        raise SystemExit("Abort.")


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
        return False
    return True


@task()
def run_demos(c, critical_demos_only: bool = False, dry_run: bool = False) -> bool:
    """Runs demo scripts and notebooks."""

    if critical_demos_only:
        demos = CRITICAL_DEMOS
    else:
        demos = get_all_demos()

    root = repo_root()
    env = dev_env()
    cwd = root / "tmp" / "demos"
    cwd.mkdir(parents=True, exist_ok=True)

    counters = {
        "script_success": 0,
        "script_total": 0,
        "notebook_success": 0,
        "notebook_total": 0,
    }
    for demo in demos:

        demo_path = root / demo
        if not demo_path.exists():
            raise SystemExit(f"Demo not found: {demo}")

        if demo_path.suffix == ".py":
            kind = "script"
            cmd = [sys.executable, str(demo_path)]

        if demo_path.suffix == ".ipynb":
            kind = "notebook"
            cmd = [
                sys.executable,
                "-m",
                "jupyter",
                "nbconvert",
                "--to",
                "notebook",
                "--execute",
                "--output-dir",
                str(cwd),
                str(demo_path),
            ]
        try:
            counters[kind + "_total"] += 1
            run(cmd, cwd=cwd, dry_run=dry_run, env=env, supress_stdout=True)
            counters[kind + "_success"] += 1
        except subprocess.CalledProcessError as e:
            print(f"FAILED: Error encountered while executing {demo}:")
            traceback.print_exception(type(e), e, e.__traceback__)

    print("Demo execution summary:")
    print(
        f"scripts executed successfully: {counters['script_success']}/{counters['script_total']}"
    )
    print(
        f"notebooks executed successfully: {counters['notebook_success']}/{counters['notebook_total']}"
    )

    return (counters["script_success"] == counters["script_total"]) and (
        counters["notebook_success"] == counters["notebook_total"]
    )


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
        )
    return sorted(demo_files)


@task
def bump_version(c, version: str, dry_run: bool = False) -> None:
    """Updates the version number in all relevant files."""
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
    message: str,
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

    # check if there are any changes (staged or not) to commit, if not, skip the commit step
    comittable_changes = capture(["git", "status", "--porcelain"], cwd=root)
    if not comittable_changes:
        print("skip commit: nothing to commit")
        promt_continue()
        return

    run(
        ["git", "add", *[file for file, _, _ in VERSION_FILES] + ["CHANGELOG.md"]],
        cwd=root,
        env=env,
        dry_run=dry_run,
    )
    run(
        ["git", "commit", "-m", message or f"Prepare OCELOT {version} release"],
        cwd=root,
        env=env,
        dry_run=True,
    )

    staged = run(
        ["git", "diff", "--cached", "--quiet"], cwd=root, env=env, dry_run=dry_run
    )
    if staged.returncode == 0:
        print("skip commit: no release metadata changes staged")
        return

    run(
        ["git", "commit", "-m", message or f"Prepare OCELOT {version} release"],
        cwd=root,
        env=env,
        dry_run=dry_run,
    )


@task()
def build(c, build_dir: str=BUILD_DIR, allow_dirty: bool = False, dry_run: bool = False) -> None:
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
    build_dir.mkdir(parents=True, exist_ok=True)

    if importlib.util.find_spec("build") is not None:
        run(
            [sys.executable, "-m", "build", "--outdir", str(build_dir)],
            cwd=root,
            dry_run=dry_run,
            env=env,
        )
    else:
        run(
            [
                sys.executable,
                "setup.py",
                "sdist",
                "--dist-dir",
                str(build_dir),
                "bdist_wheel",
                "--dist-dir",
                str(build_dir),
            ],
            cwd=root,
            dry_run=dry_run,
            env=env,
        )

    artifacts = sorted(str(path) for path in build_dir.glob("*"))
    if artifacts and importlib.util.find_spec("twine") is not None:
        run(
            [sys.executable, "-m", "twine", "check", *artifacts],
            cwd=root,
            dry_run=dry_run,
            env=env,
        )
    elif artifacts:
        print("skip twine check: twine is not installed")


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
) -> None:
    """Pushes the source and target branches and tags to the remote repository."""
    root = repo_root()

    run(
        ["git", "push", "origin", source_branch, target_branch, "--tags"],
        cwd=root,
        dry_run=dry_run,
    )


@task()
def publish(c, dry_run: bool = False) -> None:
    """Publishes the package (not implemented yet)."""
    if dry_run:
        print("Dry run: skipping publishing to PyPI & conda.")
        return
    raise NotImplementedError(
        "Publishing is not implemented yet. Please continue manually or implement it here."
    )


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
    commit_message: str,
    source_branch: str = SOURCE_BRANCH,
    target_branch: str = TARGET_BRANCH,
    dist_dir: str = BUILD_DIR,
    dry_run: bool = False,
) -> None:
    """Execute all tasks required for a full release."""
    print("Did you remember to update the CHANGELOG.md file?")
    promt_continue()
    clean(c, dry_run=dry_run)
    if not run_tests(c, dry_run=dry_run):
        print("\nWARNING: Some tests failed.")
        promt_continue()
    if not run_demos(c, dry_run=dry_run):
        print("\nWARNING: Some demos failed.")
        promt_continue()
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
    push(c, dry_run=dry_run)
    publish(c, dry_run=dry_run)
