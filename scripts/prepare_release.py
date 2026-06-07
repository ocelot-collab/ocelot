#!/usr/bin/env python3
"""Prepare an OCELOT release from the dev branch.

Typical release-prep run:

    python scripts/prepare_release.py 26.06.0 --full

The script updates version metadata, runs tests and smoke demos, builds package
artifacts, and can commit, merge dev into master, and create a git tag. It does
not publish to PyPI, Anaconda, or GitHub Releases unless the explicit optional
flags are used.
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path


DEFAULT_DEMOS = (
    "demos/ebeam/rk_track.py",
    "demos/ebeam/dba.py",
    "demos/ebeam/dba_tracking.py",
)

VERSION_FILES = (
    "setup.py",
    "ocelot/__init__.py",
    "conda-recipe/meta.yaml",
)


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def run(
    command: list[str],
    *,
    cwd: Path,
    dry_run: bool,
    env: dict[str, str] | None = None,
) -> None:
    print("+ " + " ".join(command))
    if dry_run:
        return
    subprocess.run(command, cwd=cwd, env=env, check=True)


def capture(command: list[str], *, cwd: Path) -> str:
    return subprocess.check_output(command, cwd=cwd, text=True).strip()


def require_clean_worktree(root: Path) -> None:
    status = capture(["git", "status", "--porcelain"], cwd=root)
    if status:
        raise SystemExit(
            "Release preparation should start from a clean worktree. "
            "Commit, stash, or remove local changes first, or pass --allow-dirty."
        )


def replace_once(path: Path, pattern: str, replacement: str, *, dry_run: bool) -> None:
    text = path.read_text(encoding="utf-8")
    new_text, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise SystemExit(f"Expected exactly one version match in {path}")
    if text == new_text:
        return
    print(f"update {path.relative_to(repo_root())}")
    if not dry_run:
        path.write_text(new_text, encoding="utf-8")


def update_versions(root: Path, version: str, *, dry_run: bool) -> None:
    replace_once(
        root / "setup.py",
        r"version=['\"][^'\"]+['\"]",
        f"version='{version}'",
        dry_run=dry_run,
    )
    replace_once(
        root / "ocelot/__init__.py",
        r"__version__\s*=\s*['\"][^'\"]+['\"]",
        f"__version__ = '{version}'",
        dry_run=dry_run,
    )
    replace_once(
        root / "conda-recipe/meta.yaml",
        r'version:\s*"[^"]+"',
        f'version: "{version}"',
        dry_run=dry_run,
    )
    replace_once(
        root / "conda-recipe/meta.yaml",
        r"git_tag:\s*v[^\s]+",
        f"git_tag: v{version}",
        dry_run=dry_run,
    )


def release_env(root: Path) -> dict[str, str]:
    env = os.environ.copy()
    env.setdefault("MPLBACKEND", "Agg")
    env["PYTHONPATH"] = (
        str(root)
        if not env.get("PYTHONPATH")
        else str(root) + os.pathsep + env["PYTHONPATH"]
    )
    return env


def run_tests(root: Path, args: argparse.Namespace, env: dict[str, str]) -> None:
    if args.skip_tests:
        print("skip tests")
        return
    run([sys.executable, "-m", "pytest", args.tests], cwd=root, dry_run=args.dry_run, env=env)


def run_demos(root: Path, args: argparse.Namespace, env: dict[str, str]) -> None:
    if args.skip_demos:
        print("skip demos")
        return
    demos = args.demo or list(DEFAULT_DEMOS)
    for demo in demos:
        demo_path = root / demo
        if not demo_path.exists():
            raise SystemExit(f"Demo not found: {demo}")
        run([sys.executable, demo], cwd=root, dry_run=args.dry_run, env=env)


def build_package(root: Path, args: argparse.Namespace, env: dict[str, str]) -> None:
    if args.skip_build:
        print("skip package build")
        return

    dist_dir = Path(args.dist_dir).expanduser()
    if not dist_dir.is_absolute():
        dist_dir = root / dist_dir
    dist_dir.mkdir(parents=True, exist_ok=True)

    if importlib.util.find_spec("build") is not None:
        run(
            [sys.executable, "-m", "build", "--outdir", str(dist_dir)],
            cwd=root,
            dry_run=args.dry_run,
            env=env,
        )
    else:
        run(
            [
                sys.executable,
                "setup.py",
                "sdist",
                "--dist-dir",
                str(dist_dir),
                "bdist_wheel",
                "--dist-dir",
                str(dist_dir),
            ],
            cwd=root,
            dry_run=args.dry_run,
            env=env,
        )

    artifacts = sorted(str(path) for path in dist_dir.glob("*"))
    if artifacts and importlib.util.find_spec("twine") is not None:
        run([sys.executable, "-m", "twine", "check", *artifacts], cwd=root, dry_run=args.dry_run, env=env)
    elif artifacts:
        print("skip twine check: twine is not installed")


def commit_release(root: Path, args: argparse.Namespace) -> None:
    if not args.commit:
        return
    run(["git", "add", *VERSION_FILES, "CHANGELOG.md"], cwd=root, dry_run=args.dry_run)
    if args.dry_run:
        run(
            ["git", "commit", "-m", args.commit_message or f"Prepare OCELOT {args.version} release"],
            cwd=root,
            dry_run=True,
        )
        return

    staged = subprocess.run(["git", "diff", "--cached", "--quiet"], cwd=root)
    if staged.returncode == 0:
        print("skip commit: no release metadata changes staged")
        return

    run(
        ["git", "commit", "-m", args.commit_message or f"Prepare OCELOT {args.version} release"],
        cwd=root,
        dry_run=args.dry_run,
    )


def merge_and_tag(root: Path, args: argparse.Namespace) -> None:
    tag = f"{args.tag_prefix}{args.version}"

    if args.merge_master:
        run(["git", "checkout", args.target_branch], cwd=root, dry_run=args.dry_run)
        run(
            [
                "git",
                "merge",
                "--no-ff",
                args.source_branch,
                "-m",
                f"Merge {args.source_branch} into {args.target_branch} for {args.version} release",
            ],
            cwd=root,
            dry_run=args.dry_run,
        )

    if args.tag:
        run(["git", "tag", "-a", tag, "-m", f"OCELOT {args.version}"], cwd=root, dry_run=args.dry_run)

    if args.github_release:
        run(
            [
                "gh",
                "release",
                "create",
                tag,
                "--target",
                args.target_branch,
                "--title",
                f"OCELOT {args.version}",
                "--generate-notes",
            ],
            cwd=root,
            dry_run=args.dry_run,
        )

    if args.push:
        run(
            ["git", "push", "origin", args.source_branch, args.target_branch, "--tags"],
            cwd=root,
            dry_run=args.dry_run,
        )


def parse_args() -> argparse.Namespace:
    default_dist = Path(tempfile.gettempdir()) / "ocelot-release-dist"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("version", help="Release version, for example 26.06.0")
    parser.add_argument("--source-branch", default="dev")
    parser.add_argument("--target-branch", default="master")
    parser.add_argument("--tests", default="unit_tests", help="pytest target")
    parser.add_argument("--demo", action="append", help="Demo script to run; can be passed more than once")
    parser.add_argument("--dist-dir", default=str(default_dist), help="Package artifact output directory")
    parser.add_argument("--skip-tests", action="store_true")
    parser.add_argument("--skip-demos", action="store_true")
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--commit", action="store_true")
    parser.add_argument("--commit-message")
    parser.add_argument("--merge-master", action="store_true")
    parser.add_argument("--tag", action="store_true")
    parser.add_argument("--tag-prefix", default="v")
    parser.add_argument("--github-release", action="store_true")
    parser.add_argument("--push", action="store_true")
    parser.add_argument(
        "--full",
        action="store_true",
        help="Equivalent to --commit --merge-master --tag, but still does not push or publish",
    )
    parser.add_argument("--allow-dirty", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.full:
        args.commit = True
        args.merge_master = True
        args.tag = True

    if (args.merge_master or args.push) and not args.commit:
        raise SystemExit("--merge-master and --push require --commit or --full")

    root = repo_root()
    branch = capture(["git", "branch", "--show-current"], cwd=root)
    if branch != args.source_branch and (args.commit or args.merge_master):
        raise SystemExit(f"Run release preparation from {args.source_branch}; current branch is {branch}")
    if not args.allow_dirty:
        require_clean_worktree(root)

    env = release_env(root)
    update_versions(root, args.version, dry_run=args.dry_run)
    run_tests(root, args, env)
    run_demos(root, args, env)
    build_package(root, args, env)
    commit_release(root, args)
    merge_and_tag(root, args)

    tag = f"{args.tag_prefix}{args.version}"
    print("\nRelease preparation complete.")
    print(f"Tag: {tag}")
    print("Publish steps remain explicit:")
    print("  python -m twine upload <dist-dir>/*")
    print("  conda build <recipe> && anaconda upload <package>")


if __name__ == "__main__":
    main()
