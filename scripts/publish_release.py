#!/usr/bin/env python3
"""Build and publish OCELOT release artifacts.

Build and check artifacts without uploading:

    python scripts/publish_release.py 26.06.0

Upload explicitly after checking credentials:

    python scripts/publish_release.py 26.06.0 --upload-pypi --upload-conda

PyPI upload uses Twine credentials from the usual Twine/PyPI mechanisms, such
as `TWINE_USERNAME=__token__` and `TWINE_PASSWORD=<api-token>`.
Anaconda upload uses `anaconda-client`; run `anaconda login` first if needed.
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


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


def require_tool(name: str) -> None:
    if shutil.which(name) is None:
        raise SystemExit(f"Required command is not available: {name}")


def require_python_module(name: str) -> None:
    if importlib.util.find_spec(name) is None:
        raise SystemExit(f"Required Python module is not installed: {name}")


def require_clean_tracked_worktree(root: Path) -> None:
    status = capture(["git", "status", "--porcelain", "--untracked-files=no"], cwd=root)
    if status:
        raise SystemExit("Tracked files are modified. Commit or stash them before publishing.")


def require_release_tag(root: Path, version: str, *, allow_mismatch: bool) -> None:
    tag = f"v{version}"
    head = capture(["git", "rev-parse", "HEAD"], cwd=root)
    try:
        tagged = capture(["git", "rev-parse", f"{tag}^{{}}"], cwd=root)
    except subprocess.CalledProcessError as exc:
        raise SystemExit(f"Release tag not found: {tag}") from exc
    if head != tagged and not allow_mismatch:
        raise SystemExit(
            f"HEAD is not {tag}. Check out the release tag or pass --allow-tag-mismatch."
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


def clean_directory(path: Path, *, dry_run: bool) -> None:
    if path.exists():
        print(f"clean {path}")
        if not dry_run:
            shutil.rmtree(path)
    if not dry_run:
        path.mkdir(parents=True, exist_ok=True)


def build_pypi(root: Path, args: argparse.Namespace, env: dict[str, str]) -> list[Path]:
    dist_dir = Path(args.dist_dir).expanduser()
    if not dist_dir.is_absolute():
        dist_dir = root / dist_dir

    if not args.keep_dist:
        clean_directory(dist_dir, dry_run=args.dry_run)
    elif not args.dry_run:
        dist_dir.mkdir(parents=True, exist_ok=True)

    require_python_module("build")
    run([sys.executable, "-m", "build", "--outdir", str(dist_dir)], cwd=root, dry_run=args.dry_run, env=env)

    artifacts = sorted(dist_dir.glob("*"))
    if not args.dry_run and not artifacts:
        raise SystemExit(f"No PyPI artifacts found in {dist_dir}")

    require_python_module("twine")
    run([sys.executable, "-m", "twine", "check", *map(str, artifacts)], cwd=root, dry_run=args.dry_run, env=env)
    return artifacts


def upload_pypi(root: Path, args: argparse.Namespace, artifacts: list[Path], env: dict[str, str]) -> None:
    if not args.upload_pypi:
        print("skip PyPI upload")
        return
    if not artifacts and not args.dry_run:
        raise SystemExit("No PyPI artifacts available to upload")

    command = [sys.executable, "-m", "twine", "upload"]
    if args.pypi_repository:
        command.extend(["--repository", args.pypi_repository])
    if args.pypi_repository_url:
        command.extend(["--repository-url", args.pypi_repository_url])
    command.extend(map(str, artifacts))
    run(command, cwd=root, dry_run=args.dry_run, env=env)


def build_conda(root: Path, args: argparse.Namespace, env: dict[str, str]) -> Path | None:
    require_tool("conda")

    recipe = Path(args.conda_recipe).expanduser()
    if not recipe.is_absolute():
        recipe = root / recipe

    croot = Path(args.conda_croot).expanduser()
    if not croot.is_absolute():
        croot = root / croot

    command = ["conda", "build"]
    command.extend(["--croot", str(croot)])
    for channel in args.conda_build_channel:
        command.extend(["-c", channel])
    command.append(str(recipe))
    run(command, cwd=root, dry_run=args.dry_run, env=env)

    output_command = ["conda", "build", "--output"]
    output_command.extend(["--croot", str(croot)])
    for channel in args.conda_build_channel:
        output_command.extend(["-c", channel])
    output_command.append(str(recipe))

    if args.dry_run:
        print("+ " + " ".join(output_command))
        return None

    output = capture(output_command, cwd=root)
    packages = [Path(line) for line in output.splitlines() if line.strip()]
    if not packages:
        raise SystemExit("conda-build did not report an output package")
    return packages[-1]


def upload_conda(root: Path, args: argparse.Namespace, package: Path | None, env: dict[str, str]) -> None:
    if not args.upload_conda:
        print("skip Anaconda upload")
        return
    require_tool("anaconda")

    command = ["anaconda", "upload"]
    if args.conda_upload_user:
        command.extend(["--user", args.conda_upload_user])
    for label in args.conda_label:
        command.extend(["--label", label])
    if package is not None:
        command.append(str(package))
    run(command, cwd=root, dry_run=args.dry_run, env=env)


def parse_args() -> argparse.Namespace:
    default_dist = Path(tempfile.gettempdir()) / "ocelot-release-dist"
    default_conda_croot = Path(tempfile.gettempdir()) / "ocelot-conda-bld"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("version", help="Release version, for example 26.06.0")
    parser.add_argument("--dist-dir", default=str(default_dist))
    parser.add_argument("--conda-recipe", default="conda-recipe")
    parser.add_argument("--conda-croot", default=str(default_conda_croot))
    parser.add_argument("--conda-build-channel", action="append", default=["conda-forge"])
    parser.add_argument("--conda-upload-user", default="ocelot-collab")
    parser.add_argument("--conda-label", action="append", default=["main"])
    parser.add_argument("--pypi-repository")
    parser.add_argument("--pypi-repository-url")
    parser.add_argument("--skip-pypi-build", action="store_true")
    parser.add_argument("--skip-conda-build", action="store_true")
    parser.add_argument("--upload-pypi", action="store_true")
    parser.add_argument("--upload-conda", action="store_true")
    parser.add_argument("--allow-tag-mismatch", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true")
    parser.add_argument("--keep-dist", action="store_true", help="Do not clean the PyPI dist directory before building")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = repo_root()
    env = release_env(root)

    if not args.allow_dirty:
        require_clean_tracked_worktree(root)
    require_release_tag(root, args.version, allow_mismatch=args.allow_tag_mismatch)

    pypi_artifacts: list[Path] = []
    conda_package: Path | None = None

    if args.skip_pypi_build:
        print("skip PyPI build")
    else:
        pypi_artifacts = build_pypi(root, args, env)
    upload_pypi(root, args, pypi_artifacts, env)

    if args.skip_conda_build:
        print("skip conda build")
    else:
        conda_package = build_conda(root, args, env)
    upload_conda(root, args, conda_package, env)

    print("\nPublishing workflow complete.")


if __name__ == "__main__":
    main()
