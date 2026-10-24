"""Release checks must fail before any upload or unrelated Git mutation."""
import importlib.util
import json
import subprocess
from pathlib import Path

import pytest

pytest.importorskip('invoke')
from invoke import Context

spec = importlib.util.spec_from_file_location('release_tasks', Path(__file__).resolve().parents[2] / 'tasks.py')
tasks = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tasks)


@pytest.fixture
def repository(tmp_path, monkeypatch):
    monkeypatch.setattr(tasks, 'repo_root', lambda: tmp_path)
    (tmp_path/'ocelot').mkdir()
    (tmp_path/'conda-recipe').mkdir()
    (tmp_path/'pyproject.toml').write_text('[project]\nversion = "26.06.1"\n')
    (tmp_path/'ocelot/__init__.py').write_text('__version__ = "26.06.1"\n')
    (tmp_path/'conda-recipe/meta.yaml').write_text('package:\n  version: "26.06.1"\nsource:\n  git_tag: v26.06.1\n')
    (tmp_path/'CHANGELOG.md').write_text('# Changes\n')
    def git(*args):
        return subprocess.check_output(['git', *args], cwd=tmp_path, text=True).strip()
    git('init', '-b', 'dev')
    git('config', 'user.name', 'Test')
    git('config', 'user.email', 'test@example.invalid')
    git('add', '.')
    git('commit', '-m', 'initial')
    return tmp_path, git


def test_version_commit_and_tag_validation(repository):
    root, git = repository
    tasks.bump_version(Context(), '26.10.0')
    assert 'git_tag: v26.10.0' in (root/'conda-recipe/meta.yaml').read_text()
    tasks.commit_release(Context(), '26.10.0')
    assert git('status', '--porcelain') == ''
    assert git('log', '-1', '--format=%s') == 'Prepare OCELOT 26.10.0 release'
    git('tag', 'v26.10.0')
    tasks.check_release('26.10.0')
    (root/'CHANGELOG.md').write_text('dirty')
    with pytest.raises(SystemExit, match='Commit changes'):
        tasks.check_release('26.10.0')


def test_refuse_unrelated_staged_file(repository):
    root, git = repository
    (root/'unrelated').write_text('work')
    git('add', 'unrelated')
    with pytest.raises(SystemExit, match='Unrelated'):
        tasks.commit_release(Context(), '26.10.0')


def test_bump_dry_run_leaves_source_unchanged(repository):
    root, git = repository
    tasks.bump_version(Context(), '26.10.0', dry_run=True)
    assert git('status', '--porcelain') == ''


@pytest.mark.parametrize('upload', [False, True])
def test_publishing_is_explicit_and_checks_before_upload(repository, monkeypatch, upload):
    root, git = repository
    git('tag', 'v26.06.1')
    calls = []
    def fake_run(command, **kwargs):
        calls.append(command)
        if 'build' in command:
            out = Path(command[command.index('--outdir')+1])
            (out/'ocelot_collab-26.06.1-py3-none-any.whl').write_text('wheel')
            (out/'ocelot_collab-26.06.1.tar.gz').write_text('sdist')
    monkeypatch.setattr(tasks, 'run', fake_run)
    tasks.publish(Context(), '26.06.1', upload_pypi=upload)
    assert any('check' in cmd for cmd in calls)
    assert any('upload' in cmd for cmd in calls) == upload
    if upload:
        assert calls[-1][3:5] == ['upload', '--repository']


def test_failed_build_never_uploads(repository, monkeypatch):
    root, git = repository
    git('tag', 'v26.06.1')
    calls = []
    def fail(command, **kwargs):
        calls.append(command)
        raise subprocess.CalledProcessError(1, command)
    monkeypatch.setattr(tasks, 'run', fail)
    with pytest.raises(subprocess.CalledProcessError):
        tasks.publish(Context(), '26.06.1', upload_pypi=True, upload_conda=True)
    assert not any('upload' in cmd for cmd in calls)


def test_unit_test_failure_is_nonzero(monkeypatch):
    def fail(command, **kwargs):
        raise subprocess.CalledProcessError(1, command)
    monkeypatch.setattr(tasks, 'run', fail)
    with pytest.raises(SystemExit):
        tasks.run_tests(Context())


def test_demo_copy_preserves_inputs_and_reports_failures(repository):
    root, git = repository
    demos = root/'demos'
    demos.mkdir()
    (demos/'good.py').write_text('from pathlib import Path\nPath("generated.txt").write_text("done")\n')
    (demos/'bad.py').write_text('raise RuntimeError("demo failure")\n')
    with pytest.raises(SystemExit):
        tasks.run_demos(Context(), timeout=10)
    assert not (demos/'generated.txt').exists()
    report = next((root/'tmp').glob('demo-check-*/results.json'))
    results = json.loads(report.read_text())
    assert [r['status'] for r in results] == ['failed', 'passed']
    assert 'demo failure' in Path(results[0]['log']).read_text()


def test_conda_build_and_upload_are_separate(repository, monkeypatch):
    root, git = repository
    git('tag', 'v26.06.1')
    calls = []
    def fake_run(command, **kwargs):
        calls.append(command)
        if command[0] == 'conda':
            out = Path(command[command.index('--croot')+1])/'noarch'
            out.mkdir(parents=True)
            (out/'ocelot-26.06.1-0.conda').write_text('conda')
        elif 'build' in command:
            out = Path(command[command.index('--outdir')+1])
            (out/'ocelot_collab-26.06.1-py3-none-any.whl').write_text('wheel')
    monkeypatch.setattr(tasks, 'run', fake_run)
    tasks.publish(Context(), '26.06.1', upload_conda=True)
    assert '--no-anaconda-upload' in calls[-2]
    assert calls[-1][0:2] == ['anaconda', 'upload']
    assert not any('twine' in c and 'upload' in c for c in calls)


def test_publish_dry_run_does_not_build(repository, monkeypatch):
    root, git = repository
    git('tag', 'v26.06.1')
    monkeypatch.setattr(tasks, 'run', lambda *a, **kw: pytest.fail('must not execute'))
    tasks.publish(Context(), '26.06.1', upload_pypi=True, dry_run=True)
    assert not (root/'dist').exists()


def test_local_release_never_pushes_or_uploads(monkeypatch):
    calls = []
    for name in ['clean', 'run_tests', 'run_demos', 'bump_version', 'commit_release',
                 'build', 'merge', 'tag', 'push', 'publish']:
        def record(*args, _name=name, **kwargs):
            calls.append(_name)
            return True
        monkeypatch.setattr(tasks, name, record)
    tasks.do_full_release(Context(), '26.10.0')
    assert 'push' not in calls and 'publish' not in calls
    assert calls[-1] == 'tag'


def test_push_only_requested_tag(monkeypatch):
    commands = []
    monkeypatch.setattr(tasks, 'run', lambda command, **kw: commands.append(command))
    tasks.push(Context(), version='26.10.0')
    assert commands == [['git', 'push', 'origin', 'dev', 'master', 'refs/tags/v26.10.0']]


def test_demo_timeout_is_reported(repository):
    root, git = repository
    (root/'demos').mkdir()
    (root/'demos/slow.py').write_text('import time\ntime.sleep(30)\n')
    with pytest.raises(SystemExit):
        tasks.run_demos(Context(), timeout=1)
    report = next((root/'tmp').glob('demo-check-*/results.json'))
    assert json.loads(report.read_text())[0]['status'] == 'timeout'


def test_build_ignores_stale_artifacts(repository, monkeypatch):
    root, git = repository
    (root/'dist').mkdir()
    stale = root/'dist/old.whl'
    stale.write_text('old build')
    commands = []
    def fake_run(command, **kwargs):
        commands.append(command)
        if 'build' in command:
            output = Path(command[command.index('--outdir')+1])
            (output/'new.whl').write_text('new build')
    monkeypatch.setattr(tasks, 'run', fake_run)
    tasks.build(Context(), allow_dirty=True)
    assert stale.read_text() == 'old build'
    checked = commands[-1]
    assert 'check' in checked and str(stale) not in checked
    assert checked[-1].endswith('/new.whl')
