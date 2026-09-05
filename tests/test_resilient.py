"""Crash consistency and negative recovery tests, independent of QE availability."""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import pytest
from qekit.core.errors import ErrorDeUso
from qekit.modules import resilient as r

pytestmark = pytest.mark.skipif(os.name != 'posix', reason='POSIX recovery worker')

INPUT = '''&CONTROL
calculation='scf', prefix='Si', pseudo_dir='./pp'
/
&SYSTEM
ibrav=2, celldm(1)=10.26, nat=2, ntyp=1, ecutwfc=30
/
&ELECTRONS
conv_thr=1e-10
/
ATOMIC_SPECIES
Si 28.085 Si.UPF
ATOMIC_POSITIONS crystal
Si 0 0 0
Si .25 .25 .25
K_POINTS automatic
4 4 4 0 0 0
'''

FAKE = r'''import pathlib,sys,time
p=pathlib.Path('.')
s=(p/'pw.in').read_text()
if 'failtest' in s:
 print('Error in routine: fake permanent failure',flush=True);sys.exit(1)
if 'waittest' in s:
 print('waiting',flush=True)
 while not (p/'Si.EXIT').exists(): time.sleep(.01)
restart="restart_mode     = 'restart'" in s
save=p/'out/Si.save';save.mkdir(parents=True,exist_ok=True)
(save/'charge-density.dat').write_text('charge')
(p/'out/Si.wfc1').write_text('wave')
if restart:
 (save/'data-file-schema.xml').write_text('<espresso><output><convergence_info><scf_conv><convergence_achieved>true</convergence_achieved><n_scf_steps>5</n_scf_steps></scf_conv></convergence_info></output></espresso>')
 print('convergence has been achieved\n! total energy = -22.0 Ry\nJOB DONE.',flush=True)
else:
 (save/'data-file-schema.xml').write_text('<espresso><output><convergence_info><scf_conv><convergence_achieved>false</convergence_achieved><n_scf_steps>2</n_scf_steps></scf_conv></convergence_info></output></espresso>')
 print('Maximum CPU time exceeded\nJOB DONE.',flush=True)
'''


@pytest.fixture
def job(tmp_path):
    (tmp_path / 'pp').mkdir()
    (tmp_path / 'pp/Si.UPF').write_text('pseudo')
    inp = tmp_path / 'scf.in'; inp.write_text(INPUT)
    fake = tmp_path / 'fake.py'; fake.write_text(FAKE)
    root = tmp_path / 'state'
    r.init(inp, root, [sys.executable, str(fake)], checkpoint_seconds=5, grace_seconds=.3)
    return root


def worker(root, *args):
    env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1]))
    return subprocess.Popen([sys.executable, '-m', 'qekit', 'resilient', 'run', str(root), *args],
                            env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def wait_running(root):
    for _ in range(200):
        state = json.loads((root/'state.json').read_text())
        if state['status']=='running':
            return state
        time.sleep(.02)
    raise AssertionError('worker did not start')


def test_checkpoint_then_restart_and_idempotent_success(job):
    assert r.run(job, max_segments=1)==75
    p,m,bad=r.latest(job)
    assert m['info']['status']=='checkpoint' and not bad
    assert r.run(job)==0
    count=json.loads((job/'state.json').read_text())['attempts']
    assert r.run(job)==0
    assert json.loads((job/'state.json').read_text())['attempts']==count
    assert r.status(job)['checkpoint_info']['energy_Ry']==-22.


def test_corrupt_latest_falls_back_without_using_partial(job):
    assert r.run(job,max_segments=1)==75
    older, m, _=r.latest(job)
    newer=r.publish(job,older/'work',r.load_job(job),m['info'])
    (newer/'work/out/Si.wfc1').write_text('corrupted')
    (job/'checkpoints/.partial-unfinished').mkdir()
    p,_,bad=r.latest(job)
    assert p==older and bad[0]['generation']==newer.name
    assert r.run(job)==0
    assert r.status(job)['state']['restarted_from']==older.name


def test_all_checkpoints_corrupt_refuses_fresh_start(job):
    r.run(job,max_segments=1)
    p,_,_=r.latest(job)
    (p/'work/out/Si.wfc1').unlink()
    with pytest.raises(ErrorDeUso,match='No intact checkpoint'):
        r.run(job)


def test_assets_and_binary_changes_rejected(job):
    (job/'assets/pp/Si.UPF').write_text('different physics')
    with pytest.raises(ErrorDeUso,match='fingerprint'):
        r.run(job)


def test_lock_keeps_other_workers_out(job):
    with r.locked(job):
        p=worker(job)
        _,err=p.communicate(timeout=10)
        assert p.returncode==76 and 'owns this job' in err


def test_signal_uses_exit_file_then_can_resume(job):
    # Change the frozen input before there are any checkpoints, for this test only.
    inp=job/'assets/original.in'; inp.write_text(inp.read_text().replace("calculation='scf'","title='waittest', calculation='scf'"))
    j=json.loads((job/'job.json').read_text());j['assets']=r.files_tree(job/'assets');r.atomic_json(job/'job.json',j)
    p=worker(job)
    wait_running(job)
    p.send_signal(signal.SIGTERM)
    out,err=p.communicate(timeout=10)
    assert p.returncode==75,(out,err)
    assert r.latest(job)[1]['info']['status']=='checkpoint'
    assert r.run(job)==0


def test_retry_budget_persists_and_bounds_permanent_failures(job,monkeypatch):
    inp=job/'assets/original.in';inp.write_text(inp.read_text().replace("calculation='scf'","title='failtest', calculation='scf'"))
    j=json.loads((job/'job.json').read_text());j['assets']=r.files_tree(job/'assets');j['max_failures']=2;r.atomic_json(job/'job.json',j)
    original=r.time.sleep
    monkeypatch.setattr(r.time,'sleep',lambda n:original(min(n,.01)))
    assert r.run(job)==2
    count=r.status(job)['state']['attempts']
    assert count==2 and r.run(job)==2 and r.status(job)['state']['attempts']==count


def test_partial_publication_does_not_replace_old_checkpoint(job,monkeypatch):
    r.run(job,max_segments=1)
    previous,manifest,_=r.latest(job)
    rename=r.os.rename
    def crash(source,target):
        if '.partial-' in str(source):
            raise OSError('simulated loss during publication')
        return rename(source,target)
    monkeypatch.setattr(r.os,'rename',crash)
    with pytest.raises(OSError):
        r.publish(job,previous/'work',r.load_job(job),manifest['info'])
    assert r.latest(job)[0]==previous


def test_no_job_done_success_without_convergence(job):
    r.run(job,max_segments=1)
    p,_,_=r.latest(job)
    (p/'work/pw.out').write_text('JOB DONE.\n')
    with pytest.raises(ErrorDeUso,match='Neither'):
        r.checkpoint_info(p/'work',r.load_job(job),0)


def test_service_stays_on_persistent_mount_and_stops_on_terminal_failure(job,tmp_path):
    text=r.service(job,tmp_path/'job.service','worker')
    assert 'RequiresMountsFor=' in text and 'KillMode=mixed' in text
    assert 'RestartPreventExitStatus=1 2 75' in text
    assert 'Restart=on-failure' in text


@pytest.mark.parametrize('change', ["calculation='bands'", "calculation='nscf'", "disk_io='none', calculation='scf'"])
def test_unsupported_restart_modes_rejected(tmp_path,change):
    inp=tmp_path/'in';inp.write_text(INPUT.replace("calculation='scf'",change))
    with pytest.raises(ErrorDeUso):
        r.init(inp,tmp_path/'state',sys.executable)


def test_convergence_coinciding_with_stop_request_is_success(job):
    assert r.run(job)==0
    p,_,_=r.latest(job)
    output=p/'work/pw.out'
    output.write_text(output.read_text()+'\nProgram stopped by user request\n')
    info=r.checkpoint_info(p/'work',r.load_job(job),0)
    assert info['status']=='succeeded' and info['energy_Ry']==-22.0


def test_changed_runtime_library_and_environment_refused(job, monkeypatch):
    monkeypatch.setenv('OMP_DYNAMIC', 'TRUE')
    with pytest.raises(ErrorDeUso, match='environment changed'):
        r.run(job)
    monkeypatch.delenv('OMP_DYNAMIC')
    j=r.load_job(job)
    library=next(iter(j['libraries']))
    j['libraries'][library]='wrong digest'
    r.atomic_json(job/'job.json',j)
    with pytest.raises(ErrorDeUso,match='binary changed'):
        r.run(job)


def test_changed_command_script_refused(job):
    j=r.load_job(job)
    script=next(Path(p) for p in j['binaries'] if p.endswith('fake.py'))
    script.write_text(script.read_text()+'\n# modified runtime\n')
    with pytest.raises(ErrorDeUso,match='binary changed'):
        r.run(job)


def test_supervisor_killed_child_retains_lock(job):
    inp=job/'assets/original.in'
    inp.write_text(inp.read_text().replace("calculation='scf'", "title='waittest', calculation='scf'"))
    j=json.loads((job/'job.json').read_text())
    j['assets']=r.files_tree(job/'assets')
    r.atomic_json(job/'job.json',j)
    parent=worker(job)
    child_pid=None
    try:
        for _ in range(400):
            state=json.loads((job/'state.json').read_text())
            if state.get('qe_pid'):
                child_pid=state['qe_pid']
                break
            time.sleep(.02)
        assert child_pid
        parent.kill()
        parent.wait(timeout=5)
        other=worker(job)
        _,error=other.communicate(timeout=10)
        assert other.returncode==76, error
    finally:
        if child_pid:
            try: os.killpg(child_pid, signal.SIGKILL)
            except ProcessLookupError: pass
        if parent.poll() is None:
            parent.kill()
        parent.communicate(timeout=5)
    for _ in range(100):
        try:
            with r.locked(job):
                break
        except r.BusyJob:
            time.sleep(.02)
    assert r.run(job)==0
    reports=[json.loads(p.read_text()) for p in (job/'attempts').glob('*/result.json')]
    assert any(p.get('status')=='interrupted' and p['qe_wall_seconds'] is None for p in reports)
    assert r.status(job)['state']['recovered_after_crash'] is True


def test_attempt_metrics_are_measured(job):
    assert r.run(job)==0
    for path in (job/'attempts').glob('*/result.json'):
        info=json.loads(path.read_text())
        assert info['wall_seconds'] >= info['qe_wall_seconds'] > 0
        assert info['restore_seconds'] >= 0
        assert info['checkpoint_seconds'] > 0
        assert info['checkpoint_bytes'] > 0


@pytest.mark.parametrize('corruption', [[], {'schema': 1}, 'success_flag'])
def test_corrupt_manifest_cannot_forge_success_or_displace_backup(job, corruption):
    r.run(job, max_segments=1)
    older, manifest, _ = r.latest(job)
    newer = r.publish(job, older/'work', r.load_job(job), manifest['info'])
    if corruption == 'success_flag':
        corruption = json.loads((newer/'manifest.json').read_text())
        corruption['info']['status'] = 'succeeded'
    r.atomic_json(newer/'manifest.json', corruption)
    assert r.latest(job)[0] == older
    assert r.run(job) == 0
    assert older.is_dir()
    assert newer.is_dir()  # Retained for diagnosis, never counted as a valid backup.


@pytest.mark.parametrize('options', [{'threads': 1.5}, {'keep': 2.5}, {'max_failures': True}])
def test_invalid_api_integer_options_are_rejected(tmp_path, options):
    with pytest.raises(ErrorDeUso, match='integer'):
        r.init(tmp_path/'absent.in', tmp_path/'state', **options)


def test_generated_service_is_accepted_by_systemd(job, tmp_path):
    import shutil
    if not shutil.which('systemd-analyze'):
        pytest.skip('systemd validator unavailable')
    path=tmp_path/'olla-test.service'
    r.service(job, path, 'worker')
    result=subprocess.run(['systemd-analyze','verify',str(path)],capture_output=True,text=True)
    assert result.returncode==0, result.stderr
    assert 'WorkingDirectory="' not in path.read_text()


def test_backwards_clock_does_not_select_an_old_checkpoint(job, monkeypatch):
    r.run(job, max_segments=1)
    previous, _, _ = r.latest(job)
    monkeypatch.setattr(r.time, 'time_ns', lambda: 1000)
    assert r.run(job) == 0
    latest, manifest, _ = r.latest(job)
    assert latest.name > previous.name
    assert manifest['info']['status'] == 'succeeded'


def test_missing_abandoned_attempt_does_not_lose_committed_checkpoint(job):
    r.run(job, max_segments=1)
    state=json.loads((job/'state.json').read_text())
    state.update(status='running', attempt='a'*32)
    r.atomic_json(job/'state.json',state)
    assert r.run(job)==0
    assert r.status(job)['state']['recovered_after_crash'] is True


@pytest.mark.parametrize("stop_kind", ["pause", "signal"])
def test_pause_during_restore_does_not_launch_qe(job, monkeypatch, stop_kind):
    r.run(job, max_segments=1)
    generation, _, _ = r.latest(job)
    attempts_before = r.status(job)['state']['attempts']
    restore = r.copy_durable
    def pausing_copy(source, target):
        restore(source, target)
        if Path(source) == generation/'work':
            if stop_kind == 'pause':
                r.pause(job)
            else:
                signal.raise_signal(signal.SIGTERM)
    monkeypatch.setattr(r, 'copy_durable', pausing_copy)
    def forbidden_launch(*args, **kwargs):
        raise AssertionError('QE was launched despite pending pause')
    monkeypatch.setattr(r.subprocess, 'Popen', forbidden_launch)
    assert r.run(job) == 75
    assert r.status(job)['state']['status'] == 'paused'
    assert r.latest(job)[0] == generation
    assert r.status(job)['state']['attempts'] == attempts_before


def test_pid_record_failure_reaps_qe_child(job, monkeypatch):
    inp=job/'assets/original.in'
    inp.write_text(inp.read_text().replace("calculation='scf'", "title='waittest', calculation='scf'"))
    metadata=json.loads((job/'job.json').read_text())
    metadata['assets']=r.files_tree(job/'assets')
    r.atomic_json(job/'job.json',metadata)
    save = r.atomic_json
    launch = r.subprocess.Popen
    children = []
    def capture(*args, **kwargs):
        child = launch(*args, **kwargs)
        children.append(child)
        return child
    def disk_full(path, data):
        if Path(path) == job/'state.json' and data.get('status') == 'running' and data.get('qe_pid'):
            raise OSError('simulated full disk while recording QE PID')
        return save(path, data)
    monkeypatch.setattr(r.subprocess, 'Popen', capture)
    monkeypatch.setattr(r, 'atomic_json', disk_full)
    try:
        with pytest.raises(OSError, match='full disk'):
            r.run(job)
        assert len(children) == 1
        assert children[0].poll() is not None, 'QE survived a failed supervisor startup'
        with r.locked(job):
            pass  # Reaping also released the child's inherited lock.
    finally:
        for child in children:
            if child.poll() is None:
                os.killpg(child.pid, signal.SIGKILL)
            child.wait(timeout=5)


def test_completed_status_does_not_expose_stale_pid(job):
    assert r.run(job) == 0
    assert 'qe_pid' not in r.status(job)['state']
    state=json.loads((job/'state.json').read_text())
    state['qe_pid']=999999
    r.atomic_json(job/'state.json',state)
    assert r.run(job)==0
    assert 'qe_pid' not in r.status(job)['state']
