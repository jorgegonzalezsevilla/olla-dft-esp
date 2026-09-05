"""Durable, single-worker pw.x recovery for preemptible POSIX machines.

Only cleanly stopped QE generations are committed. Each attempt writes to a
private workspace; a hard kill cannot overwrite the last committed generation.
The state directory must live on a retained, single-writer persistent disk.
"""
import contextlib
import hashlib
import io
import json
import math
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import signal
import subprocess
import sys
import time
import uuid
import xml.etree.ElementTree as ET

from qekit.core.errors import ErrorDeUso

class BusyJob(ErrorDeUso):
    """Temporary ownership conflict; service may retry after the child exits."""


SCHEMA = 1
SUPPORTED = {'scf', 'relax', 'vc-relax'}
THREAD_KEYS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def _sync_dir(path):
    fd = os.open(str(path), os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def atomic_json(path, data):
    path = Path(path)
    tmp = path.with_name('.' + path.name + '.' + uuid.uuid4().hex + '.tmp')
    try:
        with tmp.open('x', encoding='utf-8') as f:
            json.dump(data, f, indent=2, sort_keys=True, allow_nan=False)
            f.write('\n')
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
        _sync_dir(path.parent)
    finally:
        tmp.unlink(missing_ok=True)


def _load(path):
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError) as exc:
        raise ErrorDeUso(f'Cannot read durable state {path}: {exc}') from exc


@contextlib.contextmanager
def locked(root):
    if os.name != 'posix':
        raise ErrorDeUso('Resilient QE currently requires POSIX/Linux and a persistent disk.')
    import fcntl
    with (Path(root) / 'worker.lock').open('a+') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise BusyJob('Another worker or its QE child still owns this job.') from None
        # QE inherits the lock: killing only the supervisor cannot create two workers.
        yield lock.fileno()


def files_tree(root):
    if Path(root).is_symlink() or not Path(root).is_dir():
        raise ErrorDeUso('Expected a regular checkpoint directory.')
    result = {}
    for p in sorted(Path(root).rglob('*')):
        if p.is_symlink():
            raise ErrorDeUso(f'Symlinks are not accepted in recoverable state: {p}')
        if p.is_file():
            result[p.relative_to(root).as_posix()] = digest(p)
        elif not p.is_dir():
            raise ErrorDeUso(f'Not a regular checkpoint file: {p}')
    return result


def copy_durable(source, target):
    """Copy regular files and fsync before publication; never hardlink mutable QE files."""
    source, target = Path(source), Path(target)
    target.mkdir()
    for p in sorted(source.iterdir()):
        dest = target / p.name
        if p.is_symlink():
            raise ErrorDeUso(f'Symlink in checkpoint: {p}')
        if p.is_dir():
            copy_durable(p, dest)
        elif p.is_file():
            shutil.copyfile(p, dest)
            with dest.open('rb') as f:
                os.fsync(f.fileno())
        else:
            raise ErrorDeUso(f'Unsupported checkpoint entry: {p}')
    _sync_dir(target)


def _namelists(text):
    from ase.io.espresso import read_fortran_namelist
    from ase.io.espresso_namelist.namelist import Namelist
    params, cards = read_fortran_namelist(io.StringIO(text))
    return Namelist(params), cards


def _positive(value, name):
    if not math.isfinite(float(value)) or float(value) <= 0:
        raise ErrorDeUso(f'{name} must be positive and finite.')


def runtime_environment():
    # Thread counts and the temporary directory are set by the worker itself.
    return {k: v for k, v in os.environ.items()
            if k.startswith(('LD_', 'OMP_', 'MKL_', 'OPENBLAS_', 'OMPI_', 'I_MPI_'))
            and k not in THREAD_KEYS}


def runtime_libraries(binaries):
    libraries = {}
    for name in binaries:
        with Path(name).open('rb') as stream:
            if stream.read(4) != b'\x7fELF':
                continue
        result = subprocess.run(['ldd', name], capture_output=True, text=True, timeout=15)
        if 'not found' in result.stdout:
            raise ErrorDeUso('A required runtime library is missing.')
        for item in re.findall(r'(/[^\s()]+)', result.stdout):
            path = Path(item).resolve()
            if path.is_file():
                libraries[str(path)] = digest(path)
    return libraries


def init(input_path, state, pw_cmd='pw.x', checkpoint_seconds=900,
         grace_seconds=300, max_failures=3, threads=1, keep=2, runtime_id=None):
    if os.name != 'posix':
        raise ErrorDeUso('Resilient QE requires POSIX/Linux.')
    _positive(checkpoint_seconds, 'checkpoint-seconds')
    _positive(grace_seconds, 'grace-seconds')
    for value, name, minimum in ((max_failures, 'max-failures', 1), (threads, 'threads', 1), (keep, 'keep', 2)):
        if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
            raise ErrorDeUso(f'{name} must be an integer >= {minimum}.')
    input_path = Path(input_path).resolve()
    text = input_path.read_text()
    params, cards = _namelists(text)
    control = params.get('control', {})
    calculation = control.get('calculation', 'scf').lower()
    if calculation not in SUPPORTED:
        raise ErrorDeUso('Resilient mode supports pw.x scf, relax and vc-relax; other engines/modes require separate recovery protocols.')
    prefix = control.get('prefix', 'pwscf')
    if not re.fullmatch(r'[A-Za-z0-9_-]+', prefix):
        raise ErrorDeUso('prefix must contain only letters, digits, underscore or hyphen.')
    if control.get('restart_mode', 'from_scratch') != 'from_scratch':
        raise ErrorDeUso('Initialize from an original from_scratch input, not an unverified restart directory.')
    if control.get('disk_io', 'low') in ('none', 'nowf', 'minimal'):
        raise ErrorDeUso('disk_io must retain wavefunctions (low, medium or high).')
    pseudo_dir = Path(os.path.expandvars(os.path.expanduser(str(control.get('pseudo_dir', os.environ.get('ESPRESSO_PSEUDO', '~/espresso/pseudo'))))))
    if not pseudo_dir.is_absolute():
        pseudo_dir = input_path.parent / pseudo_dir
    try:
        i = next(i for i, line in enumerate(cards) if line.upper().startswith('ATOMIC_SPECIES'))
        n = int(params['system']['ntyp'])
        pseudos = [line.split()[2] for line in cards[i + 1:i + 1 + n]]
        if len(pseudos) != n or any(Path(p).name != p for p in pseudos):
            raise ValueError('invalid pseudopotential names')
    except (KeyError, IndexError, StopIteration, ValueError) as exc:
        raise ErrorDeUso('Cannot read ATOMIC_SPECIES and ntyp.') from exc
    for name in pseudos:
        if not (pseudo_dir / name).is_file():
            raise ErrorDeUso(f'Missing pseudopotential: {pseudo_dir / name}')
    command = shlex.split(pw_cmd) if isinstance(pw_cmd, str) else list(pw_cmd)
    if not command or any(x in command for x in ('-in', '-inp', '-input', '-i')):
        raise ErrorDeUso('pw-cmd must be an executable/launcher with parallel flags, without input arguments.')
    executable = shutil.which(command[0])
    if not executable:
        raise ErrorDeUso(f'Executable not found: {command[0]}')
    command[0] = str(Path(executable).resolve())
    # Freeze launcher and QE binaries; require absolute paths for custom MPI executables.
    binaries = {command[0]: digest(command[0])}
    for i, arg in enumerate(command[1:], 1):
        resolved = str(Path(arg).resolve()) if Path(arg).is_file() else shutil.which(arg)
        if resolved and Path(resolved).is_file():
            command[i] = str(Path(resolved).resolve())
            binaries[command[i]] = digest(command[i])
    if Path(command[0]).name in ('mpirun', 'mpiexec', 'orterun', 'srun') and len(binaries) < 2:
        raise ErrorDeUso('The QE executable behind MPI/Slurm must be resolvable.')
    root = Path(state).resolve()
    if root.exists():
        raise ErrorDeUso(f'State directory already exists: {root}; use run/status to resume.')
    root.parent.mkdir(parents=True, exist_ok=True)
    staging = root.with_name('.' + root.name + '.init-' + uuid.uuid4().hex)
    staging.mkdir()
    try:
        assets = staging / 'assets'
        (assets / 'pp').mkdir(parents=True)
        (assets / 'original.in').write_text(text)
        for name in set(pseudos):
            shutil.copyfile(pseudo_dir / name, assets / 'pp' / name)
        for p in assets.rglob('*'):
            if p.is_file():
                with p.open('rb') as f:
                    os.fsync(f.fileno())
        _sync_dir(assets / 'pp')
        _sync_dir(assets)
        job = dict(schema=SCHEMA, calculation=calculation, prefix=prefix, command=command,
                   binaries=binaries, libraries=runtime_libraries(binaries),
                   architecture=platform.machine(), runtime_id=runtime_id or platform.platform(),
                   environment=runtime_environment(), assets=files_tree(assets), checkpoint_seconds=float(checkpoint_seconds),
                   grace_seconds=float(grace_seconds), max_failures=max_failures,
                   keep=keep, threads=threads, created=time.time())
        atomic_json(staging / 'job.json', job)
        for folder in ('checkpoints', 'attempts'):
            (staging / folder).mkdir()
        atomic_json(staging / 'state.json', dict(status='pending', failures=0, attempts=0))
        _sync_dir(staging)
        os.rename(staging, root)
        _sync_dir(root.parent)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return job


def load_job(root):
    root = Path(root)
    job = _load(root / 'job.json')
    if job.get('schema') != SCHEMA:
        raise ErrorDeUso('Unsupported recovery schema.')
    if files_tree(root / 'assets') != job['assets']:
        raise ErrorDeUso('Input/pseudopotential fingerprint changed; recovery refused.')
    if job['architecture'] != platform.machine() or job['environment'] != runtime_environment():
        raise ErrorDeUso('Runtime architecture or environment changed; recovery refused.')
    for path, expected in {**job['binaries'], **job['libraries']}.items():
        if not Path(path).is_file() or digest(path) != expected:
            raise ErrorDeUso(f'QE/launcher binary changed or missing: {path}; restore the pinned environment.')
    return job


def checkpoint_info(work, job, returncode):
    """Classify clean stop versus physical success, never using JOB DONE alone."""
    work = Path(work)
    output = (work / 'pw.out').read_text(errors='replace')
    if returncode != 0 or 'JOB DONE.' not in output or re.search(r'Error in routine|convergence NOT achieved', output, re.I):
        raise ErrorDeUso('QE did not stop cleanly or failed convergence; checkpoint not published.')
    clean_stop = bool(re.search(r'Maximum (?:CPU|wall) time exceeded|stopped by user|stop requested|user request', output, re.I))
    save = work / 'out' / (job['prefix'] + '.save')
    try:
        tree = ET.parse(save / 'data-file-schema.xml')
    except (OSError, ET.ParseError) as exc:
        raise ErrorDeUso('Checkpoint XML missing or incomplete.') from exc
    values = {node.tag.split('}')[-1]: node.text for node in tree.iter()}
    charge = list(save.glob('charge-density.*'))
    wave = list((work / 'out').rglob('*wfc*'))
    if not any(p.is_file() and p.stat().st_size for p in charge) or not any(p.is_file() and p.stat().st_size for p in wave):
        raise ErrorDeUso('Charge density or wavefunctions missing from checkpoint.')
    scf = next((node for node in tree.iter() if node.tag.split('}')[-1] == 'scf_conv'), None)
    converged = scf is not None and any(node.tag.split('}')[-1] == 'convergence_achieved' and (node.text or '').strip().lower() == 'true' for node in scf)
    if job['calculation'] != 'scf':
        opt = next((node for node in tree.iter() if node.tag.split('}')[-1] == 'opt_conv'), None)
        converged = converged and opt is not None and any(node.tag.split('}')[-1] == 'convergence_achieved' and (node.text or '').strip().lower() == 'true' for node in opt)
    energy = None
    # A stop request can coincide with the final converged iteration.
    # XML convergence plus a finite energy remains a completed calculation.
    if converged:
        matches = re.findall(r'!\s+total energy\s*=\s*([-+\d.EeDd]+)\s+Ry', output)
        if matches:
            energy = float(matches[-1].replace('D', 'E').replace('d', 'e'))
        if energy is None or not math.isfinite(energy):
            raise ErrorDeUso('Converged calculation has no finite total energy.')
        status = 'succeeded'
    elif clean_stop:
        status = 'checkpoint'
    else:
        raise ErrorDeUso('Neither a completed calculation nor a recognized clean QE stop.')
    return {'status': status, 'energy_Ry': energy, 'n_scf_steps': values.get('n_scf_steps'),
            'qe_version': next((n.attrib.get('VERSION') for n in tree.iter() if n.tag.split('}')[-1] == 'creator'), None)}


def publish(root, work, job, info):
    root, work = Path(root), Path(work)
    if files_tree(root / 'assets') != job['assets']:
        raise ErrorDeUso('Assets changed while QE ran.')
    size = sum(p.stat().st_size for p in work.rglob('*') if p.is_file())
    if shutil.disk_usage(root).free < size + 16 * 1024 * 1024:
        raise ErrorDeUso('Insufficient disk space for an atomic checkpoint; previous generations retained.')
    # Wall clocks can move backwards after replacement/reboot. Publication order cannot.
    previous_ids = [int(p.name.split('-')[0]) for p in (root/'checkpoints').iterdir()
                    if re.fullmatch(r'[0-9]{20}-[0-9a-f]{8}', p.name)]
    generation_id = max([time.time_ns()-1, *previous_ids]) + 1
    name = f'{generation_id:020d}-{uuid.uuid4().hex[:8]}'
    stage = root / 'checkpoints' / ('.partial-' + name)
    stage.mkdir()
    copy_durable(work, stage / 'work')
    manifest = dict(schema=SCHEMA, job_hash=digest(root / 'job.json'), created=time.time(),
                    info=info, files=files_tree(stage / 'work'))
    atomic_json(stage / 'manifest.json', manifest)
    _sync_dir(stage)
    target = stage.with_name(name)
    os.rename(stage, target)
    _sync_dir(target.parent)
    return target


def _verified_generation(path, root):
    if path.is_symlink() or not path.is_dir():
        raise ErrorDeUso('Invalid generation directory.')
    manifest = _load(path / 'manifest.json')
    if (manifest['schema'] != SCHEMA or manifest['job_hash'] != digest(root / 'job.json')
            or not manifest['files'] or files_tree(path / 'work') != manifest['files']):
        raise ErrorDeUso('Fingerprint mismatch')
    expected = checkpoint_info(path / 'work', _load(root / 'job.json'), 0)
    if manifest['info'] != expected:
        raise ErrorDeUso('Checkpoint completion metadata does not match QE evidence.')
    return manifest


def latest(root):
    root = Path(root)
    bad = []
    generations = sorted((root / 'checkpoints').glob('[0-9]*'), reverse=True)
    for p in generations:
        try:
            return p, _verified_generation(p, root), bad
        except (ErrorDeUso, KeyError, TypeError, ValueError, OSError) as exc:
            bad.append({'generation': p.name, 'reason': str(exc)})
    if generations:
        raise ErrorDeUso('No intact checkpoint remains; refusing to overwrite or silently restart from scratch.')
    return None, None, bad


def _prune(root, job, protected):
    # Corrupt metadata must not count toward the verified backup retention minimum.
    good = []
    for p in sorted((root / 'checkpoints').glob('[0-9]*'), reverse=True):
        try:
            _verified_generation(p, root)
            good.append(p)
        except (ErrorDeUso, KeyError, TypeError, ValueError, OSError):
            pass
    for p in good[job['keep']:]:
        if p != protected:
            shutil.rmtree(p)
    for p in (root / 'checkpoints').glob('.partial-*'):
        if p.is_dir() and not p.is_symlink():
            shutil.rmtree(p)
    _sync_dir(root / 'checkpoints')


def _cleanup_attempts(root):
    for p in (root / 'attempts').iterdir():
        if not re.fullmatch(r'[0-9a-f]{32}', p.name) or p.is_symlink():
            raise ErrorDeUso('Unexpected attempt directory; cleanup refused.')
        w = p / 'work'
        if w.is_dir():
            for name in ('pw.in', 'pw.out', 'pw.stderr'):
                if (w / name).is_file():
                    shutil.copyfile(w / name, p / name)
            shutil.rmtree(w)


def _write_input(root, work, job, restarting):
    params, cards = _namelists((root / 'assets' / 'original.in').read_text())
    params.setdefault('control', {})
    params['control'].update(restart_mode='restart' if restarting else 'from_scratch',
                             max_seconds=job['checkpoint_seconds'], outdir='./out', wfcdir='./out',
                             pseudo_dir=str(root / 'assets' / 'pp'))
    (work / 'pw.in').write_text(params.to_string() + '\n' + '\n'.join(cards) + '\n')
    (work / (job['prefix'] + '.EXIT')).unlink(missing_ok=True)
    (work / 'out' / (job['prefix'] + '.EXIT')).unlink(missing_ok=True)


def pause(root):
    root = Path(root).resolve()
    if not (root / 'job.json').is_file():
        raise ErrorDeUso('Not a resilient job.')
    atomic_json(root / 'PAUSE.json', {'requested': time.time()})


def status(root):
    root = Path(root).resolve()
    job = load_job(root)
    generation, manifest, bad = latest(root)
    return {'state': _load(root / 'state.json'), 'calculation': job['calculation'],
            'checkpoint': generation.name if generation else None,
            'checkpoint_info': manifest['info'] if manifest else None,
            'rejected_generations': bad, 'pause_requested': (root / 'PAUSE.json').exists()}


def run(root, max_segments=0, resume=False):
    if isinstance(max_segments, bool) or not isinstance(max_segments, int) or max_segments < 0:
        raise ErrorDeUso('max-segments must be a nonnegative integer.')
    root = Path(root).resolve()
    with locked(root) as lock_fd:
        job = load_job(root)
        if resume:
            (root / 'PAUSE.json').unlink(missing_ok=True)
            _sync_dir(root)
        state = _load(root / 'state.json')
        # Holding the inherited lock proves no previous QE child still owns this job.
        state.pop('qe_pid', None)
        generation, manifest, bad = latest(root)
        if manifest and manifest['info']['status'] == 'succeeded':
            state.update(status='succeeded', checkpoint=generation.name, result=manifest['info'])
            atomic_json(root / 'state.json', state)
            return 0
        if state.get('status') == 'running':
            state['failures'] = state.get('failures', 0) + 1
            state['recovered_after_crash'] = True
            # Downtime and interrupted compute cannot be separated reliably after reboot.
            state['lost_compute_seconds'] = None
            if state.get('attempt') and not re.fullmatch(r'[0-9a-f]{32}', state['attempt']):
                raise ErrorDeUso('Invalid attempt reference in durable state.')
            if state.get('attempt') and (root/'attempts'/state['attempt']).is_dir():
                atomic_json(root/'attempts'/state['attempt']/'result.json',
                            {'status': 'interrupted', 'wall_seconds': None, 'qe_wall_seconds': None,
                             'restore_seconds': None, 'checkpoint_seconds': None, 'checkpoint_bytes': None})
        state['rejected_generations'] = bad
        _cleanup_attempts(root)
        requested = [False]
        old_handlers = {}
        def stop(signum, frame):
            requested[0] = True
        for sig in (signal.SIGTERM, signal.SIGINT):
            old_handlers[sig] = signal.signal(sig, stop)
        completed_segments = 0
        try:
            while True:
                if (root / 'PAUSE.json').exists() or requested[0]:
                    state['status'] = 'paused'
                    atomic_json(root / 'state.json', state)
                    return 75
                if state.get('failures', 0) >= job['max_failures']:
                    state['status'] = 'failed'
                    atomic_json(root / 'state.json', state)
                    return 2
                attempt_started = time.monotonic()
                restore_started = attempt_started
                attempt = root / 'attempts' / uuid.uuid4().hex
                attempt.mkdir()
                work = attempt / 'work'
                if generation:
                    copy_durable(generation / 'work', work)
                else:
                    work.mkdir()
                restore_seconds = time.monotonic() - restore_started
                _write_input(root, work, job, bool(generation))
                _sync_dir(attempt)
                _sync_dir(attempt.parent)
                # Restoration can be slow; honor a stop received while copying.
                if requested[0] or (root / 'PAUSE.json').exists():
                    state.update(status='paused', updated=time.time())
                    atomic_json(root / 'state.json', state)
                    return 75
                state.pop('qe_pid', None)
                state.update(status='running', attempts=state.get('attempts', 0) + 1,
                             attempt=attempt.name, restarted_from=generation.name if generation else None,
                             started=time.time())
                atomic_json(root / 'state.json', state)
                env = dict(os.environ)
                env.update({key: str(job['threads']) for key in THREAD_KEYS})
                env['ESPRESSO_TMPDIR'] = str(work / 'out')
                started = time.monotonic()
                stop_at = None
                forced = False
                with (work / 'pw.in').open() as fin, (work / 'pw.out').open('w') as fout, (work / 'pw.stderr').open('w') as ferr:
                    child = subprocess.Popen(job['command'], stdin=fin, stdout=fout, stderr=ferr,
                                             cwd=work, env=env, start_new_session=True, pass_fds=(lock_fd,))
                    try:
                        # Every operation after spawn belongs inside child cleanup.
                        state['qe_pid'] = child.pid
                        atomic_json(root / 'state.json', state)
                        while child.poll() is None:
                            stopping = requested[0] or (root / 'PAUSE.json').exists()
                            elapsed = time.monotonic() - started
                            if (stopping or elapsed >= job['checkpoint_seconds']) and stop_at is None:
                                (work / (job['prefix'] + '.EXIT')).touch()
                                stop_at = time.monotonic()
                            if stop_at is not None and time.monotonic() - stop_at > job['grace_seconds']:
                                os.killpg(child.pid, signal.SIGKILL)
                                forced = True
                                break
                            time.sleep(.1)
                        code = child.wait()
                    except BaseException:
                        with contextlib.suppress(ProcessLookupError):
                            os.killpg(child.pid, signal.SIGKILL)
                        child.wait()
                        raise
                    finally:
                        state.pop('qe_pid', None)
                    fout.flush()
                    ferr.flush()
                    os.fsync(fout.fileno())
                    os.fsync(ferr.fileno())
                qe_wall_seconds = time.monotonic() - started
                metrics = dict(restore_seconds=restore_seconds, qe_wall_seconds=qe_wall_seconds,
                               checkpoint_seconds=None, checkpoint_bytes=None)
                try:
                    if forced:
                        raise ErrorDeUso('Grace period expired; previous checkpoint retained.')
                    info = checkpoint_info(work, job, code)
                    tick = time.monotonic()
                    metrics['checkpoint_bytes'] = sum(p.stat().st_size for p in work.rglob('*') if p.is_file())
                    generation = publish(root, work, job, info)
                    metrics['checkpoint_seconds'] = time.monotonic() - tick
                    state.update(status=info['status'], checkpoint=generation.name, failures=0,
                                 result=info, updated=time.time())
                    atomic_json(root / 'state.json', state)
                    atomic_json(attempt / 'result.json', {'returncode': code, 'info': info, 'checkpoint': generation.name,
                                                         'wall_seconds': time.monotonic()-attempt_started, **metrics})
                    with contextlib.suppress(OSError):
                        _prune(root, job, generation)
                    completed_segments += 1
                    print(json.dumps({'checkpoint': generation.name, **info}), flush=True)
                    # Only logs survive outside committed checkpoints, bounding disk use.
                    with contextlib.suppress(OSError):
                        _cleanup_attempts(root)
                    if info['status'] == 'succeeded':
                        return 0
                    if max_segments and completed_segments >= max_segments:
                        return 75
                except (ErrorDeUso, OSError) as exc:
                    state.update(status='retrying', failures=state.get('failures', 0) + 1, error=str(exc))
                    atomic_json(root / 'state.json', state)
                    atomic_json(attempt / 'result.json', {'returncode': code, 'error': str(exc),
                                                         'wall_seconds': time.monotonic()-attempt_started, **metrics})
                    if requested[0] or (root / 'PAUSE.json').exists():
                        return 75
                    if state['failures'] < job['max_failures']:
                        time.sleep(min(5, state['failures']))
        finally:
            for sig, previous in old_handlers.items():
                signal.signal(sig, previous)


def service(root, output, user):
    """Generate, never install, a service for a retained disk attached to one worker."""
    root = Path(root).resolve()
    job = load_job(root)
    if not re.fullmatch(r'[a-z_][a-z0-9_-]*', user):
        raise ErrorDeUso('Provide an existing unprivileged Linux service user.')
    def quote(s, command_arg=False):
        s = str(s)
        if any(c in s for c in '\n\r\x00'):
            raise ErrorDeUso('Invalid service path.')
        s = s.replace('\\', '\\\\').replace('"', '\\"').replace('%', '%%')
        if command_arg:
            s = s.replace('$', '$$')
        return '"' + s + '"'
    command = ' '.join(quote(x, True) for x in (sys.executable, '-m', 'qekit', 'resilient', 'run', str(root)))
    if '\\' in str(root) or str(root) != str(root).rstrip():
        raise ErrorDeUso('Service state path cannot contain backslashes or trailing whitespace.')
    working_directory = str(root).replace('%', '%%')
    environment = '\n'.join('Environment=' + quote(k + '=' + v) for k, v in sorted(job['environment'].items()))
    text = f'''[Unit]
Description=Olla-DFT recoverable Quantum ESPRESSO job
After=local-fs.target
RequiresMountsFor={quote(root)}
StartLimitIntervalSec=0

[Service]
Type=simple
User={user}
WorkingDirectory={working_directory}
{environment}
ExecStart={command}
Restart=on-failure
RestartSec=30
RestartPreventExitStatus=1 2 75
KillMode=mixed
TimeoutStopSec=25
UMask=0077

[Install]
WantedBy=multi-user.target
'''
    Path(output).write_text(text)
    return text


def cli(args):
    if args.action == 'init':
        if not args.state:
            raise ErrorDeUso('init requires --state on a retained persistent disk.')
        init(args.target, args.state, args.pw_cmd, args.checkpoint_seconds,
             args.grace_seconds, args.max_failures, args.threads, args.keep, args.runtime_id)
        print(f'Recoverable job created: {Path(args.state).resolve()}')
        return 0
    if args.action == 'run':
        if args.max_segments < 0:
            raise ErrorDeUso('max-segments must be nonnegative.')
        try:
            return run(args.target, args.max_segments, args.resume)
        except BusyJob as exc:
            print(str(exc), file=sys.stderr)
            return 76
    if args.action == 'pause':
        pause(args.target)
        return 0
    if args.action == 'status':
        print(json.dumps(status(args.target), indent=2))
        return 0
    if args.action == 'service':
        if not args.output or not args.user:
            raise ErrorDeUso('service requires --output and --user.')
        service(args.target, args.output, args.user)
        print(args.output)
        return 0
    raise ErrorDeUso('Unknown resilient action.')
