#!/usr/bin/env python3
"""Opt-in real QE crash/recovery comparison. Produces auditable evidence, no cloud resources."""
import argparse
import hashlib
import json
import os
import re
from pathlib import Path
import signal
import subprocess
import sys
import time
from qekit.modules import resilient as r


LIVE = []


def launch(root, *args):
    with (root/'controller.log').open('a') as log:
        process = subprocess.Popen([sys.executable, '-m', 'qekit', 'resilient', 'run', str(root), *args],
                                   stdout=log, stderr=subprocess.STDOUT)
    LIVE.append((root, process))
    return process


def cleanup():
    # A failed test or timeout must not leave a costly compute process running.
    for root, process in LIVE:
        if process.poll() is not None:
            continue
        process.terminate()
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            state=json.loads((root/'state.json').read_text())
            if state.get('qe_pid'):
                try: os.killpg(state['qe_pid'], signal.SIGKILL)
                except ProcessLookupError: pass
            process.kill()
            process.wait(timeout=5)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--pw',required=True);ap.add_argument('--out',required=True);a=ap.parse_args()
    out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=False)
    r.init(a.input,out/'baseline',a.pw,checkpoint_seconds=600,grace_seconds=30)
    baseline=launch(out/'baseline');assert baseline.wait(timeout=120)==0
    r.init(a.input,out/'recovered',a.pw,checkpoint_seconds=60,grace_seconds=30)
    root=out/'recovered'
    first=launch(root,'--max-segments','1')
    for _ in range(6000):
        first_state=json.loads((root/'state.json').read_text())
        attempt_log=root/'attempts'/first_state.get('attempt','missing')/'work/pw.out'
        if attempt_log.exists() and re.search(r'iteration\s*#\s*2', attempt_log.read_text(errors='replace')):
            first.send_signal(signal.SIGTERM)
            break
        time.sleep(.01)
    assert first.wait(timeout=120)==75
    previous,manifest,_=r.latest(root)
    assert previous is not None, 'Input stopped before a recoverable checkpoint was available'
    before=hashlib.sha256((previous/'manifest.json').read_bytes()).hexdigest()
    # Kill both the supervisor and its compute process after the attempt is durable.
    child=launch(root)
    work=None
    for _ in range(6000):
        state=json.loads((root/'state.json').read_text())
        if state['status']=='running':
            work=root/'attempts'/state['attempt']/'work'
            pids=[str(state['qe_pid'])] if state.get('qe_pid') and state['qe_pid'] != first_state.get('qe_pid') else []
            if pids:
                break
        time.sleep(.01)
    else:
        raise AssertionError('QE child never started')
    os.kill(child.pid,signal.SIGKILL)
    for pid in pids:
        try:os.killpg(int(pid),signal.SIGKILL)
        except ProcessLookupError:pass
    assert child.wait(timeout=10)<0
    # Simulate partial/corrupt writes in the crashed attempt, never the backup.
    if work:
        (work/'out').mkdir(exist_ok=True)
        (work/'out'/'torn-write').write_bytes(b'incomplete')
    assert hashlib.sha256((previous/'manifest.json').read_bytes()).hexdigest()==before
    retry=launch(root);assert retry.wait(timeout=180)==0
    result=r.status(root);reference=r.status(out/'baseline')
    delta=abs(result['checkpoint_info']['energy_Ry']-reference['checkpoint_info']['energy_Ry'])
    assert delta<1e-7
    assert result['state'].get('recovered_after_crash') is True
    logs='\n'.join(p.read_text(errors='replace') for p in (root/'attempts').glob('*/pw.out'))
    assert 'Calculation restarted' in logs or 'Starting wfcs from file' in logs
    evidence={'baseline':reference,'recovered':result,'energy_difference_Ry':delta,
              'kill':'SIGKILL supervisor and QE process group','checkpoint_unchanged_after_kill':True,
              'resumed_from_files':True,'cloud_preemption_tested':False}
    (out/'evidence.json').write_text(json.dumps(evidence,indent=2)+'\n')
    print(json.dumps({'energy_difference_Ry':delta,'attempts':result['state']['attempts'],'evidence':str(out/'evidence.json')}))

if __name__=='__main__':
    try:
        main()
    finally:
        cleanup()
