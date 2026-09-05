# Recovery contract: Olla-DFT and Olla-Lungo

Olla-DFT owns the recovery engine in `qekit.modules.resilient`. Olla-Lungo
owns cost estimation, orchestration and comparison of observed results. It
must call this engine rather than maintain a second restart implementation.
The interface is provisional until the 1.2.0 release review is complete.

## Python API and CLI

```python
from qekit.modules import resilient
resilient.init(input_path, state, pw_cmd="/opt/qe/bin/pw.x",
               checkpoint_seconds=900, grace_seconds=300, max_failures=3,
               threads=1, keep=2, runtime_id="immutable-image-build-id")
code = resilient.run(state, max_segments=0, resume=False)
report = resilient.status(state)
```

Use keyword arguments. `init` runs once, copies input/UPF assets and refuses an
existing state directory. It does not start QE. `run` uses the frozen job; it
must run on the main thread because it handles signals. `status` verifies
assets/runtime/checksums and returns a dictionary. `pause` requests a clean
persistent pause; `run(..., resume=True)` clears that manual pause. SIGTERM
requests a clean stop without a persistent pause, allowing recovery on boot.
No function provisions resources or sends data to the cloud.

CLI equivalent:

```sh
olla-dft resilient init scf.in --state /mnt/olla/job-001 \
  --pw-cmd '/opt/qe/bin/pw.x' --runtime-id IMAGE_BUILD_ID
olla-dft resilient run /mnt/olla/job-001
olla-dft resilient status /mnt/olla/job-001
```

Return codes: 0 = verified physical success; 75 = stopped/paused/segment limit;
2 = exhausted retry budget or invalid state; 76 = another worker/QE child
holds the lock (CLI only). Python raises `BusyJob` for 76 and `ErrorDeUso` for
invalid input/state. A nonzero code must never count as scientific success.
`max_segments=0` has no segment limit; positive values bound clean segments
in this invocation, not lifetime spending. Permanent failure attempts persist
across worker restarts. Very short segments can spend all their time starting
QE: callers must enforce their own wall-time and monetary limits.

## Identity and storage

Schema 1 records original input and UPF hashes, command and fixed MPI flags,
thread count, binary and linked-library hashes, architecture, relevant runtime
environment and an image/build label. The label is provenance, not proof that
a VM actually booted that image. Restore the pinned image and paths on a
replacement VM. No checkpoint format migration is implicit. Dynamic plugins,
GPU drivers and filesystem guarantees are outside the library-hash check;
these require a fixed environment and separate validation.

The state root must be a retained POSIX filesystem attached to exactly one VM
for writing, supporting flock, fsync and atomic rename. Local SSD and bucket
FUSE mounts do not satisfy this contract. Mount the persistent disk before
starting the worker; retain it when deleting/replacing the VM.

Each attempt restores a committed snapshot into a private writable directory.
Only return code 0, a recognized clean QE stop, valid XML and nonempty charge
and wavefunction files can become a checkpoint. Success additionally requires
SCF convergence, finite energy, and ionic convergence for relax/vc-relax.
A stop marker coinciding with convergence does not override verified success.
Two verified generations are retained by default. Partial publications are
ignored; a corrupt newest generation falls back to an older intact one. If
all committed generations fail validation, recovery refuses to start fresh.
The surviving checkpoint can precede the interruption; unfinished work may be
lost. Neither filesystem checksums nor this protocol prove physical accuracy.

## Results and measurement

`status()['checkpoint_info']` contains `status`, `energy_Ry`, `n_scf_steps` and
`qe_version`; no checkpoint means null. `status()['state']` includes attempt
count, retry failures, last generation and whether an abrupt crash was seen.
The committed directory is `checkpoints/<generation>/work/`; original outputs
and XML remain available there. Archived attempt inputs/stdout/stderr are in
`attempts/<id>/`; abandoned mutable workspaces are removed under the lock.

`attempts/<id>/result.json` records `wall_seconds`, `qe_wall_seconds`,
`restore_seconds`, `checkpoint_seconds`, `checkpoint_bytes`, return code and
success info or error. Times use a monotonic clock. Publication includes copy,
hashing and syncing; restore includes copying. Total attempt time excludes
later retention cleanup. An abrupt crash records unknown timings as null;
never treat missing measurements as zero. Logs and rejected corrupt
generations are retained for diagnosis and need an external retention policy.

To compare scientific equivalence, require identical original input/assets,
matching runtime/parallelization, both runs physically converged and explicit
tolerances and units for energy, every force component and every stress
component. Missing requested observables must fail the comparison. Generate
both runs with the same `tprnfor`/`tstress` values; do not silently enable them
only on one side. A silicon energy-only test cannot certify force/stress
recovery, different materials or every supported calculation mode.

## Evidence boundary

Local QE 7.4 silicon SCF validation killed both supervisor and QE, corrupted an
abandoned workspace and recovered the last committed checkpoint. The printed
final energy matched an uninterrupted run: -22.83929159 Ry (difference 0 Ry
at printed precision). This does not establish bitwise internal equivalence,
cloud availability, savings, MPI restart support on every deployment, or
recovery after loss of the persistent disk itself. The opt-in reproduction
script is `validation/resilience/exercise.py`. Google Cloud replacement and
billing must be measured separately before making a cost claim.

Additional installed-wheel tests passed for displaced-silicon `relax` and
`vc-relax`. Their XML comparison includes energy, every force/stress component,
atom ordering, final positions and cell vectors. `relax` differences were zero.
For `vc-relax`, maximum differences were 1.88e-11 Ry (energy), 2.27e-8 Ry/bohr
(force), 1.09e-9 Ry/bohr^3 (stress), 1.06e-7 bohr (positions), and 2.79e-7 bohr
(cell). These passed fixed tolerances of 1e-7, 1e-6, 1e-8, 1e-6 and 1e-6 in
the respective units. Reproduce using `validation/resilience/compare_geometry.py`
on the two independent completed job directories. These are small serial
cases, not a certification of all materials, MPI layouts or QE builds.

### Small lifecycle fixes after focused Fable review

A pause or SIGTERM received during restoration is rechecked before publishing
`running` or incrementing the launch count. The private unlaunched workspace
is removed by the existing cleanup on the next invocation; repeated pauses
do not accumulate writable workspaces. The last committed generation is intact.

After Popen succeeds, recording the PID is inside the child-cleanup exception
handler. Failure to persist that record kills/reaps QE without requiring another
disk write. A normal completed/checkpoint/failed result no longer retains
`qe_pid`. An abrupt crash before the next state write can still leave historical
`running`/PID metadata: status alone is not proof that that process is alive.
Acquiring the inherited worker lock, not a PID lookup, governs recovery.
