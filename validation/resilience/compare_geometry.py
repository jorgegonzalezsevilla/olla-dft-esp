#!/usr/bin/env python3
"""Opt-in restart equivalence gate including final geometry (QE XML atomic units)."""
import argparse
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET
from qekit.modules import resilient


def observables(root):
    job = resilient.load_job(root)
    generation, manifest, _ = resilient.latest(root)
    if manifest is None or manifest['info']['status'] != 'succeeded':
        raise ValueError('A verified completed calculation is required.')
    xml = generation/'work/out'/(job['prefix']+'.save')/'data-file-schema.xml'
    output = ET.parse(xml).getroot().find('output')
    structure = output.find('atomic_structure')

    def values(node):
        if node is None:
            raise ValueError('Missing required observable; enable tprnfor/tstress in both original inputs.')
        result = [float(v) for text in node.itertext() for v in text.split()]
        if not result or not all(math.isfinite(v) for v in result):
            raise ValueError('Missing or nonfinite observable.')
        return result

    # QE XML uses Hartree atomic units. Energy, force and stress gain a factor 2;
    # positions/cell already use bohr. Do not apply text-output kbar conversions.
    result = {
        'energy_Ry': [2*v for v in values(output.find('total_energy/etot'))],
        'forces_Ry_bohr': [2*v for v in values(output.find('forces'))],
        'stress_Ry_bohr3': [2*v for v in values(output.find('stress'))],
        'positions_bohr': values(structure.find('atomic_positions')),
        'cell_bohr': values(structure.find('cell')),
        'species': [a.attrib['name'] for a in structure.find('atomic_positions')],
    }
    nat = len(result['species'])
    if (len(result['energy_Ry']) != 1 or len(result['cell_bohr']) != 9
            or len(result['stress_Ry_bohr3']) != 9 or nat < 1
            or len(result['positions_bohr']) != 3*nat or len(result['forces_Ry_bohr']) != 3*nat):
        raise ValueError('Observable dimensions do not match the atomic structure.')
    return job, result, generation


def compare(baseline, candidate):
    baseline, candidate = Path(baseline).resolve(), Path(candidate).resolve()
    if baseline == candidate or baseline.samefile(candidate):
        raise ValueError('Independent jobs are required; a job cannot validate itself.')
    ja, a, ga = observables(baseline)
    jb, b, gb = observables(candidate)
    if ga.samefile(gb):
        raise ValueError('Independent checkpoint generations are required.')
    for key in ('assets', 'calculation', 'command', 'binaries', 'libraries', 'architecture', 'environment', 'threads'):
        if ja[key] != jb[key]:
            raise ValueError('Restart comparison requires identical physics/runtime: '+key)
    if a['species'] != b['species'] or any(len(a[k]) != len(b[k]) for k in a):
        raise ValueError('Species, atom ordering or observable dimensions differ.')
    delta = {k: max(abs(x-y) for x, y in zip(a[k], b[k])) for k in a if k != 'species'}
    tolerances = dict(energy_Ry=1e-7, forces_Ry_bohr=1e-6, stress_Ry_bohr3=1e-8,
                      positions_bohr=1e-6, cell_bohr=1e-6)
    return dict(calculation=ja['calculation'], baseline_generation=str(ga), candidate_generation=str(gb), source='QE XML, Hartree converted to Rydberg',
                max_absolute_differences=delta, tolerances=tolerances,
                passed=all(delta[k] <= tolerances[k] for k in delta), baseline=a, recovered=b)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline')
    parser.add_argument('candidate')
    args = parser.parse_args()
    try:
        result = compare(args.baseline, args.candidate)
        print(json.dumps(result, indent=2, allow_nan=False))
        return 0 if result['passed'] else 2
    except (ValueError, OSError, ET.ParseError, AttributeError, KeyError, resilient.ErrorDeUso) as exc:
        print(json.dumps({'passed': False, 'error': str(exc)}))
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
