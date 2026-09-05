"""Regression cases found by the September 2026 benchmark audit."""
import numpy as np
import pytest
from ase.build import bulk
from qekit.core.qeout import QEResult
from qekit.core.errors import ErrorDeUso, FaltanDatos
from qekit.modules import bands, inputgen


def structure(empty=False, converged=True, fermi=-0.8):
    energies = [[-2., -1.], [-1.8, -.8]]
    if empty:
        energies = [row + [1. + i] for i, row in enumerate(energies)]
    return bands.BandStructure(QEResult(
        eigenvalues=np.array([energies]), nbnd=len(energies[0]), nelec=4,
        nspin=1, kpoints_frac=np.zeros((2, 3)), fermi=fermi,
        occupations_kind='fixed', converged=converged))


@pytest.mark.parametrize('fermi', [-.8, None])
def test_occupied_only_is_insufficient_not_metal(fermi):
    bs = structure(fermi=fermi)
    info = bands.analyze_gap(bs)
    assert not info.is_metal and info.gap is None
    assert 'No hay bandas de conducción' in bands.gap_report(bs)
    assert 'METÁLICO' not in bands.gap_report(bs)


def test_gap_still_detects_insulator_and_crossing():
    assert bands.analyze_gap(structure(empty=True)).gap == pytest.approx(1.8)
    bs = structure(empty=True, fermi=0.)
    bs.result.eigenvalues[0, :, 1] = [-.2, .2]
    assert bands.analyze_gap(bs).is_metal


@pytest.mark.parametrize('converged,empty,expected', [(False, True, 2), (True, False, 2), (True, True, 0)])
def test_gap_cli_status_and_convergence_warning(monkeypatch, capsys, converged, empty, expected):
    from qekit import cli
    monkeypatch.setattr(bands, 'load', lambda *a, **k: structure(empty=empty, converged=converged))
    assert cli.main(['gap', '.']) == expected
    text = capsys.readouterr().out
    assert ('no convergió' in text) == (not converged)


@pytest.mark.parametrize('field,value', [
    ('kgrid', (0, -2, 4)), ('kgrid', (4, 4)), ('kgrid', (4., 4, 4)),
    ('ecutwfc', -30), ('ecutwfc', 0), ('ecutrho', -240),
    ('ecutrho', float('nan')), ('ecutwfc', float('inf')),
    ('kspacing', -.2), ('kspacing_nscf', float('nan')),
])
def test_invalid_generation_is_rejected_before_writing(tmp_path, field, value):
    target = tmp_path / 'must-not-exist'
    opts = inputgen.GenOptions(preset='scf', outdir=str(target), **{field: value})
    with pytest.raises(ErrorDeUso):
        inputgen.generate(bulk('Si', 'diamond', 5.43), opts)
    assert not target.exists()


def test_nonfinite_bands_are_not_a_metal():
    bs = structure()
    bs.result.eigenvalues[0, 0, 0] = np.nan
    with pytest.raises(FaltanDatos):
        bands.analyze_gap(bs)
