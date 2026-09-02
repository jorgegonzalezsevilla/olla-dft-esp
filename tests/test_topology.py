from pathlib import Path

import numpy as np
import pytest

from qekit.cli import ALIASES, COMMAND_GROUPS, _DISPATCH, build_parser, main
from qekit.core.errors import ErrorDeUso
from qekit.modules import topology


SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)


def _qwz_terms(mass):
    """Qi-Wu-Zhang: sin(kx)sx + sin(ky)sy + (m+cos kx+cos ky)sz."""
    R = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0],
                  [0, 1, 0], [0, -1, 0]])
    H = np.array([
        mass * SZ,
        -0.5j * SX + 0.5 * SZ,
        +0.5j * SX + 0.5 * SZ,
        -0.5j * SY + 0.5 * SZ,
        +0.5j * SY + 0.5 * SZ,
    ])
    return H, R, np.ones(len(R))


def _write_hr(path, H, R, deg):
    path = Path(path)
    nw, nr = H.shape[1], len(R)
    lines = ["# test Hamiltonian", f"{nw:12d}", f"{nr:12d}",
             "".join(f"{int(d):5d}" for d in deg)]
    for ir in range(nr):
        for n in range(nw):
            for m in range(nw):
                z = H[ir, m, n]
                lines.append(
                    f"{R[ir, 0]:5d}{R[ir, 1]:5d}{R[ir, 2]:5d}"
                    f"{m + 1:5d}{n + 1:5d}{z.real:14.8f}{z.imag:14.8f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


@pytest.mark.parametrize("mass, expected", [(-1.0, -1), (1.0, 1), (3.0, 0)])
def test_qwz_known_chern_phases(tmp_path, mass, expected):
    model = _write_hr(tmp_path / "QWZ_hr.dat", *_qwz_terms(mass))
    run = topology.analyze(model, occupied=1, grid=(30, 30))
    assert run.chern == expected
    assert run.chern_raw == pytest.approx(expected, abs=1e-12)
    assert run.direct_gap > 1.9
    assert run.min_overlap > 0.99


def test_non_abelian_result_is_invariant_under_local_gauge():
    n = 16
    base = np.empty((n, n, 4, 2), complex)
    for i in range(n):
        for j in range(n):
            x, y = 2 * np.pi * i / n, 2 * np.pi * j / n
            H = (np.sin(x) * SX + np.sin(y) * SY
                 + (-1 + np.cos(x) + np.cos(y)) * SZ)
            _e, v = np.linalg.eigh(H)
            base[i, j, :2, 0] = v[:, 0]
            base[i, j, 2:, 0] = 0
            base[i, j, :2, 1] = 0
            base[i, j, 2:, 1] = v[:, 0]
    flux0, chern0, wilson0, _ = topology.invariants_from_vectors(base)

    rng = np.random.default_rng(90210)
    gauged = base.copy()
    for i in range(n):
        for j in range(n):
            z = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
            q, r = np.linalg.qr(z)
            q = q @ np.diag(np.diag(r) / np.abs(np.diag(r)))
            gauged[i, j] = gauged[i, j] @ q
    flux1, chern1, wilson1, _ = topology.invariants_from_vectors(gauged)

    assert chern0 == pytest.approx(-2.0, abs=1e-12)
    assert chern1 == pytest.approx(chern0, abs=1e-12)
    assert flux1 == pytest.approx(flux0, abs=1e-12)
    circular_difference = (wilson1 - wilson0 + 0.5) % 1.0 - 0.5
    assert circular_difference == pytest.approx(0.0, abs=1e-12)


def test_gap_closing_is_rejected(tmp_path):
    model = _write_hr(tmp_path / "critical_hr.dat", *_qwz_terms(0.0))
    with pytest.raises(ErrorDeUso, match="no está aislado"):
        topology.analyze(model, occupied=1, grid=(40, 40))


def test_occupation_is_explicit_and_model_resolution_is_unambiguous(tmp_path):
    model = _write_hr(tmp_path / "WANNIER_hr.dat", *_qwz_terms(3.0))
    with pytest.raises(ErrorDeUso, match="exactamente una"):
        topology.analyze(model)
    run = topology.analyze(tmp_path, fermi=0.0, grid=(12, 12))
    assert run.model_path == str(model.resolve())
    assert run.occupied == 1

    model.unlink()
    _write_hr(tmp_path / "a_hr.dat", *_qwz_terms(3.0))
    _write_hr(tmp_path / "b_hr.dat", *_qwz_terms(3.0))
    with pytest.raises(ErrorDeUso, match="varios modelos"):
        topology.resolve_model(tmp_path)


def test_export_and_cli(tmp_path):
    model = _write_hr(tmp_path / "QWZ_hr.dat", *_qwz_terms(-1.0))
    out = tmp_path / "result"
    code = main(["topology", str(model), "--occupied", "1", "--grid",
                 "12x12", "--no-plot", "-o", str(out)])
    assert code == 0
    assert (out / "TOPOLOGY.txt").is_file()
    assert "Chern entero:       -1" in (out / "TOPOLOGY.txt").read_text()
    curvature = np.loadtxt(out / "TOPOLOGY_curvature.dat")
    wilson = np.loadtxt(out / "TOPOLOGY_wilson.dat")
    assert curvature.shape == (12 * 12, 3)
    assert wilson.shape == (12, 2)

    run = topology.analyze(model, occupied=1, grid=(12, 12))
    figures = topology.plot(run, str(out / "figure"), formats="png", dpi=72)
    assert len(figures) == 1
    assert Path(figures[0]).is_file()


def test_cli_catalog_covers_every_command_once(capsys):
    parser = build_parser()
    commands = set(parser._subparsers._group_actions[0].choices) - set(ALIASES)
    grouped = [name for _title, names in COMMAND_GROUPS for name in names]
    assert len(grouped) == len(set(grouped))
    assert set(grouped) == commands == set(_DISPATCH)
    assert set(ALIASES.values()) <= commands
    parser.print_help()
    help_text = capsys.readouterr().out
    assert "COMANDO ..." in help_text
    assert "Estructura electrónica" in help_text
    assert "topology" in help_text


def test_el_reporte_no_mezcla_idiomas(tmp_path):
    """Salía "Gap indirect" entre líneas en español."""
    model = _write_hr(tmp_path / "QWZ_hr.dat", *_qwz_terms(-1.0))
    run = topology.analyze(model, occupied=1, grid=(12, 12))
    rep = topology.report(run)
    assert "Gap indirecto:" in rep
    assert "Gap directo mínimo:" in rep
    assert "indirect " not in rep and "Gap indirect:" not in rep
