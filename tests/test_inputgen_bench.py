"""Cambios motivados por olla-dft-bench: malla explícita y mixing_beta por tipo de ocupación."""
from qekit.core import structure
from qekit.modules import inputgen


def _pw(insulator):
    atoms = structure.load("examples/ZnO.cif")
    pseudos = {"Zn": {"filename": "Zn.UPF"}, "O": {"filename": "O.UPF"}}
    return inputgen.build_pw_input(atoms, pseudos, "scf", "zno", "/pp",
                                   30.0, 240.0, inputgen._kgrid_card((6, 6, 4)), insulator=insulator)


def test_mixing_beta_depende_de_las_ocupaciones():
    assert "0.7" in [l for l in _pw(True).splitlines() if "mixing_beta" in l][0]
    assert "0.4" in [l for l in _pw(False).splitlines() if "mixing_beta" in l][0]


def test_kgrid_explicita_anula_kspacing():
    opts = inputgen.GenOptions(kgrid=(5, 5, 3), kspacing=0.9)
    assert tuple(opts.kgrid) == (5, 5, 3)
    card = inputgen._kgrid_card((5, 5, 3))
    assert card.split()[-6:-3] == ["5", "5", "3"]
