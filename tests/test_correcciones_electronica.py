"""Regresiones de los módulos de estructura electrónica (v0.35).

Cada test documenta un fallo concreto que se vio y se corrigió; si vuelve a
fallar, el mensaje del docstring dice qué se rompió.
"""

import pytest

from qekit.core import wfc
from qekit.modules import unfold


# ----------------------------------------------------------------------
# buscar_wfc: nunca mezclar canales de espín
# ----------------------------------------------------------------------
def _toca(carpeta, *nombres):
    carpeta.mkdir(parents=True, exist_ok=True)
    for n in nombres:
        (carpeta / n).write_bytes(b"")


def test_buscar_wfc_en_lsda_devuelve_solo_el_canal_pedido(tmp_path):
    """Con lsda, `wfc*.dat` ordenado por dígitos intercalaba up1, dw1, up2...

    La mitad de los puntos k del desdoblamiento salían del canal contrario.
    """
    save = tmp_path / "si.save"
    _toca(save, "wfcup1.dat", "wfcdw1.dat", "wfcup2.dat", "wfcdw2.dat",
          "wfcup10.dat", "wfcdw10.dat")
    up = [f.name for f in wfc.buscar_wfc(save)]
    assert up == ["wfcup1.dat", "wfcup2.dat", "wfcup10.dat"]
    dw = [f.name for f in wfc.buscar_wfc(save, spin="dw")]
    assert dw == ["wfcdw1.dat", "wfcdw2.dat", "wfcdw10.dat"]
    assert [f.name for f in wfc.buscar_wfc(save, spin="down")] == dw
    assert wfc.es_lsda(save)


def test_buscar_wfc_sin_espin_ignora_el_selector(tmp_path):
    save = tmp_path / "si.save"
    _toca(save, "wfc1.dat", "wfc2.dat", "wfc3.dat", "wfc12.dat",
          "charge-density.dat")
    esperado = ["wfc1.dat", "wfc2.dat", "wfc3.dat", "wfc12.dat"]
    assert [f.name for f in wfc.buscar_wfc(save)] == esperado
    assert [f.name for f in wfc.buscar_wfc(save, spin="dw")] == esperado
    assert not wfc.es_lsda(save)


def test_buscar_wfc_rechaza_un_canal_que_no_existe(tmp_path):
    save = tmp_path / "si.save"
    _toca(save, "wfc1.dat")
    with pytest.raises(ValueError):
        wfc.buscar_wfc(save, spin="lateral")


def test_el_aviso_de_lsda_del_desdoblamiento_dice_que_canal_es():
    txt = unfold.aviso_lsda("up", "dw")
    assert "lsda" in txt
    assert "'up'" in txt
    assert "--spin dw" in txt


def test_el_reporte_del_desdoblamiento_lleva_el_aviso_de_lsda():
    import numpy as np
    d = unfold.Desdoblado(kpath=np.zeros((2, 3)), distancias=np.zeros(2),
                          energias=np.zeros((2, 3)), pesos=np.ones((2, 3)),
                          M=np.eye(3, dtype=int), ncel=1, spin="up")
    d.avisos.append(unfold.aviso_lsda("up", "dw"))
    rep = unfold.report(d)
    assert "espín polarizado" in rep and "--spin dw" in rep


# ----------------------------------------------------------------------
# dos.load: un pdos_atm# con otra malla se salta, pero se AVISA
# ----------------------------------------------------------------------
def _pdos_falso(carpeta, nombre, n, l="s"):
    """Un archivo de projwfc.x con `n` puntos: E, ldos, (2l+1) pdos."""
    import numpy as np
    nm = 2 * {"s": 0, "p": 1, "d": 2}[l] + 1
    e = np.linspace(-5.0, 5.0, n)
    cols = [e, np.full(n, 0.5)] + [np.full(n, 0.5 / nm)] * nm
    np.savetxt(carpeta / nombre, np.column_stack(cols),
               header="E (eV)  ldos(E)  pdos(E)")


def test_una_pdos_con_otra_malla_se_salta_y_se_avisa(tmp_path):
    """Antes el `continue` era silencioso y la PDOS salía sin ese orbital."""
    from qekit.modules import dos
    _pdos_falso(tmp_path, "si.pdos.pdos_atm#1(Si)_wfc#1(s)", 101)
    _pdos_falso(tmp_path, "si.pdos.pdos_atm#1(Si)_wfc#2(p)", 101, l="p")
    _pdos_falso(tmp_path, "si.pdos.pdos_atm#2(Si)_wfc#1(s)", 57)   # vieja
    dd = dos.load(str(tmp_path))
    assert ("Si", "s") in dd.projected and ("Si", "p") in dd.projected
    assert len(dd.avisos) == 1
    assert "SALTADO" in dd.avisos[0]
    assert "pdos_atm#2(Si)_wfc#1(s)" in dd.avisos[0]
    assert "57 puntos" in dd.avisos[0]
    rep = dos.report(dd)
    assert "AVISO:" in rep and "pdos_atm#2(Si)_wfc#1(s)" in rep


def test_sin_mallas_distintas_no_hay_aviso(tmp_path):
    from qekit.modules import dos
    _pdos_falso(tmp_path, "si.pdos.pdos_atm#1(Si)_wfc#1(s)", 101)
    _pdos_falso(tmp_path, "si.pdos.pdos_atm#2(Si)_wfc#1(s)", 101)
    dd = dos.load(str(tmp_path))
    assert dd.avisos == []
    assert "AVISO" not in dos.report(dd)


# ----------------------------------------------------------------------
# effmass: la ventana por omisión cabe dentro del régimen parabólico
# ----------------------------------------------------------------------
def _bandas_parabolicas(nk=201, kmax=0.5, m_e=0.5, m_h=0.4, gap=1.0):
    """Un camino 1D fino y dos bandas parabólicas exactas en torno a k = 0."""
    import numpy as np
    from qekit.core.qeout import QEResult
    from qekit.modules import bands as bm
    from qekit.modules.effmass import HBAR2_OVER_ME
    k = np.linspace(-kmax, kmax, nk)
    vb = -HBAR2_OVER_ME * k ** 2 / (2 * m_h)
    cb = gap + HBAR2_OVER_ME * k ** 2 / (2 * m_e)
    res = QEResult()
    res.nbnd, res.nelec, res.nspin = 2, 2.0, 1
    res.eigenvalues = np.stack([vb, cb], axis=1)[None, :, :]
    res.kpoints_cart = np.column_stack([k, np.zeros(nk), np.zeros(nk)])
    res.kpoints_frac = res.kpoints_cart.copy()
    res.fermi = gap / 2
    return bm.BandStructure(result=res, kdist=k - k[0], labels=[], breaks=[])


def test_la_ventana_por_omision_esta_dentro_del_regimen_parabolico():
    """--window 0.15 (±0.15, tramo 0.30) contra un límite de 0.12 avisaba
    SIEMPRE, incluso con un camino finísimo: se comparaba un semiancho con
    una extensión total. Ahora las dos son extensiones y el valor por
    omisión cabe dentro del límite.
    """
    from qekit.modules import effmass as em
    assert 2 * em.WINDOW_DEFAULT <= em.PARABOLIC_MAX + 1e-12
    bs = _bandas_parabolicas()
    run = em.from_bands(bs)                       # ventana por omisión
    assert not run.is_metal and len(run.fits) == 2
    for f in run.fits:
        assert f.warning == "", f.warning
        assert f.window <= em.PARABOLIC_MAX + 1e-6
    masas = {f.carrier: abs(f.mass) for f in run.fits}
    assert masas["electrón"] == pytest.approx(0.5, rel=1e-6)
    assert masas["hueco"] == pytest.approx(0.4, rel=1e-6)


def test_una_ventana_mas_ancha_que_el_limite_si_avisa():
    from qekit.modules import effmass as em
    bs = _bandas_parabolicas()
    run = em.from_bands(bs, window=0.2)         # tramo 0.40 > 0.12
    assert run.fits
    for f in run.fits:
        assert "fuera del régimen parabólico" in f.warning \
            or "límite parabólico" in f.warning
        assert f"{f.window:.3f}" in f.warning


def test_el_parser_de_effmass_deja_la_ventana_al_modulo():
    """El valor por omisión vive en effmass.WINDOW_DEFAULT, no en cli.py."""
    from qekit import cli
    parser = cli.build_parser()
    args = parser.parse_args(["effmass", "Si.cif", "--bands-dir", "x"])
    assert args.window is None
    args = parser.parse_args(["effmass", "Si.cif", "--window", "0.03"])
    assert args.window == pytest.approx(0.03)


# ----------------------------------------------------------------------
# gen --soc comprueba que los pseudos sean totalmente relativistas
# ----------------------------------------------------------------------
def _upf_falso(carpeta, nombre, element="Si", rel="scalar"):
    """Un UPF v2 mínimo con la cabecera que mira el selector de pseudos."""
    cab = (f'element="{element}" pseudo_type="NC" relativistic="{rel}" '
           f'functional="PBE" z_valence="4.000000e+00" '
           f'wfc_cutoff="4.000000e+01" rho_cutoff="1.600000e+02" '
           f'mesh_size="1000"')
    (carpeta / nombre).write_text(
        f'<UPF version="2.0.1">\n  <PP_HEADER {cab}/>\n</UPF>\n')


def _generar_soc(tmp_path, rel):
    from ase.build import bulk
    from qekit.modules import inputgen
    pdir = tmp_path / "pseudos"
    pdir.mkdir()
    _upf_falso(pdir, "Si.pbe.UPF", rel=rel)
    opts = inputgen.GenOptions(preset="scf", outdir=str(tmp_path / "calc"),
                               pseudo_dir=str(pdir), kspacing=0.55, soc=True)
    return inputgen.generate(bulk("Si", "diamond", 5.4073), opts)


def test_gen_soc_se_niega_con_pseudos_escalares(tmp_path):
    """`gen --soc` escribía lspinorb sobre pseudos escalares sin decir nada:
    pw.x no falla, da un desdoblamiento espín-órbita de cero que parece un
    resultado. `sweep.check_soc_pseudos` existía pero nadie lo llamaba.
    """
    from qekit.core.errors import ErrorDeUso
    with pytest.raises(ErrorDeUso, match="TOTALMENTE RELATIVISTAS"):
        _generar_soc(tmp_path, rel="scalar")
    assert not (tmp_path / "calc" / "scf.in").exists()


def test_gen_soc_pasa_con_pseudos_totalmente_relativistas(tmp_path):
    rep = _generar_soc(tmp_path, rel="full")
    assert "lspinorb" in rep or "Espín-órbita" in rep
    assert (tmp_path / "calc").exists()


# ----------------------------------------------------------------------
# La malla que se escribe es uniforme centrada en Γ, no Monkhorst-Pack
# ----------------------------------------------------------------------
def test_la_malla_se_describe_como_lo_que_es():
    """`_kgrid_card` escribe `n1 n2 n3 0 0 0`: centrada en Γ, sin
    desplazamiento. Con n par una malla Monkhorst-Pack no contiene Γ, así
    que llamarla así en el texto de métodos era una afirmación falsa.
    """
    from qekit.core import kpoints
    from qekit.modules import inputgen
    assert inputgen._kgrid_card((4, 4, 4)).strip().endswith("4 4 4 0 0 0")
    assert "Malla Monkhorst-Pack automática" not in kpoints.__doc__
    assert "centrada en Γ" in kpoints.__doc__
    # y el texto de métodos de la ficha dice lo mismo
    import inspect
    from qekit.modules import datasheet
    fuente = inspect.getsource(datasheet)
    assert "malla de \"\n        f\"Monkhorst-Pack" not in fuente
    assert "malla uniforme centrada" in fuente


# ----------------------------------------------------------------------
# hubbard.prepare con la sintaxis de QE >= 7.1
# ----------------------------------------------------------------------
def _nio():
    from ase.build import bulk
    return bulk("NiO", "rocksalt", a=4.17)


def test_hubbard_prepare_legacy_sigue_escribiendo_lda_plus_u(tmp_path):
    from qekit.modules import hubbard
    hubbard.prepare(_nio(), outdir=str(tmp_path), pseudo_dir="/no/existe",
                    kspacing=0.6)
    scf = (tmp_path / "scf.in").read_text()
    assert "lda_plus_u" in scf and "Hubbard_U(" in scf
    assert "U_projection_type = 'ortho-atomic'" in scf
    assert "\nHUBBARD" not in scf


def test_hubbard_prepare_card_escribe_la_tarjeta_hubbard(tmp_path):
    """`prepare` escribía siempre lda_plus_u, que pw.x >= 7.1 rechaza,
    aunque `gen --hubbard-style card` ya sabía escribir la tarjeta.
    """
    from qekit.modules import hubbard
    hubbard.prepare(_nio(), outdir=str(tmp_path), pseudo_dir="/no/existe",
                    kspacing=0.6, hubbard_style="card",
                    proyeccion="atomic")
    scf = (tmp_path / "scf.in").read_text()
    assert "lda_plus_u" not in scf and "Hubbard_U(" not in scf
    assert "U_projection_type" not in scf     # obsoleto en QE >= 7.1
    assert "HUBBARD (atomic)" in scf         # la proyección va en la tarjeta
    assert "U Ni-3d" in scf


def test_hubbard_prepare_rechaza_una_sintaxis_inventada(tmp_path):
    from qekit.core.errors import ErrorDeUso
    from qekit.modules import hubbard
    with pytest.raises(ErrorDeUso, match="hubbard_style"):
        hubbard.prepare(_nio(), outdir=str(tmp_path),
                        pseudo_dir="/no/existe", hubbard_style="xml")


def test_el_parser_de_hubbard_acepta_hubbard_style():
    from qekit import cli
    args = cli.build_parser().parse_args(
        ["hubbard", "NiO.cif", "--hubbard-style", "card"])
    assert args.hubbard_style == "card"
    args = cli.build_parser().parse_args(["hubbard", "NiO.cif"])
    assert args.hubbard_style == "legacy"


# ----------------------------------------------------------------------
# wannier: la cabecera de la DOS interpolada dice su normalización
# ----------------------------------------------------------------------
def test_la_cabecera_de_la_dos_de_wannier_declara_que_no_lleva_espin():
    """La cabecera decía "estados/eV/celda" a secas, pero la curva integra a
    num_wann (sin el 2 de espín): comparada con dos.x salía la mitad.
    """
    from qekit.modules import wannier
    assert "sin factor de espín" in wannier.DOS_UNIDADES
    assert "num_wann" in wannier.DOS_UNIDADES
    assert "num_wann" in wannier.dos_interpolada.__doc__
