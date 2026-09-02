"""Selector de pseudopotenciales, TDDFPT y transporte balístico.

Lo que se puede comprobar sin correr Quantum ESPRESSO: la lógica de
selección (que es donde estaban los errores silenciosos), el formato de
los inputs, y la lectura y el análisis de resultados ya calculados.
"""

import gzip
from pathlib import Path

import numpy as np
import pytest

from qekit.core import pseudo as ps
from qekit.core.errors import ErrorDeUso
from qekit.modules import ballistic as bl, pseudos as pz, tddft as td

DATOS = Path(__file__).parent / "datos"
PSEUDOS_SISTEMA = Path("/usr/share/espresso/pseudo")


def _hay_pseudos():
    return PSEUDOS_SISTEMA.is_dir() and \
        len(list(PSEUDOS_SISTEMA.glob("*.UPF"))) > 5


# ======================================================================
# Selector de pseudopotenciales
# ======================================================================
def _upf(tmp_path, nombre, **campos):
    """Un UPF v2 mínimo con los atributos que el selector mira."""
    at = {"element": campos.get("element", "Si"),
          "pseudo_type": campos.get("tipo", "NC"),
          "relativistic": campos.get("rel", "scalar"),
          "functional": campos.get("func", "PBE"),
          "z_valence": f"{campos.get('zval', 4.0):.6e}",
          "wfc_cutoff": f"{campos.get('ecut', 40.0):.6e}",
          "rho_cutoff": f"{campos.get('rho', 160.0):.6e}",
          "mesh_size": "1000"}
    cab = " ".join(f'{k}="{v}"' for k, v in at.items())
    cuerpo = f'<UPF version="2.0.1">\n  <PP_HEADER {cab}/>\n'
    if campos.get("gipaw"):
        cuerpo += "  <PP_GIPAW><PP_GIPAW_CORE_ORBITALS/></PP_GIPAW>\n"
    cuerpo += "</UPF>\n"
    f = tmp_path / nombre
    f.write_text(cuerpo)
    return f


def test_lee_lo_que_el_selector_necesita(tmp_path):
    f = _upf(tmp_path, "Si.a.UPF", tipo="US", rel="full", func="PBESOL",
             zval=12.0, ecut=55.0, gipaw=True)
    c = pz.leer(f)
    assert c.elemento == "Si" and c.tipo == "US"
    assert c.relativista == "full" and c.funcional == "PBESOL"
    assert c.z_valence == 12.0 and c.ecutwfc == 55.0
    assert c.gipaw is True


def test_optica_descarta_los_que_no_conservan_la_norma(tmp_path):
    cands = [pz.leer(_upf(tmp_path, "Si.us.UPF", tipo="US")),
             pz.leer(_upf(tmp_path, "Si.nc.UPF", tipo="NC")),
             pz.leer(_upf(tmp_path, "Si.paw.UPF", tipo="PAW"))]
    ev = pz.evaluar(cands, "optics")
    buenos = [c for c in ev if c.ok]
    assert [c.tipo for c in buenos] == ["NC"]
    assert all("epsilon.x" in c.descartado for c in ev if not c.ok)


def test_soc_exige_totalmente_relativista_en_elementos_pesados(tmp_path):
    cands = [pz.leer(_upf(tmp_path, "Au.s.UPF", element="Au", rel="scalar")),
             pz.leer(_upf(tmp_path, "Au.f.UPF", element="Au", rel="full"))]
    ev = pz.evaluar(cands, "soc")
    buenos = [c for c in ev if c.ok]
    assert len(buenos) == 1 and buenos[0].relativista == "full"


def test_soc_es_indulgente_con_los_elementos_ligeros(tmp_path):
    """En el oxigeno no hay pseudos relativistas y no hacen falta."""
    cands = [pz.leer(_upf(tmp_path, "O.s.UPF", element="O", rel="scalar"))]
    ev = pz.evaluar(cands, "soc")
    assert ev[0].ok
    assert any("despreciable" in n for n in ev[0].notas)


def test_xanes_exige_gipaw(tmp_path):
    cands = [pz.leer(_upf(tmp_path, "Si.sin.UPF")),
             pz.leer(_upf(tmp_path, "Si.con.UPF", gipaw=True))]
    ev = pz.evaluar(cands, "xanes")
    buenos = [c for c in ev if c.ok]
    assert len(buenos) == 1 and buenos[0].gipaw


def test_exigir_funcional_descarta_los_demas(tmp_path):
    cands = [pz.leer(_upf(tmp_path, "Si.pbe.UPF", func="PBE")),
             pz.leer(_upf(tmp_path, "Si.pz.UPF", func="PZ"))]
    ev = pz.evaluar(cands, "general", funcional="PBE")
    buenos = [c for c in ev if c.ok]
    assert len(buenos) == 1 and buenos[0].funcional == "PBE"


@pytest.mark.parametrize("a,b", [
    ("PBE", "SLA PW PBX PBC"),
    ("PZ", "SLA PZ NOGX NOGC"),
    ("pbe", "PBE"),
])
def test_el_mismo_funcional_escrito_de_dos_formas(a, b):
    """QE escribe unas veces el nombre corto y otras las cuatro piezas."""
    assert pz._mismo_funcional(a, b)


def test_funcionales_distintos_no_se_confunden():
    assert not pz._mismo_funcional("PBE", "PZ")
    assert not pz._mismo_funcional("PBE", "PBESOL")


def test_coherencia_detecta_funcionales_mezclados(tmp_path):
    a = pz.leer(_upf(tmp_path, "Ni.UPF", element="Ni", func="PBE"))
    b = pz.leer(_upf(tmp_path, "O.UPF", element="O", func="BLYP"))
    avisos = pz.coherencia({"Ni": a, "O": b})
    assert any("FUNCIONALES DISTINTOS" in x for x in avisos)


def test_coherencia_avisa_al_mezclar_NC_con_ultrasuave(tmp_path):
    a = pz.leer(_upf(tmp_path, "A.UPF", element="Si", tipo="NC"))
    b = pz.leer(_upf(tmp_path, "B.UPF", element="O", tipo="US"))
    avisos = pz.coherencia({"Si": a, "O": b})
    assert any("ecutrho" in x for x in avisos)


def test_coherencia_avisa_de_cutoffs_muy_dispares(tmp_path):
    a = pz.leer(_upf(tmp_path, "A.UPF", element="Si", ecut=30.0))
    b = pz.leer(_upf(tmp_path, "B.UPF", element="O", ecut=100.0))
    avisos = pz.coherencia({"Si": a, "O": b})
    assert any("decide el coste" in x for x in avisos)


def test_ninguno_sirve_lo_dice_con_el_motivo(tmp_path):
    _upf(tmp_path, "Si.us.UPF", tipo="US")
    with pytest.raises(ErrorDeUso, match="ninguno sirve"):
        pz.elegir("Si", str(tmp_path), tarea="optics")


def test_tarea_desconocida(tmp_path):
    with pytest.raises(ErrorDeUso, match="Opciones"):
        pz.evaluar([], "noexiste")


def test_carpeta_inexistente_dice_como_arreglarlo():
    with pytest.raises(ErrorDeUso, match="config set pseudo_dir"):
        pz.candidatos("Si", "/no/existe/en/ningun/sitio")


# --- integracion con resolve ------------------------------------------
def test_forzar_un_pseudo_manda(tmp_path):
    _upf(tmp_path, "Si.a.UPF")
    _upf(tmp_path, "Si.b.UPF")
    r = ps.resolve(["Si"], str(tmp_path), forzados={"Si": "Si.b.UPF"})
    assert r["Si"]["filename"] == "Si.b.UPF"


def test_forzar_algo_que_no_esta_avisa(tmp_path):
    _upf(tmp_path, "Si.a.UPF")
    with pytest.raises(ErrorDeUso, match="no esta en la carpeta"):
        ps.resolve(["Si"], str(tmp_path), forzados={"Si": "noexiste.UPF"})


def test_resolve_arregla_la_mezcla_de_funcionales(tmp_path):
    """El error silencioso que motivo el modulo: por orden alfabetico
    salia un Ni de PBE junto a un O de BLYP."""
    _upf(tmp_path, "Ni.aaa-blyp.UPF", element="Ni", func="BLYP")
    _upf(tmp_path, "Ni.zzz-pbe.UPF", element="Ni", func="PBE")
    _upf(tmp_path, "O.aaa-pbe.UPF", element="O", func="PBE")
    r = ps.resolve(["Ni", "O"], str(tmp_path))
    elegidos = {k: pz.leer(Path(tmp_path) / v["filename"]).funcional
                for k, v in r.items()}
    assert len(set(elegidos.values())) == 1, elegidos


def test_lo_forzado_no_se_toca_al_arreglar_funcionales(tmp_path):
    _upf(tmp_path, "Ni.blyp.UPF", element="Ni", func="BLYP")
    _upf(tmp_path, "Ni.pbe.UPF", element="Ni", func="PBE")
    _upf(tmp_path, "O.pbe.UPF", element="O", func="PBE")
    r = ps.resolve(["Ni", "O"], str(tmp_path),
                   forzados={"Ni": "Ni.blyp.UPF"})
    assert r["Ni"]["filename"] == "Ni.blyp.UPF"


@pytest.mark.skipif(not _hay_pseudos(), reason="sin tabla de pseudos")
def test_pseudos_reales_del_sistema():
    r = ps.resolve(["Ni", "O"], str(PSEUDOS_SISTEMA))
    funcs = {pz.leer(PSEUDOS_SISTEMA / v["filename"]).funcional
             for v in r.values()}
    assert len(funcs) == 1
    r = ps.resolve(["Ni"], str(PSEUDOS_SISTEMA), tarea="soc")
    assert "rel" in r["Ni"]["filename"]


@pytest.mark.skipif(not _hay_pseudos(), reason="sin tabla de pseudos")
def test_relativista_se_reconoce_en_upf_v1():
    """'Fully-Relativistic' con guion: reconocerlo mal deja el filtro
    de espin-orbita sin efecto."""
    assert ps.relativistic(
        PSEUDOS_SISTEMA / "Ni.rel-pbe-nd-rrkjus.UPF") == "full"


# ======================================================================
# TDDFPT
# ======================================================================
def test_input_de_lanczos_tiene_lo_imprescindible():
    txt = td.build_lanczos_input("Si", itermax=800, ipol=4)
    assert "itermax = 800" in txt
    assert "ipol = 4" in txt and "n_ipol = 3" in txt
    assert "&lr_input" in txt and "&lr_control" in txt


def test_polarizacion_simple_usa_una_sola_cadena():
    txt = td.build_lanczos_input("Si", ipol=2)
    assert "n_ipol = 1" in txt


def test_tamm_dancoff_y_rpa_se_declaran():
    assert "ltammd = .true." in td.build_lanczos_input("Si", ltammd=True)
    assert "lrpa = .true." in td.build_lanczos_input("Si", lrpa=True)


def test_el_espectro_pide_las_energias_en_eV():
    txt = td.build_spectrum_input("Si", emin=1.0, emax=10.0)
    assert "units = 1" in txt
    assert "start = 1.0" in txt and "end = 10.0" in txt


def test_extrapolacion_desconocida():
    with pytest.raises(ErrorDeUso, match="Opciones"):
        td.build_spectrum_input("Si", extrapolation="magia")


def test_metodo_desconocido():
    from ase.build import bulk
    with pytest.raises(ErrorDeUso, match="lanczos"):
        td.prepare(bulk("Si", "diamond", a=5.43), metodo="montecarlo")


def test_davidson_convierte_las_energias_a_Ry():
    txt = td.build_davidson_input("Si", emin=0.0, emax=13.605693)
    assert "finish = 1.000000" in txt


def _run_sintetico(gap=2.0, e_exciton=1.6):
    """Un espectro con un pico dentro del gap: la firma de un excitón."""
    e = np.linspace(0.0, 8.0, 800)
    s = (0.6 * np.exp(-0.5 * ((e - e_exciton) / 0.08) ** 2)
         + 1.0 * np.exp(-0.5 * ((e - (gap + 0.8)) / 0.5) ** 2))
    r = td.TddftRun(energias=e, total=s, gap_ip=gap)
    r.picos = td._picos(r)
    td._avisar(r)
    return r


def test_detecta_el_exciton_bajo_el_gap():
    r = _run_sintetico(gap=2.0, e_exciton=1.6)
    assert r.onset < 2.0
    assert any("excitón" in a or "exciton" in a for a in r.avisos)


def test_sin_exciton_lo_dice_tambien():
    r = _run_sintetico(gap=2.0, e_exciton=2.0)
    assert any("adiabático" in a or "adiabatico" in a for a in r.avisos)


def test_los_picos_salen_ordenados_por_altura():
    r = _run_sintetico()
    alturas = [h for _, h in r.picos]
    assert alturas == sorted(alturas, reverse=True)
    assert r.picos[0][1] == pytest.approx(1.0, abs=1e-6)


def test_collect_sin_datos_avisa(tmp_path):
    with pytest.raises(ErrorDeUso, match="turbo_spectrum"):
        td.collect(tmp_path)


def test_el_scf_de_tddfpt_apaga_la_simetria(tmp_path):
    """Sin nosym/noinv, turbo_*.x se planta DESPUES del scf."""
    from ase.build import molecule
    import numpy as _np
    m = molecule("H2")
    m.set_cell(_np.eye(3) * 12.0); m.center(); m.pbc = True
    td.prepare(m, outdir=str(tmp_path), pseudo_dir=str(tmp_path))
    txt = (tmp_path / "scf.in").read_text()
    assert "nosym" in txt and ".true." in txt
    assert "noinv" in txt


# ======================================================================
# Transporte balístico
# ======================================================================
def test_el_cuanto_de_conductancia():
    """G0 = 2e^2/h = 7.748e-5 S; un canal perfecto son 12.906 kOhm."""
    assert bl.G0 == pytest.approx(7.748e-5, rel=1e-3)
    assert bl.R0 == pytest.approx(12906.4, rel=1e-3)


def test_geometria_con_z_torcido_se_rechaza():
    from ase import Atoms
    a = Atoms("Al", positions=[[0, 0, 0]],
              cell=[[10, 0, 0], [0, 10, 0], [1, 0, 5]], pbc=True)
    problemas = bl.comprobar_geometria(a)
    assert any("z" in p for p in problemas)


def test_geometria_correcta_pasa():
    from ase import Atoms
    a = Atoms("Al", positions=[[0, 0, 0]],
              cell=[[10, 0, 0], [0, 10, 0], [0, 0, 4.5]], pbc=True)
    assert bl.comprobar_geometria(a) == []


def test_celdas_del_plano_distintas_se_rechazan(tmp_path):
    from ase import Atoms
    e = Atoms("Al", positions=[[0, 0, 0]],
              cell=[[10, 0, 0], [0, 10, 0], [0, 0, 4.5]], pbc=True)
    d = Atoms("Al", positions=[[0, 0, 0]],
              cell=[[12, 0, 0], [0, 12, 0], [0, 0, 9.0]], pbc=True)
    with pytest.raises(ErrorDeUso, match="misma celda en el plano"):
        bl.prepare(e, outdir=str(tmp_path), dispersor=d,
                   pseudo_dir=str(tmp_path))


def test_limites_z_en_unidades_de_alat():
    from ase import Atoms
    a = Atoms("Al2", positions=[[0, 0, 1.0], [0, 0, 3.0]],
              cell=[[10, 0, 0], [0, 10, 0], [0, 0, 5.0]], pbc=True)
    zmin, zmax = bl.limites_z(a)
    assert zmin == pytest.approx(0.1)
    assert zmax == pytest.approx(0.3)


def test_input_de_pwcond_sin_dispersor_no_pide_transmision():
    txt = bl.build_cond_input("elec", ikind=0)
    assert "ikind=0" in txt
    assert "tran_file" not in txt
    assert txt.strip().split("\n")[-1].strip().isdigit()


def test_input_de_pwcond_con_dispersor_pide_transmision():
    txt = bl.build_cond_input("elec", ikind=1, prefixs="disp")
    assert "prefixs='disp'" in txt and "tran_file" in txt


def test_conductancia_cuantizada_se_reconoce():
    e = np.linspace(-1, 1, 21)
    r = bl.CondRun(energias=e, transmision=np.full_like(e, 2.0), ikind=1)
    r.G_fermi = 2.0
    assert r.G_siemens == pytest.approx(2 * bl.G0)
    assert r.R_ohm == pytest.approx(bl.R0 / 2)
    assert "cuantización" in bl.report(r)


def test_transmision_mayor_que_los_canales_es_imposible():
    e = np.linspace(-1, 1, 5)
    r = bl.CondRun(energias=e, transmision=np.full_like(e, 2.5),
                   canales=np.ones_like(e), ikind=1)
    bl._avisar(r)
    assert any("imposible" in a for a in r.avisos)


def test_transmision_negativa_avisa():
    e = np.linspace(-1, 1, 5)
    r = bl.CondRun(energias=e, transmision=np.full_like(e, -0.5), ikind=1)
    bl._avisar(r)
    assert any("NEGATIVAS" in a for a in r.avisos)


def test_lee_el_hilo_de_aluminio_calculado():
    """Hilo monoatomico de Al: 1 canal de la banda s, 3 donde entran
    las p degeneradas."""
    d = DATOS / "balistico_al"
    if not (d / "cond.out.gz").exists():
        pytest.skip("falta el calculo de prueba")
    import tempfile
    with tempfile.TemporaryDirectory() as t:
        t = Path(t)
        (t / "cond.out").write_bytes(
            gzip.decompress((d / "cond.out.gz").read_bytes()))
        (t / "cond.in").write_bytes((d / "cond.in").read_bytes())
        run = bl.collect(t)
    assert run.ikind == 0
    assert run.canales is not None
    assert int(run.canales.max()) == 3
    assert int(run.canales.min()) == 0
    assert any("NO es la conductancia" in a for a in run.avisos)


def test_collect_sin_nada_avisa(tmp_path):
    with pytest.raises(ErrorDeUso, match="pwcond"):
        bl.collect(tmp_path)


def test_lee_las_excitaciones_del_etileno():
    """Etileno con TDDFPT: la pi->pi* es la primera brillante y va
    polarizada a lo largo de un solo eje."""
    d = DATOS / "tddft_c2h4"
    if not (d / "CH2.eigen").exists():
        pytest.skip("falta el calculo de prueba")
    run = td.collect(d, metodo="davidson")
    assert len(run.excitaciones) == 6
    e0, f0 = run.excitaciones[0]
    assert 6.0 < e0 < 7.0                  # LDA subestima la experimental
    assert f0 > 0.02                       # es brillante
    # las energias suben
    assert all(a[0] <= b[0] for a, b in
               zip(run.excitaciones, run.excitaciones[1:]))
    # la mayoria son oscuras
    assert sum(1 for _, f in run.excitaciones if f < 1e-3) >= 3


def test_las_polarizaciones_no_ensucian_los_espectros():
    """El fallo real: los tripletes por estado acabaron en `componentes`,
    junto a arrays de miles de puntos, y la exportacion reventaba."""
    d = DATOS / "tddft_c2h4"
    if not (d / "CH2.eigen").exists():
        pytest.skip("falta el calculo de prueba")
    run = td.collect(d, metodo="davidson")
    assert len(run.polarizaciones) == len(run.excitaciones)
    for v in run.componentes.values():
        assert len(v) == len(run.energias)


def test_exportar_davidson_no_revienta(tmp_path):
    d = DATOS / "tddft_c2h4"
    if not (d / "CH2.eigen").exists():
        pytest.skip("falta el calculo de prueba")
    run = td.collect(d, metodo="davidson")
    salidas = td.export(run, str(tmp_path))
    assert any("EXCITACIONES" in s for s in salidas)
    a = np.loadtxt(tmp_path / "TDDFT_EXCITACIONES.dat")
    assert a.shape == (6, 2)


def test_el_eigen_esta_en_rydberg_y_se_convierte():
    """0.4774 Ry = 6.50 eV; leerlo como eV daria 0.48 eV, absurdo."""
    d = DATOS / "tddft_c2h4"
    if not (d / "CH2.eigen").exists():
        pytest.skip("falta el calculo de prueba")
    crudo = np.loadtxt(d / "CH2.eigen", comments="#")
    run = td.collect(d, metodo="davidson")
    assert run.excitaciones[0][0] == pytest.approx(
        crudo[0, 0] * td.RY_EV, rel=1e-6)
