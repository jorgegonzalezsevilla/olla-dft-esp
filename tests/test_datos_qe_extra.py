"""Pruebas de los módulos que leen datos de cálculos ya hechos.

Diagnóstico, auditoría, base de datos y casco convexo. Casi todo se prueba
con datos sintéticos de respuesta conocida, que es lo único que permite
verificar de verdad un casco convexo o un clasificador de convergencia.
"""

from pathlib import Path

import pytest


# ----------------------------------------------------------------------
# Casco convexo: respuestas exactas calculables a mano
# ----------------------------------------------------------------------
def test_casco_binario_contra_calculo_a_mano():
    """A, AB y B en el casco; A3B y AB3 exactamente 0.05 eV/át por encima.

    La recta A-AB vale -0.15 en x=0.25, y A3B está en -0.10: la diferencia
    es 0.05. Igual por simetría del otro lado.
    """
    from qekit.modules import thermo
    filas = [("A", {"A": 1}, 0.0), ("B", {"B": 1}, 0.0),
             ("AB", {"A": 1, "B": 1}, -0.60),
             ("A3B", {"A": 3, "B": 1}, -0.40),
             ("AB3", {"A": 1, "B": 3}, -0.40)]
    r = thermo.from_table(filas)
    d = {f.nombre: f for f in r.fases}
    assert d["AB"].e_form == pytest.approx(-0.30)
    assert d["A3B"].e_form == pytest.approx(-0.10)
    assert d["AB"].en_casco and d["A"].en_casco and d["B"].en_casco
    assert not d["A3B"].en_casco and not d["AB3"].en_casco
    assert d["A3B"].e_hull == pytest.approx(0.05, abs=1e-9)
    assert d["AB3"].e_hull == pytest.approx(0.05, abs=1e-9)


def test_una_fase_mas_estable_redefine_el_casco():
    """Si A3B baja lo suficiente, entra al casco y empuja a AB fuera de
    la recta directa A-B: el casco tiene que recalcularse, no fijarse."""
    from qekit.modules import thermo
    filas = [("A", {"A": 1}, 0.0), ("B", {"B": 1}, 0.0),
             ("AB", {"A": 1, "B": 1}, -0.60),
             ("A3B", {"A": 3, "B": 1}, -1.60)]      # E_f = -0.40
    r = thermo.from_table(filas)
    d = {f.nombre: f for f in r.fases}
    assert d["A3B"].e_form == pytest.approx(-0.40)
    assert d["A3B"].en_casco
    # la recta A3B(-0.40 en 0.25) -> B(0 en 1) vale -0.40*(1-0.5)/(1-0.25)
    # = -0.2667 en x=0.5; AB está en -0.30, o sea por debajo => estable
    assert d["AB"].en_casco
    assert d["AB"].e_hull == pytest.approx(0.0, abs=1e-9)


def test_casco_marca_inestable_lo_que_esta_por_encima():
    from qekit.modules import thermo
    filas = [("A", {"A": 1}, 0.0), ("B", {"B": 1}, 0.0),
             ("AB", {"A": 1, "B": 1}, +0.40)]        # E_f = +0.20
    r = thermo.from_table(filas)
    d = {f.nombre: f for f in r.fases}
    assert d["AB"].e_form == pytest.approx(0.20)
    assert not d["AB"].en_casco
    assert d["AB"].e_hull == pytest.approx(0.20, abs=1e-9)


def test_casco_ternario():
    """Con tres elementos entra el camino de scipy. Un compuesto ternario
    muy estable debe quedar en el casco y uno flojo por encima."""
    from qekit.modules import thermo
    filas = [("A", {"A": 1}, 0.0), ("B", {"B": 1}, 0.0),
             ("C", {"C": 1}, 0.0),
             ("ABC", {"A": 1, "B": 1, "C": 1}, -3.0),      # -1.0 eV/át
             ("A2BC", {"A": 2, "B": 1, "C": 1}, -0.4)]     # -0.1 eV/át
    r = thermo.from_table(filas)
    d = {f.nombre: f for f in r.fases}
    assert d["ABC"].e_form == pytest.approx(-1.0)
    assert d["ABC"].en_casco
    assert d["A2BC"].e_hull is not None and d["A2BC"].e_hull > 0.1


def test_casco_exige_las_referencias_elementales():
    from qekit.modules import thermo
    r = thermo.from_table([("AB", {"A": 1, "B": 1}, -1.0)])
    assert set(r.faltan_ref) == {"A", "B"}
    assert "referencias elementales" in thermo.report(r)


# ----------------------------------------------------------------------
# Clasificador de convergencia SCF
# ----------------------------------------------------------------------
def _stdout_falso(tmp_path, accuracies, beta=0.4, convergio=False):
    txt = ["     Program PWSCF v.6.6 starts"]
    for i, a in enumerate(accuracies, start=1):
        txt += [f"     iteration # {i:3d}     ecut=    60.00 Ry     "
                f"beta= {beta:.2f}",
                "     total energy              =     -15.8000000 Ry",
                f"     estimated scf accuracy    <    {a:.8f} Ry"]
    txt.append("     convergence has been achieved" if convergio
               else "     convergence NOT achieved after 20 iterations")
    f = tmp_path / "scf.out"
    f.write_text("\n".join(txt))
    return f


def test_detecta_oscilacion_de_carga(tmp_path):
    """El error sube y baja: hay que mezclar MENOS, no más."""
    from qekit.modules import diagnose
    acc = [1e-1, 5e-2, 2e-1, 3e-2, 1e-1, 2e-2, 8e-2, 1e-2]
    h = diagnose.read_scf_history(_stdout_falso(tmp_path, acc, beta=0.7))
    assert h.patologia == "oscilacion"
    assert "mixing_mode = 'local-TF'" in h.consejo
    assert "lo EMPEORA" in h.consejo


def test_detecta_convergencia_lenta(tmp_path):
    """El error baja siempre pero despacio: aquí sí hay que mezclar MÁS."""
    from qekit.modules import diagnose
    acc = [1e-1 * (0.45 ** i) for i in range(10)]
    h = diagnose.read_scf_history(_stdout_falso(tmp_path, acc, beta=0.2))
    assert h.patologia == "lenta"
    assert "mixing_beta" in h.consejo
    assert "electron_maxstep" in h.consejo


def test_detecta_estancamiento(tmp_path):
    from qekit.modules import diagnose
    acc = [1e-2, 9e-3, 8.5e-3, 8.3e-3, 8.2e-3, 8.2e-3, 8.1e-3, 8.1e-3]
    h = diagnose.read_scf_history(_stdout_falso(tmp_path, acc))
    assert h.patologia == "estancada"


def test_no_diagnostica_nada_si_convergio(tmp_path):
    from qekit.modules import diagnose
    acc = [1e-1, 1e-3, 1e-6, 1e-10]
    h = diagnose.read_scf_history(_stdout_falso(tmp_path, acc,
                                                convergio=True))
    assert h.converged is True
    assert h.patologia == ""


# ----------------------------------------------------------------------
# Auditoría
# ----------------------------------------------------------------------
class _FalsoResultado:
    def __init__(self, **kw):
        self.functional = kw.get("functional", "PBE")
        self.pseudo_files = kw.get("pseudo_files", {"Si": "Si.upf"})
        self.ecutwfc = kw.get("ecutwfc", 60.0)
        self.ecutrho = kw.get("ecutrho", 480.0)
        self.smearing = kw.get("smearing", "")
        self.degauss = kw.get("degauss", None)
        self.occupations_kind = kw.get("occupations_kind", "fixed")
        self.nspin = kw.get("nspin", 1)
        self.converged = kw.get("converged", True)
        self.calculation = kw.get("calculation", "scf")
        self.kgrid = kw.get("kgrid", (8, 8, 8))
        self.volume = kw.get("volume", 40.0)
        self.symbols = kw.get("symbols", ["Si", "Si"])
        self.total_energy = kw.get("total_energy", -215.0)
        self.homo = self.lumo = None
        self.pressure = self.max_force = None
        self.total_magnetization = None
        self.n_scf_steps = 10
        self.wall_time = 1.0
        self.xml_path = "x.xml"

    @property
    def fingerprint(self):
        return (self.functional, tuple(sorted(self.pseudo_files.items())),
                self.ecutwfc, self.ecutrho, self.smearing, self.degauss,
                self.occupations_kind, self.nspin)


def _run(path, **kw):
    from qekit.modules.audit import RunInfo
    return RunInfo(path=path, result=_FalsoResultado(**kw))


def test_auditoria_acepta_un_conjunto_homogeneo():
    from qekit.modules import audit
    a = audit.audit([_run("a"), _run("b"), _run("c")])
    assert a["comparables"]
    assert "COMPARABLES" in audit.report(a)


def test_auditoria_nombra_el_parametro_que_difiere():
    from qekit.modules import audit
    a = audit.audit([_run("a"), _run("b", ecutwfc=80.0)])
    assert not a["comparables"]
    claves = [c for c, _ in a["difieren"]]
    assert "ecutwfc" in claves
    rep = audit.report(a)
    assert "NO COMPARABLES" in rep and "ecutwfc" in rep


def test_auditoria_detecta_pseudos_distintos():
    """El caso más traicionero: mismo elemento, otro pseudopotencial."""
    from qekit.modules import audit
    a = audit.audit([_run("a"),
                     _run("b", pseudo_files={"Si": "Si_OTRO.upf"})])
    assert not a["comparables"]
    assert "pseudos" in [c for c, _ in a["difieren"]]


def test_auditoria_no_marca_un_nscf_como_no_convergido():
    """Un nscf no tiene ciclo SCF: exigirle convergencia es falso positivo."""
    from qekit.modules import audit
    a = audit.audit([_run("bandas", calculation="nscf", converged=False)])
    assert not a["no_convergidos"]
    assert len(a["sin_energia"]) == 1


def test_auditoria_si_marca_un_scf_no_convergido():
    from qekit.modules import audit
    a = audit.audit([_run("malo", converged=False)])
    assert len(a["no_convergidos"]) == 1
    assert "NO CONVERGIERON" in audit.report(a)


def test_densidad_de_k_es_comparable_entre_celdas():
    """Dos celdas distintas con el MISMO espaciado de k dan la misma
    densidad.

    Si la celda crece 8x en volumen (2x lineal), los vectores recíprocos se
    encogen a la mitad, así que el mismo espaciado se consigue con la MITAD
    de puntos por dirección — no con el doble, que es el error intuitivo.
    La densidad puntos/Å⁻³ captura eso; el número crudo de puntos no.
    """
    from qekit.modules import audit
    r1 = _FalsoResultado(kgrid=(8, 8, 8), volume=40.0)
    r2 = _FalsoResultado(kgrid=(4, 4, 4), volume=320.0)
    assert audit.kdensity(r1) == pytest.approx(audit.kdensity(r2), rel=1e-9)
    # y una malla claramente más pobre da una densidad menor
    r3 = _FalsoResultado(kgrid=(4, 4, 4), volume=40.0)
    assert audit.kdensity(r3) < audit.kdensity(r1)


# ----------------------------------------------------------------------
# Base de datos
# ----------------------------------------------------------------------
def test_base_de_datos_indexa_y_actualiza(tmp_path):
    from qekit.modules import audit
    db = tmp_path / "q.db"
    nuevos, act = audit.index([_run("/x/a"), _run("/x/b")], db)
    assert (nuevos, act) == (2, 0)
    nuevos, act = audit.index([_run("/x/a")], db)      # misma ruta
    assert (nuevos, act) == (0, 1)
    filas = audit.query("SELECT * FROM calculos", db)
    assert len(filas) == 2
    assert filas[0]["energia_por_atomo_eV"] == pytest.approx(-107.5)


def test_base_de_datos_solo_admite_select(tmp_path):
    from qekit.modules import audit
    db = tmp_path / "q.db"
    audit.index([_run("/x/a")], db)
    with pytest.raises(ValueError, match="solo se admiten consultas SELECT"):
        audit.query("DELETE FROM calculos", db)


def test_base_de_datos_exporta_json(tmp_path):
    import json

    from qekit.modules import audit
    db = tmp_path / "q.db"
    audit.index([_run("/x/a")], db)
    out = audit.export_json(db, tmp_path / "c.json")
    doc = json.loads(Path(out).read_text())
    assert doc["n"] == 1 and doc["calculos"][0]["formula"] == "Si2"
