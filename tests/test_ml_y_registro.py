"""Pruebas del bloque de ML, procedencia, registro de incidencias y
recomendador.

Lo importante aquí no es que MACE funcione —eso depende de una dependencia
opcional— sino que QEkit se comporte bien SIN ella, que nunca mezcle
energías de un potencial aprendido con energías DFT, y que el registro de
incidencias capture lo necesario para reproducir un fallo.
"""

import json
import sys
from pathlib import Path

import pytest


# ----------------------------------------------------------------------
# El MLIP es una dependencia OPCIONAL
# ----------------------------------------------------------------------
def test_sin_torch_el_mensaje_explica_que_instalar(monkeypatch):
    """Si falta el paquete, no debe salir un ImportError pelado."""
    from qekit.modules import mlip
    import builtins
    real = builtins.__import__

    def sin_mace(name, *a, **kw):
        if name.startswith("mace"):
            raise ImportError("No module named 'mace'")
        return real(name, *a, **kw)

    monkeypatch.setattr(builtins, "__import__", sin_mace)
    with pytest.raises(ImportError, match="pip install torch"):
        mlip.calculator("mace")


def test_qekit_importa_sin_los_paquetes_de_ml():
    """El resto de QEkit no puede depender de torch ni de mace."""
    import subprocess
    codigo = (
        "import sys\n"
        "for m in list(sys.modules):\n"
        "    if m.split('.')[0] in ('torch','mace','chgnet','matgl'):\n"
        "        del sys.modules[m]\n"
        "import qekit.cli, qekit.modules.mlip\n"
        "assert 'torch' not in sys.modules, 'importar QEkit arrastra torch'\n"
        "print('ok')\n")
    r = subprocess.run([sys.executable, "-c", codigo], capture_output=True,
                       text=True)
    assert r.returncode == 0, r.stderr
    assert "ok" in r.stdout


def test_modelo_desconocido_se_rechaza():
    from qekit.modules import mlip
    with pytest.raises(ValueError, match="modelo desconocido"):
        mlip.calculator("inventado")


# ----------------------------------------------------------------------
# Procedencia: una energía MLIP no se mezcla con una DFT
# ----------------------------------------------------------------------
def test_marca_de_procedencia_se_escribe_y_se_lee(tmp_path):
    from qekit.modules import mlip
    run = mlip.MlipRun(modelo="mace", detalle="MACE-MP-0 (PBE)",
                       pasos=5, fmax_final=0.004, cambio_volumen=1.2)
    f = mlip.write_provenance(run, tmp_path / "relajado.cif")
    assert Path(f).exists()
    d = mlip.read_provenance(tmp_path)
    assert d["origen"] == "mlip:mace"
    assert "NO es comparable" in d["aviso"]


def test_audit_separa_energias_mlip_de_las_dft():
    """Es el agujero que abre el MLIP y que audit tiene que cerrar."""
    from qekit.modules import audit
    from tests.test_datos_qe_extra import _run
    a_dft = _run("/x/dft")
    a_ml = _run("/x/mlip")
    a_ml.origen = "mlip:mace"
    a = audit.audit([a_dft, a_ml])
    assert not a["comparables"]
    claves = [c for c, _ in a["difieren"]]
    assert "origen" in claves
    assert "potencial aprendido" in audit.report(a)


def test_casco_se_niega_a_mezclar_origenes():
    from qekit.modules import thermo
    from tests.test_datos_qe_extra import _run
    a, b = _run("/x/a"), _run("/x/b")
    b.origen = "mlip:mace"
    res = thermo.from_runs([a, b])
    assert res.fases == []
    assert any("potenciales aprendidos" in w for w in res.warnings)


def test_la_base_registra_el_origen(tmp_path):
    from qekit.modules import audit
    from tests.test_datos_qe_extra import _run
    x = _run("/x/ml")
    x.origen = "mlip:mace"
    db = tmp_path / "q.db"
    audit.index([x], db)
    fila = audit.query("SELECT origen, huella FROM calculos", db)[0]
    assert fila["origen"] == "mlip:mace"
    assert "mlip:mace" in fila["huella"]


# ----------------------------------------------------------------------
# Registro de incidencias
# ----------------------------------------------------------------------
def test_incidencia_captura_lo_necesario_para_reproducir(tmp_path):
    from qekit.core import provenance
    from qekit.modules import feedback as fb
    provenance.record_argv(["olla-dft", "optics", "x.cif", "--scissor", "0.65"])
    try:
        raise ValueError("algo se rompió")
    except ValueError as exc:
        inc = fb.registrar("no debería fallar", exc=exc, dir_=tmp_path)
    assert inc.tipo == "error"
    assert "--scissor 0.65" in inc.comando          # el comando exacto
    assert "ValueError" in inc.excepcion
    assert "algo se rompió" in inc.traceback        # la traza completa
    assert inc.versiones["qekit"] and inc.versiones["python"]
    assert "numpy" in inc.versiones                 # versiones de deps
    guardado = json.loads(
        (tmp_path / inc.id / "incidencia.json").read_text())
    assert guardado["id"] == inc.id


def test_incidencia_manual_sin_excepcion(tmp_path):
    from qekit.modules import feedback as fb
    inc = fb.registrar("la salida de optics no dice las unidades",
                       dir_=tmp_path)
    assert inc.tipo == "manual"
    assert not inc.traceback
    assert fb.listar(tmp_path)[0].id == inc.id


def test_adjuntar_copia_el_archivo(tmp_path):
    from qekit.modules import feedback as fb
    origen = tmp_path / "estructura.cif"
    origen.write_text("datos")
    inc = fb.registrar("con adjunto", adjuntos=[str(origen)],
                       dir_=tmp_path)
    assert inc.adjuntos == ["estructura.cif"]
    assert (tmp_path / inc.id / "estructura.cif").read_text() == "datos"


def test_cerrar_incidencia(tmp_path):
    from qekit.modules import feedback as fb
    inc = fb.registrar("algo", dir_=tmp_path)
    assert fb.cerrar(inc.id, nota="arreglado en 0.11", dir_=tmp_path)
    guardada = fb.listar(tmp_path)[0]
    assert guardada.estado == "cerrada" and "0.11" in guardada.nota
    assert not fb.cerrar("noexiste", dir_=tmp_path)


def test_estadisticas_senalan_el_subcomando_problematico(tmp_path):
    from qekit.core import provenance
    from qekit.modules import feedback as fb
    for cmd in (["qekit", "optics", "a"], ["qekit", "optics", "b"],
                ["qekit", "bands", "c"]):
        provenance.record_argv(cmd)
        fb.registrar("x", dir_=tmp_path)
    st = fb.estadisticas(tmp_path)
    assert st["total"] == 3
    assert list(st["por_comando"])[0] == "optics"    # el más frecuente
    assert st["por_comando"]["optics"] == 2


def test_exportar_lleva_encabezado_explicativo(tmp_path):
    from qekit.modules import feedback as fb
    fb.registrar("algo", dir_=tmp_path)
    out = fb.exportar(tmp_path / "inc.json", dir_=tmp_path)
    doc = json.loads(Path(out).read_text())
    assert "reproducir el fallo" in doc["que_es"]
    assert doc["incidencias"] and doc["estadisticas"]["total"] == 1


# ----------------------------------------------------------------------
# Parser de mallas (bug encontrado por el propio registro)
# ----------------------------------------------------------------------
def test_malla_exige_tres_numeros():
    """'1x2' llegaba al fondo y reventaba con IndexError."""
    from qekit.cli import _malla
    assert _malla("8x8x8") == (8, 8, 8)
    assert _malla("8,8,8") == (8, 8, 8)
    assert _malla(None) is None
    with pytest.raises(ValueError, match="TRES numeros"):
        _malla("1x2")
    with pytest.raises(ValueError, match="enteros"):
        _malla("8x8xocho")
    with pytest.raises(ValueError, match="positiva"):
        _malla("0x8x8")


# ----------------------------------------------------------------------
# Recomendador
# ----------------------------------------------------------------------
def _fila(formula="Si2", ecutwfc=60.0, ecutrho=480.0, kdensity=200.0,
          n_scf=12, natoms=2, convergido=1):
    return {"formula": formula, "ecutwfc": ecutwfc, "ecutrho": ecutrho,
            "kdensity": kdensity, "n_scf": n_scf, "natoms": natoms,
            "convergido": convergido}


def test_recomendador_toma_el_cutoff_maximo_no_la_media():
    """Un cutoff bajo que funcionó en un sistema no garantiza otro."""
    from qekit.modules import recommend as rc
    filas = [_fila(ecutwfc=40.0), _fila(ecutwfc=60.0), _fila(ecutwfc=80.0)]
    sug = {s.campo: s for s in rc.sugerir(filas, ["Si"])}
    assert sug["ecutwfc"].valor == 80.0
    assert sug["ecutwfc"].n_casos == 3


def test_recomendador_ignora_los_no_convergidos():
    from qekit.modules import recommend as rc
    filas = [_fila(ecutwfc=60.0), _fila(ecutwfc=200.0, convergido=0)]
    sug = {s.campo: s for s in rc.sugerir(filas, ["Si"])}
    assert sug["ecutwfc"].valor == 60.0


def test_recomendador_no_inventa_sin_historial():
    from qekit.modules import recommend as rc
    sug = rc.sugerir([_fila(formula="Si2")], ["B", "N"])
    assert sug[0].campo == "(sin historial)"
    assert "pseudopotencial" in sug[0].razon
    assert "no una predicción" in sug[0].razon


def test_recomendador_marca_la_confianza():
    from qekit.modules import recommend as rc
    uno = rc.sugerir([_fila()], ["Si"])
    assert {s.confianza for s in uno if s.n_casos == 1} == {"baja"}
    muchos = rc.sugerir([_fila() for _ in range(10)], ["Si"])
    assert [s for s in muchos if s.campo == "ecutwfc"][0].confianza == "alta"
    texto = rc.report(uno, ["Si"], 1)
    assert "UN SOLO CASO" in texto


def test_recomendador_avisa_en_losas():
    from qekit.modules import recommend as rc
    sug = rc.sugerir([_fila() for _ in range(4)], ["Si"], es_losa=True)
    mixing = [s for s in sug if s.campo == "mixing_beta"]
    assert mixing and "oscilación de carga" in mixing[0].razon
    assert mixing[0].n_casos == 0        # es regla general, no historial
