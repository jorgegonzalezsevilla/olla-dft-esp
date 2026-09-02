"""Que esto corra en Linux, en macOS y en Windows.

La prueba que justifica el archivo es
`test_todo_lo_que_imprime_qekit_se_puede_transliterar`: recorre los informes
de VERDAD que generan los módulos y comprueba que cada carácter fuera de
ASCII o está en la tabla de transliteración o es una vocal acentuada que
cp1252 admite. Sin ella, alguien añade un ⟨símbolo bonito⟩ a un informe y
QEkit vuelve a morir en la consola de Windows sin que nadie lo note hasta
que un usuario lo sufre.
"""

import ast
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from qekit.core import consola, plataforma


# ----------------------------------------------------------------------
# Codificación: el fallo que mataba comandos enteros en Windows
# ----------------------------------------------------------------------
def test_transliterar_deja_todo_en_ascii():
    muestra = "Å α ε² Ω → ← ① ✓ ± × · — ≈ ≤ ħ ∫ ⁻¹ ₂ │ └ áéíóúñ ¿¡"
    salida = consola.transliterar(muestra)
    salida.encode("ascii")            # lanza si quedó algo fuera


def test_transliterar_no_toca_lo_que_ya_era_ascii():
    t = "kappa_L = 100.7 W/m/K  (RTA, fc3)"
    assert consola.transliterar(t) == t


def test_transliterar_conserva_el_significado():
    assert consola.transliterar("5.43 Å") == "5.43 A"
    assert consola.transliterar("κ_L") == "kappa_L"
    assert consola.transliterar("① → ②") == "(1) -> (2)"
    assert consola.transliterar("✓ bien") == "ok bien"
    assert consola.transliterar("energía") == "energia"


def _todo_el_texto_de_los_informes():
    """Los informes de verdad, no una lista escrita a mano."""
    from qekit.modules import recipes, wizard, docs
    trozos = [recipes.listar(), wizard.report_catalogo(), plataforma.informe()]
    for r in recipes.RECETAS:
        trozos.append(recipes.report(r))
        trozos.append(recipes.script(r))
    for m in wizard.METAS:
        trozos.append(wizard.report_meta(m, "x.cif"))
    trozos.append("\n".join(g[1] for g in docs.GRUPOS))
    return "\n".join(trozos)


def test_todo_lo_que_imprime_qekit_se_puede_transliterar():
    texto = _todo_el_texto_de_los_informes()
    limpio = consola.transliterar(texto)
    malos = sorted({c for c in limpio if ord(c) > 127})
    # lo que sobreviva tiene que caber al menos en cp1252 (acentos, ¿, ñ)
    problemas = []
    for c in malos:
        try:
            c.encode("cp1252")
        except UnicodeEncodeError:
            problemas.append(f"{c!r} (U+{ord(c):04X})")
    assert not problemas, (
        "estos caracteres no están en la tabla de transliteración y tampoco "
        "caben en cp1252: " + ", ".join(problemas))


def test_una_consola_ascii_no_mata_el_comando(capsys):
    """El escalón 2: si no cabe, se transliteria; nunca se revienta."""
    import io

    class Duro(io.StringIO):
        encoding = "ascii"

        def write(self, t):
            t.encode("ascii")
            return super().write(t)

    duro = Duro()
    envuelto = consola._Transliterando(duro)
    envuelto.write("El parámetro de red es 5.43 Å (α = 90°)\n")
    assert "5.43 A" in duro.getvalue()
    assert "parametro" in duro.getvalue()


@pytest.mark.parametrize("comando", [
    ["recetas", "mecanicas"],
    ["sistema"],
    ["wizard", "--list"],
    ["--version"],
])
def test_ningun_comando_revienta_en_una_consola_cp1252(comando, tmp_path):
    """Lo que se rompía de verdad, comprobado lanzando el programa entero."""
    raiz = Path(__file__).resolve().parents[1]
    env = dict(os.environ, PYTHONPATH=str(raiz), PYTHONIOENCODING="cp1252",
               QEKIT_CONFIG_DIR=str(tmp_path))
    r = subprocess.run([sys.executable, "-m", "qekit.cli"] + comando,
                       capture_output=True, env=env, cwd=str(tmp_path))
    err = r.stderr.decode("utf-8", "replace")
    assert "UnicodeEncodeError" not in err, err[-400:]
    assert r.returncode == 0, err[-400:]


def test_la_bandera_ascii_funciona_este_donde_este(tmp_path):
    raiz = Path(__file__).resolve().parents[1]
    env = dict(os.environ, PYTHONPATH=str(raiz),
               QEKIT_CONFIG_DIR=str(tmp_path))
    for argv in (["--ascii", "recetas"], ["recetas", "--ascii"]):
        r = subprocess.run([sys.executable, "-m", "qekit.cli"] + argv,
                           capture_output=True, env=env, cwd=str(tmp_path))
        assert r.returncode == 0, r.stderr.decode()[-300:]
        r.stdout.decode("utf-8").encode("ascii")   # lanza si quedó algo


# ----------------------------------------------------------------------
# Dónde van las cosas en cada sistema
# ----------------------------------------------------------------------
def test_la_variable_de_entorno_manda_sobre_todo(monkeypatch, tmp_path):
    monkeypatch.delenv("OLLA_DFT_CONFIG_DIR", raising=False)
    monkeypatch.setenv("QEKIT_CONFIG_DIR", str(tmp_path / "otro"))
    assert plataforma.dir_config() == tmp_path / "otro"
    # La variable nueva manda sobre la heredada.
    monkeypatch.setenv("OLLA_DFT_CONFIG_DIR", str(tmp_path / "nuevo"))
    assert plataforma.dir_config() == tmp_path / "nuevo"


def test_config_en_windows(monkeypatch, tmp_path):
    monkeypatch.delenv("QEKIT_CONFIG_DIR", raising=False)
    monkeypatch.delenv("OLLA_DFT_CONFIG_DIR", raising=False)
    monkeypatch.setattr(plataforma, "WINDOWS", True)
    monkeypatch.setattr(plataforma, "MACOS", False)
    monkeypatch.setenv("APPDATA", str(tmp_path / "AppData" / "Roaming"))
    d = plataforma.dir_config()
    assert d.name == "olla-dft" and "Roaming" in str(d)
    assert plataforma.dirs_config_heredados()[1].name == "QEkit"


def test_config_en_macos(monkeypatch):
    monkeypatch.delenv("QEKIT_CONFIG_DIR", raising=False)
    monkeypatch.delenv("OLLA_DFT_CONFIG_DIR", raising=False)
    monkeypatch.setattr(plataforma, "WINDOWS", False)
    monkeypatch.setattr(plataforma, "MACOS", True)
    d = plataforma.dir_config()
    assert d.parts[-3:] == ("Library", "Application Support", "olla-dft")


def test_config_en_linux_respeta_xdg(monkeypatch, tmp_path):
    monkeypatch.delenv("QEKIT_CONFIG_DIR", raising=False)
    monkeypatch.delenv("OLLA_DFT_CONFIG_DIR", raising=False)
    monkeypatch.setattr(plataforma, "WINDOWS", False)
    monkeypatch.setattr(plataforma, "MACOS", False)
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "cfg"))
    assert plataforma.dir_config() == tmp_path / "cfg" / "olla-dft"


def test_datos_grandes_respetan_xdg_y_variable(monkeypatch, tmp_path):
    monkeypatch.setattr(plataforma, "WINDOWS", False)
    monkeypatch.setattr(plataforma, "MACOS", False)
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "datos"))
    monkeypatch.delenv("QEKIT_DATA_DIR", raising=False)
    monkeypatch.delenv("OLLA_DFT_DATA_DIR", raising=False)
    assert plataforma.dir_data() == tmp_path / "datos" / "olla-dft"
    monkeypatch.setenv("QEKIT_DATA_DIR", str(tmp_path / "modelo-externo"))
    assert plataforma.dir_data() == tmp_path / "modelo-externo"


def test_la_configuracion_vieja_se_migra(tmp_path, monkeypatch):
    """Quien ya usaba QEkit no puede perder sus ajustes al actualizar."""
    vieja = tmp_path / "viejo"
    vieja.mkdir()
    (vieja / "config.ini").write_text("[qekit]\npseudo_dir = /mis/pseudos\n")
    nueva = tmp_path / "nuevo"
    from qekit import config as qcfg
    monkeypatch.setattr(qcfg, "CONFIG_DIR", nueva)
    monkeypatch.setattr(qcfg, "CONFIG_FILE", nueva / "config.ini")
    monkeypatch.setattr(plataforma, "dirs_config_heredados", lambda: [vieja])
    assert qcfg.load()["pseudo_dir"] == "/mis/pseudos"
    assert (vieja / "config.ini").exists()      # no se borra el original


# ----------------------------------------------------------------------
# Binarios y lanzadores
# ----------------------------------------------------------------------
def test_en_windows_se_busca_el_exe_primero(monkeypatch):
    monkeypatch.setattr(plataforma, "WINDOWS", True)
    assert plataforma.nombres_ejecutable("pw") == ["pw.exe", "pw.x"]


def test_fuera_de_windows_manda_el_punto_x(monkeypatch):
    monkeypatch.setattr(plataforma, "WINDOWS", False)
    assert plataforma.nombres_ejecutable("pw") == ["pw.x", "pw.exe"]
    assert plataforma.nombres_ejecutable("pw.x")[0] == "pw.x"
    assert plataforma.nombres_ejecutable("pw.exe")[0] == "pw.x"


def test_sin_mpi_se_corre_en_serie_en_vez_de_fallar(monkeypatch):
    from qekit.core import runner
    monkeypatch.setattr(plataforma, "lanzador_mpi", lambda: "")
    monkeypatch.setattr(runner.qcfg, "load",
                        lambda: {"pw_cmd": "pw.x", "nproc": "4",
                                 "mpi_cmd": ""})
    assert runner.build_command(nproc=4) == ["pw.x"]


def test_un_lanzador_ya_puesto_se_respeta(monkeypatch):
    from qekit.core import runner
    monkeypatch.setattr(runner.qcfg, "load",
                        lambda: {"pw_cmd": "srun -n 8 pw.x", "nproc": "4",
                                 "mpi_cmd": "mpirun -np {n}"})
    assert runner.build_command() == ["srun", "-n", "8", "pw.x"]


# ----------------------------------------------------------------------
# Los guiones generados
# ----------------------------------------------------------------------
def test_los_guiones_se_escriben_con_finales_posix(tmp_path):
    """Un .sh con CRLF falla en WSL con «bad interpreter: /bin/bash^M»."""
    f = plataforma.escribir_script(tmp_path / "x.sh", "#!/bin/bash\necho a\n")
    crudo = f.read_bytes()
    assert b"\r\n" not in crudo


def test_el_conteo_de_nucleos_sirve_tambien_en_macos():
    trozo = plataforma.cuenta_nucleos_shell()
    assert "nproc" in trozo and "hw.ncpu" in trozo


def test_el_par_de_guiones_deja_sh_y_py_valido(tmp_path):
    hechos = plataforma.escribir_par_de_guiones(
        tmp_path, [("pw.x", "1.in", "1.out"), ("ph.x", "2.in", "2.out")])
    assert len(hechos) == 2
    sh = (tmp_path / "correr.sh").read_text()
    py = (tmp_path / "correr.py").read_text()
    assert "mpiexec" in sh and "mpirun" in sh      # los dos lanzadores
    ast.parse(py)                                   # es Python válido
    assert "pw.x" in py and "ph.x" in py
    assert ".exe" in py                             # sabe buscar el de Windows


def test_el_run_py_del_barrido_es_python_valido(tmp_path):
    from qekit.modules import sweep
    f = sweep.write_run_py(["a", "b"], tmp_path / "run.py", 2, "pw.x")
    texto = Path(f).read_text()
    ast.parse(texto)
    assert "['a', 'b']" in texto
    assert "mpiexec" in texto                       # no solo mpirun


def test_el_barrido_escribe_los_dos_lanzadores(tmp_path):
    from qekit.core.runner import Job
    from qekit.modules import sweep
    jobs = [Job(name="a", directory=tmp_path / "a"),
            Job(name="b", directory=tmp_path / "b")]
    for j in jobs:
        Path(j.directory).mkdir()
    sweep.write_run_script(jobs, tmp_path / "run.sh")
    assert (tmp_path / "run.sh").exists()
    assert (tmp_path / "run.py").exists()
    ast.parse((tmp_path / "run.py").read_text())


def test_run_sh_del_barrido_propaga_fallos(tmp_path):
    from qekit.core.runner import Job
    from qekit.modules import sweep

    d = tmp_path / "a"
    d.mkdir()
    (d / "pw.in").write_text("&CONTROL\n/\n")
    falso = tmp_path / "pw-falso"
    falso.write_text("#!/bin/sh\nexit 7\n")
    falso.chmod(0o755)
    sweep.write_run_script([Job(name="a", directory=d)],
                           tmp_path / "run.sh", nproc=1)
    env = dict(os.environ, PW=str(falso), NPROC="1")
    bash = shutil.which("bash")
    if bash is None:
        pytest.skip("esta comprobación necesita Bash")
    r = subprocess.run([bash, str(tmp_path / "run.sh")], env=env,
                       cwd=str(tmp_path), capture_output=True, text=True)
    assert r.returncode != 0


def test_ningun_guion_generado_cablea_mpirun_a_secas():
    """Si un módulo vuelve a escribir 'mpirun' sin comprobar que existe,
    los usuarios de Windows y de máquinas sin MPI se quedan sin correrlo."""
    raiz = Path(__file__).resolve().parents[1] / "qekit"
    malos = []
    for f in raiz.rglob("*.py"):
        txt = f.read_text(encoding="utf-8", errors="replace")
        for i, linea in enumerate(txt.split("\n"), 1):
            if "mpirun" not in linea or linea.strip().startswith("#"):
                continue
            # vale si en el mismo archivo se comprueba que existe
            if "command -v mpirun" in txt or "shutil.which" in txt \
                    or "lanzador_mpi" in txt or "mpi_cmd" in txt:
                continue
            malos.append(f"{f.name}:{i}")
    assert not malos, malos
