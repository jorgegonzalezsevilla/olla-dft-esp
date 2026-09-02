"""Ejecución de barridos: paralelo, reanudación y presupuesto de tiempo."""

import time

import pytest

from qekit.cli import _duracion
from qekit.core import runner
from qekit.core.errors import ErrorDeUso


# ----------------------------------------------------------------------
# reparto de recursos
# ----------------------------------------------------------------------
def test_reparto_divide_los_hilos_entre_los_calculos(monkeypatch):
    monkeypatch.setattr(runner, "nucleos", lambda: 8)
    par, n, aviso = runner.reparto(4)
    assert (par, n) == (4, 2), "8 hilos entre 4 calculos = 2 cada uno"
    assert not aviso


def test_reparto_respeta_un_nproc_explicito(monkeypatch):
    monkeypatch.setattr(runner, "nucleos", lambda: 8)
    par, n, _ = runner.reparto(2, nproc=3)
    assert (par, n) == (2, 3)


def test_reparto_avisa_de_la_sobresuscripcion(monkeypatch):
    monkeypatch.setattr(runner, "nucleos", lambda: 4)
    _, _, aviso = runner.reparto(4, nproc=4)
    assert aviso and "16" in aviso and "4 hilo" in aviso


def test_reparto_nunca_baja_de_un_proceso(monkeypatch):
    monkeypatch.setattr(runner, "nucleos", lambda: 2)
    par, n, _ = runner.reparto(8)
    assert par == 8 and n == 1


def test_nucleos_es_positivo():
    assert runner.nucleos() >= 1


# ----------------------------------------------------------------------
# duraciones
# ----------------------------------------------------------------------
@pytest.mark.parametrize("texto,segundos", [
    ("90m", 5400), ("2h", 7200), ("1h30m", 5400), ("3600", 3600),
    ("45s", 45), ("1.5h", 5400), ("2 h", 7200),
])
def test_duracion(texto, segundos):
    assert _duracion(texto) == pytest.approx(segundos)


def test_duracion_vacia():
    assert _duracion(None) is None


@pytest.mark.parametrize("malo", ["mañana", "-5m", "0", "dos horas"])
def test_duracion_rechaza_lo_que_no_es(malo):
    with pytest.raises(ErrorDeUso):
        _duracion(malo)


# ----------------------------------------------------------------------
# reanudación: qué cuenta como "ya hecho"
# ----------------------------------------------------------------------
def _job(tmp_path, texto_salida):
    d = tmp_path / "c1"
    d.mkdir()
    (d / "pw.out").write_text(texto_salida)
    return runner.Job(name="prueba", directory=d)


def test_sin_salida_no_esta_hecho(tmp_path):
    d = tmp_path / "vacia"; d.mkdir()
    assert not runner.Job(name="x", directory=d).is_done()


def test_salida_a_medias_no_esta_hecha(tmp_path):
    job = _job(tmp_path, "     iteration #  3     ecut= 30.00 Ry\n")
    assert not job.is_done()
    assert not job.is_done(estricto=False)


def test_job_done_sin_xml_no_cuenta_como_hecho(tmp_path):
    """Sin XML no se puede comprobar la convergencia: no vale por terminado."""
    job = _job(tmp_path, "     JOB DONE.\n")
    assert job.is_done(estricto=False), "en modo laxo basta la marca"
    assert not job.is_done(), "en modo estricto hace falta el XML"


def test_calculo_no_convergido_no_cuenta_como_hecho(tmp_path, monkeypatch):
    """Un scf que agota electron_maxstep escribe JOB DONE igualmente."""
    class Falso:
        converged = False
    monkeypatch.setattr(runner.qeout, "read_xml", lambda *a, **k: Falso())
    job = _job(tmp_path, "     convergence NOT achieved\n     JOB DONE.\n")
    assert not job.is_done(), "terminar no es lo mismo que converger"


def test_calculo_convergido_si_cuenta(tmp_path, monkeypatch):
    class Falso:
        converged = True
    monkeypatch.setattr(runner.qeout, "read_xml", lambda *a, **k: Falso())
    job = _job(tmp_path, "     JOB DONE.\n")
    assert job.is_done()


def test_convergencia_desconocida_no_cuenta_como_hecho(tmp_path, monkeypatch):
    class Falso:
        converged = None
    monkeypatch.setattr(runner.qeout, "read_xml", lambda *a, **k: Falso())
    job = _job(tmp_path, "     JOB DONE.\n")
    assert not job.is_done()


def _run_one_simulado(tmp_path, monkeypatch, returncode, converged):
    d = tmp_path / "ejecucion"
    d.mkdir()
    (d / "pw.in").write_text("&CONTROL\n/\n")
    job = runner.Job(name="simulado", directory=d)

    def falso_run(*args, stdout, **kwargs):
        stdout.write("     JOB DONE.\n")
        stdout.flush()

        class Proc:
            pass

        proc = Proc()
        proc.returncode = returncode
        return proc

    class Resultado:
        total_energy = -10.0

    Resultado.converged = converged
    monkeypatch.setattr(runner.subprocess, "run", falso_run)
    monkeypatch.setattr(runner.qeout, "read_xml", lambda *a, **k: Resultado())
    return runner.run_one(job, ["pw.x"])


def test_codigo_no_cero_siempre_es_fallo_aunque_haya_job_done(
        tmp_path, monkeypatch):
    res = _run_one_simulado(tmp_path, monkeypatch, returncode=7,
                            converged=True)
    assert not res.ok
    assert "código 7" in res.error


def test_run_one_rechaza_convergencia_desconocida(tmp_path, monkeypatch):
    res = _run_one_simulado(tmp_path, monkeypatch, returncode=0,
                            converged=None)
    assert not res.ok
    assert "no confirma la convergencia" in res.error


# ----------------------------------------------------------------------
# run_all: orden, paralelismo y presupuesto
# ----------------------------------------------------------------------
def _jobs(tmp_path, n):
    fuera = []
    for i in range(n):
        d = tmp_path / f"j{i}"; d.mkdir()
        (d / "pw.in").write_text("&CONTROL\n/\n")
        fuera.append(runner.Job(name=f"punto {i}", directory=d,
                                meta={"i": i}))
    return fuera


def _falso_run_one(retardos):
    """run_one de mentira que duerme lo que le digan y siempre va bien."""
    def _f(job, cmd, timeout=None, rehacer=False):
        t = retardos[job.meta["i"]]
        time.sleep(t)
        r = runner.JobResult(job=job, ok=True, seconds=t)
        return r
    return _f


def test_los_resultados_vuelven_en_orden_aunque_terminen_al_reves(
        tmp_path, monkeypatch):
    jobs = _jobs(tmp_path, 4)
    # el primero tarda mucho y el ultimo nada: terminan al reves
    monkeypatch.setattr(runner, "run_one",
                        _falso_run_one([0.30, 0.20, 0.10, 0.01]))
    monkeypatch.setattr(runner, "check_available", lambda *a, **k: "pw.x")
    res = runner.run_all(jobs, paralelo=4, verbose=False)
    assert [r.job.meta["i"] for r in res] == [0, 1, 2, 3]


def test_el_paralelo_tarda_menos_que_la_serie(tmp_path, monkeypatch):
    jobs = _jobs(tmp_path, 4)
    monkeypatch.setattr(runner, "run_one", _falso_run_one([0.25] * 4))
    monkeypatch.setattr(runner, "check_available", lambda *a, **k: "pw.x")
    t0 = time.time(); runner.run_all(jobs, paralelo=1, verbose=False)
    serie = time.time() - t0
    t0 = time.time(); runner.run_all(jobs, paralelo=4, verbose=False)
    par = time.time() - t0
    assert par < serie * 0.6, f"serie {serie:.2f}s, paralelo {par:.2f}s"


def test_el_presupuesto_deja_de_lanzar(tmp_path, monkeypatch):
    jobs = _jobs(tmp_path, 6)
    monkeypatch.setattr(runner, "run_one", _falso_run_one([0.12] * 6))
    monkeypatch.setattr(runner, "check_available", lambda *a, **k: "pw.x")
    res = runner.run_all(jobs, paralelo=1, verbose=False, presupuesto=0.25)
    lanzados = [r for r in res if r.ok]
    sin_lanzar = [r for r in res if "presupuesto" in (r.error or "")]
    assert lanzados, "alguno tiene que haberse lanzado"
    assert sin_lanzar, "y el resto tiene que quedarse sin lanzar"
    assert len(lanzados) + len(sin_lanzar) == 6


def test_el_presupuesto_no_mata_lo_que_ya_corre(tmp_path, monkeypatch):
    """Matar un pw.x a media iteracion deja la carpeta irreanudable."""
    jobs = _jobs(tmp_path, 3)
    monkeypatch.setattr(runner, "run_one", _falso_run_one([0.4, 0.4, 0.4]))
    monkeypatch.setattr(runner, "check_available", lambda *a, **k: "pw.x")
    res = runner.run_all(jobs, paralelo=3, verbose=False, presupuesto=0.05)
    # los tres arrancaron a la vez: los tres tienen que haber terminado
    assert all(r.ok for r in res)


def test_run_all_usa_la_ruta_resuelta_del_ejecutable(tmp_path, monkeypatch):
    jobs = _jobs(tmp_path, 1)
    visto = []

    def falso_run_one(job, cmd, timeout=None, rehacer=False):
        visto.append(cmd)
        return runner.JobResult(job=job, ok=True)

    monkeypatch.setattr(runner, "build_command", lambda *a, **k: ["pw.x"])
    monkeypatch.setattr(runner, "check_available",
                        lambda *a, **k: "/opt/qe/bin/pw.x")
    monkeypatch.setattr(runner, "run_one", falso_run_one)
    runner.run_all(jobs, verbose=False)
    assert visto == [["/opt/qe/bin/pw.x"]]
