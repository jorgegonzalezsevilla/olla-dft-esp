"""Pruebas puras del proyecto reproducible: workflow, campañas, resultados y dashboard."""

import json
from pathlib import Path

import pytest

from qekit.cli import build_parser
from qekit.modules import (campaign, compare, dashboard, environment,
                            onboarding, project, quality, results, tuning,
                            uncertainty, validation)


def _project(tmp_path):
    root, data = project.init(tmp_path / "demo", "Demo")
    source = root / "estructura.cif"
    source.write_text("data_Si\n_cell_length_a 5.4\n", encoding="utf-8")
    project.add_source(root, data, source)
    project.save(root, data)
    return root, data, source


def test_parser_no_confunde_project_command_con_subcomando():
    args = build_parser().parse_args(
        ["project", "plan", "dos", "--command", "olla-dft info x.cif"])
    assert args.command == "project"
    assert args.task_commands == ["olla-dft info x.cif"]


def test_project_workflow_calidad_dashboard_y_snapshot(tmp_path):
    root, data, source = _project(tmp_path)
    tasks = project.plan(root, data, "dos")
    project.save(root, data)
    assert [task["id"] for task in tasks] == ["info", "gen-dos", "dos"]
    assert project.run(root, data, execute=False)[0][1] is None

    gate = quality.evaluate(root, data)
    assert gate["fails"] == 0
    assert gate["verdict"] == "revisar"

    html = dashboard.generate(root, data)
    snapshot = project.export_snapshot(root, data)
    html_text = html.read_text(encoding="utf-8")
    assert html.is_file() and "Puerta de calidad" in html_text
    assert "Qué hacer ahora" in html_text and "Campañas" in html_text
    assert "data-command" in html_text and "task-filter" in html_text
    assert "theme-select" in html_text and "table-wrap" in html_text
    assert "aria-live=\"polite\"" in html_text and "fallbackCopy" in html_text
    exported = json.loads(snapshot.read_text(encoding="utf-8"))
    assert exported["source_state"]["sources"] == 1
    assert exported["source_state"]["changed_sources"] == []

    source.write_text(source.read_text(encoding="utf-8") + "# edit\n", encoding="utf-8")
    assert quality.evaluate(root, data)["fails"] == 1


def test_cache_detecta_fuente_cambiada_y_persiste_invalidacion(tmp_path):
    root, data, source = _project(tmp_path)
    tasks = project.plan(root, data, "custom", ["olla-dft info estructura.cif"])
    task = tasks[0]
    task["status"] = "succeeded"
    task["input_fingerprint"] = project._task_fingerprint(root, data, task)
    project.save(root, data)

    cached = project.run(root, data, execute=False)
    assert cached[0][2] == "omitida: caché válida"

    source.write_text(source.read_text(encoding="utf-8") + "# cambio\n",
                      encoding="utf-8")
    pending = project.run(root, data, execute=False)
    assert pending[0][1] is None
    assert pending[0][0]["status"] == "pending"
    assert pending[0][0]["cache_invalidated"] is True
    _root, loaded = project.load(root)
    assert loaded["tasks"][0]["status"] == "pending"


def test_bloqueo_de_entorno_se_captura_y_detecta_deriva(tmp_path):
    root, data, _ = _project(tmp_path)
    target = environment.write(root)
    current = environment.verify(root)
    assert target.is_file() and current["ok"]

    payload = json.loads(target.read_text(encoding="utf-8"))
    payload["qekit_version"] = "0.0.0-test"
    target.write_text(json.dumps(payload), encoding="utf-8")
    changed = environment.verify(root)
    assert not changed["ok"] and "qekit_version" in changed["changed"]


def test_diff_compara_fuentes_tareas_y_campanas_sin_modificar(tmp_path):
    left = {"name": "izq", "sources": [{"path": "a.cif", "sha256": "1"}],
            "tasks": [{"id": "a", "status": "pending", "command": "olla-dft info a.cif"}],
            "campaigns": []}
    right = {"name": "der", "sources": [{"path": "a.cif", "sha256": "2"}],
             "tasks": [{"id": "a", "status": "succeeded", "command": "olla-dft info a.cif"},
                       {"id": "b", "status": "pending", "command": "olla-dft info b.cif"}],
             "campaigns": [{"id": "grid"}]}
    result = project.diff(left, right)
    assert result["source_changes"] == [{"path": "a.cif", "change": "hash_changed"}]
    assert {item["change"] for item in result["task_changes"]} == {"modified", "added"}
    assert result["campaigns_left"] == 0 and result["campaigns_right"] == 1


def test_inicio_guiado_prepara_workflow_y_entorno(tmp_path):
    structure = Path(__file__).parent / "datos" / "Si_relajado.cif"
    result = onboarding.guide(tmp_path / "guided", structure_path=structure,
                              goal="gap", interactive=False)
    assert result["created"] and result["goal"] == "gap"
    assert [task["id"] for task in result["tasks"]] == ["info", "gen-bands", "bands"]
    assert (result["root"] / ".qekit" / environment.LOCK_NAME).is_file()
    assert any(item["code"] == "structure.geometry" for item in result["validation"])


def test_inicio_guiado_relax_no_se_confunde_con_scf(tmp_path):
    structure = Path(__file__).parent / "datos" / "Si_relajado.cif"
    result = onboarding.guide(tmp_path / "relax", structure_path=structure,
                              goal="relax", interactive=False)
    commands = [task["command"] for task in result["tasks"]]
    assert any("--preset relax" in command for command in commands)
    assert not any("--preset scf" in command for command in commands)


def test_inicio_guiado_en_usa_archivo_de_idioma_equivalente(tmp_path):
    structure = Path(__file__).parent / "datos" / "Si_relajado.cif"
    result = onboarding.guide(tmp_path / "guided-en", structure_path=structure,
                              goal="gap", interactive=False, language="en")
    text = onboarding.report(result)
    assert result["language"] == "en"
    assert "Recommended next step" in text and "Tasks prepared" in text


def test_plan_personalizado_es_idempotente_y_con_prefijo(tmp_path):
    root, data, _ = _project(tmp_path)
    assert project.plan(root, data, "scf", None)
    data["tasks"] = []
    first = project.plan(root, data, "custom", ["olla-dft info estructura.cif"],
                         task_prefix="review-abc")
    second = project.plan(root, data, "custom", ["olla-dft info estructura.cif"],
                          task_prefix="review-abc")
    assert first[0]["id"] == "review-abc-1"
    assert second[0]["id"] == "review-abc-1"
    assert len(data["tasks"]) == 1


def test_comparador_no_resta_corridas_incompatibles(monkeypatch):
    from tests.test_datos_qe_extra import _run

    runs = [_run("a", total_energy=-10.0), _run("b", total_energy=-9.0)]
    monkeypatch.setattr(compare.audit, "collect", lambda paths: runs)
    result = compare.compare(["a", "b"])
    assert result["comparable_energy"]
    assert result["runs"][1]["delta_energia_eV"] == pytest.approx(1.0)
    runs[1] = _run("b", total_energy=-9.0, ecutwfc=80.0)
    result = compare.compare(["a", "b"])
    assert not result["comparable_energy"]
    assert result["runs"][1]["delta_energia_eV"] is None
    assert "NO se restan" in compare.report(result)


def test_tuning_recomienda_confirmar_o_extender(tmp_path):
    data = tmp_path / "CONVERGENCIA.dat"
    data.write_text("# valor E dE\n30 -10 4\n40 -10.1 1.5\n50 -10.2 0.4\n",
                    encoding="utf-8")
    result = tuning.analyze(data, threshold=1.0)
    assert result["status"] == "confirm"
    assert result["recommended_value"] > 50
    assert "CONFIRMAR" in tuning.report(result)


def test_validacion_avanzada_lee_estructura_y_detecta_colision(tmp_path):
    root, data, source = _project(tmp_path)
    cif = (Path(__file__).parent / "datos" / "Si_relajado.cif").read_text(encoding="utf-8")
    source.write_text(cif, encoding="utf-8")
    project.add_source(root, data, source)
    data["tasks"] = [
        {"id": "a", "command": "olla-dft info estructura.cif", "outputs": ["out"]},
        {"id": "b", "command": "olla-dft info estructura.cif", "outputs": ["out"]},
    ]
    checks = validation.check(root, data)
    assert any(item["code"] == "structure.geometry" and item["level"] == "ok"
               for item in checks)
    assert any(item["code"] == "task.output_collision" for item in checks)


def test_busqueda_db_parametrizada_no_necesita_sql(tmp_path):
    from qekit.modules import audit
    from tests.test_datos_qe_extra import _run

    db = tmp_path / "history.db"
    audit.index([_run("a"), _run("b", calculation="relax")], db)
    rows = audit.search(db, formula="Si", calculation="relax", limit=5)
    assert len(rows) == 1
    assert rows[0]["calculation"] == "relax"


def test_resultados_normalizados_son_idempotentes_y_guardan_unidades(tmp_path,
                                                                      monkeypatch):
    from tests.test_datos_qe_extra import _run

    output = tmp_path / "scf"
    output.mkdir()
    (output / "pw.xml").write_text("version-1", encoding="utf-8")
    run = _run(str(output), total_energy=-20.0)
    monkeypatch.setattr(results.audit, "collect", lambda paths: [run])
    db = tmp_path / "results.sqlite3"
    first = results.ingest([output], db, tag="scf-1")
    second = results.ingest([output], db, tag="scf-1")
    rows = results.list_results(db)
    assert first["inserted"] == 1 and second["existing"] == 1
    assert len(rows) == 1
    assert rows[0]["status"] == "converged"
    assert rows[0]["metrics"]["energy_total"]["unit"] == "eV"
    assert rows[0]["metrics"]["energy_per_atom"]["value"] == pytest.approx(-10.0)


def test_resultado_modificado_crea_nueva_evidencia(tmp_path, monkeypatch):
    from tests.test_datos_qe_extra import _run

    output = tmp_path / "scf"
    output.mkdir()
    xml = output / "pw.xml"
    xml.write_text("version-1", encoding="utf-8")
    monkeypatch.setattr(results.audit, "collect",
                        lambda paths: [_run(str(output), total_energy=-20.0)])
    db = tmp_path / "results.sqlite3"
    results.ingest([output], db)
    xml.write_text("version-2", encoding="utf-8")
    results.ingest([output], db)
    assert len(results.list_results(db)) == 2


def test_campana_crea_puntos_independientes_y_validos(tmp_path):
    root, data, _ = _project(tmp_path)
    record = campaign.create(
        root, data, "malla",
        "olla-dft gen {structure} --preset scf --outdir artifacts/{id}",
        ["ecut=30,40", "k=2x2x2,4x4x4"], goal="convergencia")
    project.save(root, data)
    assert record["points"] == 4
    tasks = [task for task in data["tasks"] if task.get("campaign_id") == "malla"]
    assert len(tasks) == 4
    assert all(not task["depends_on"] for task in tasks)
    assert all(task["command"].startswith("olla-dft gen") for task in tasks)
    assert all(task["outputs"] == [f"artifacts/point-{i:03d}"]
               for i, task in enumerate(tasks, 1))
    assert campaign.status(data, "malla")["counts"]["pending"] == 4


def test_campana_adaptativa_añade_recomendacion(tmp_path):
    root, data, _ = _project(tmp_path)
    convergence = tmp_path / "CONVERGENCIA.dat"
    convergence.write_text("# v e d\n30 -1 4\n40 -2 0.2\n", encoding="utf-8")
    record = campaign.create(
        root, data, "cutoff",
        "olla-dft gen {structure} --ecutwfc {ecutwfc} --outdir artifacts/{id}",
        ["ecutwfc=30,40"], convergence_file=convergence, adaptive=True)
    assert record["axes"]["ecutwfc"][-1] > 40


def test_campana_extend_aplica_un_punto_nuevo_sin_rehacer_los_previos(tmp_path):
    root, data, _ = _project(tmp_path)
    record = campaign.create(
        root, data, "cutoff", "olla-dft gen {structure} --ecutwfc {ecutwfc} "
        "--outdir artifacts/{id}", ["ecutwfc=30,40"])
    convergence = tmp_path / "CONVERGENCIA.dat"
    convergence.write_text("# v e d\n30 -1 4\n40 -2 0.2\n", encoding="utf-8")
    result = campaign.extend(root, data, record["id"], convergence)
    assert result["extended"] and result["points_added"] == 1
    assert data["campaigns"][0]["points"] == 3
    assert len({task["parameters"]["ecutwfc"] for task in data["tasks"]}) == 3


def test_project_ingest_actualiza_manifiesto_y_calidad(tmp_path, monkeypatch):
    root, data, _ = _project(tmp_path)
    output = root / ".qekit" / "artifacts" / "scf"
    output.mkdir(parents=True)
    (output / "pw.xml").write_text("xml", encoding="utf-8")
    from tests.test_datos_qe_extra import _run
    monkeypatch.setattr(results.audit, "collect",
                        lambda paths: [_run(str(output))])
    record = results.ingest_project(root, data)
    assert record["inserted"] == 1
    _root, loaded = project.load(root)
    gate = quality.evaluate(root, loaded)
    assert any(check.code == "results.indexed" for check in gate["checks"])


def test_migracion_de_manifiesto_y_lock_no_pierde_datos(tmp_path):
    root = tmp_path / "legacy"
    directory = root / ".qekit"
    directory.mkdir(parents=True)
    target = directory / "project.json"
    target.write_text(json.dumps({
        "schema_version": 1, "name": "legacy", "sources": [],
        "tasks": [{"id": "old", "status": "pending"}],
    }), encoding="utf-8")
    _root, data = project.load(root)
    assert data["schema_version"] == project.SCHEMA_VERSION
    assert data["tasks"][0]["id"] == "old"
    assert data["metadata"]["migrations"][0]["from"] == 1
    assert json.loads(target.read_text(encoding="utf-8"))["schema_version"] == 2
    assert (directory / ".project.lock").is_file()


def test_run_paralelo_ejecuta_puntos_independientes(monkeypatch, tmp_path):
    root, data, _ = _project(tmp_path)
    data["tasks"] = [
        {"id": "a", "command": "olla-dft info estructura.cif", "depends_on": [],
         "outputs": [], "status": "pending"},
        {"id": "b", "command": "olla-dft info estructura.cif", "depends_on": [],
         "outputs": [], "status": "pending"},
    ]

    class Completed:
        returncode = 0
        stdout = "ok"
        stderr = ""

    monkeypatch.setattr(project.subprocess, "run",
                        lambda *args, **kwargs: Completed())
    result = project.run(root, data, execute=True, parallel=2)
    assert len(result) == 2 and all(item[0]["status"] == "succeeded" for item in result)


def test_run_registra_reintento_y_recupera_fallo_transitorio(monkeypatch, tmp_path):
    root, data, _ = _project(tmp_path)
    data["tasks"] = [{"id": "a", "command": "olla-dft info estructura.cif",
                       "depends_on": [], "outputs": [], "status": "pending"}]
    calls = []

    class Completed:
        stdout = ""
        stderr = ""

        def __init__(self, code):
            self.returncode = code

    def fake_run(*args, **kwargs):
        calls.append(kwargs)
        return Completed(1 if len(calls) == 1 else 0)

    monkeypatch.setattr(project.subprocess, "run", fake_run)
    result = project.run(root, data, execute=True, retries=1, timeout=12)
    task = result[0][0]
    assert task["status"] == "succeeded" and task["retry_count"] == 1
    assert [attempt["returncode"] for attempt in task["attempts"]] == [1, 0]
    assert all(item["timeout"] == 12 for item in calls)


def test_cancelacion_cooperativa_y_reanudacion(tmp_path):
    root, data, _ = _project(tmp_path)
    data["tasks"] = [
        {"id": "a", "command": "olla-dft info estructura.cif",
         "depends_on": [], "outputs": [], "status": "pending"},
        {"id": "b", "command": "olla-dft gen estructura.cif",
         "depends_on": ["a"], "outputs": [], "status": "pending"},
    ]
    marker = project.cancel(root, "pausa del usuario")
    result = project.run(root, data, execute=True)
    assert marker.is_file() and all(item[0]["status"] == "cancelled" for item in result)
    assert project.status(root, data)["counts"]["cancelled"] == 2
    assert project.resume(root, data) == 2
    assert not marker.exists() and all(task["status"] == "pending" for task in data["tasks"])


def test_cancelacion_personalizada_resuelve_ruta_relativa_al_proyecto(tmp_path):
    root, data, _ = _project(tmp_path)
    custom = project.cancel(root, "pausa", "markers/STOP.json")
    assert custom == root / "markers" / "STOP.json"
    assert project.run(root, data, execute=True, cancel_file="markers/STOP.json") == []
    project.resume(root, data, "markers/STOP.json")
    assert not custom.exists()


def test_resultados_revision_humana_y_migracion_de_esquema(tmp_path, monkeypatch):
    from tests.test_datos_qe_extra import _run
    output = tmp_path / "scf"
    output.mkdir()
    (output / "pw.xml").write_text("xml", encoding="utf-8")
    monkeypatch.setattr(results.audit, "collect", lambda paths: [_run(str(output))])
    db = tmp_path / "results.sqlite3"
    results.ingest([output], db)
    row = results.list_results(db)[0]
    reviewed = results.review(db, row["id"], "accepted", "revisado por QA")
    assert reviewed["review"]["status"] == "accepted"
    assert reviewed["review"]["note"] == "revisado por QA"


def test_incertidumbre_es_determinista(tmp_path):
    result = uncertainty.propagate(lambda values: values[0] ** 2 + values[1],
                                    [2.0, 3.0], [0.1, 0.2])
    assert result["value"] == pytest.approx(7.0)
    assert result["uncertainty"] == pytest.approx((0.4 ** 2 + 0.2 ** 2) ** 0.5,
                                                   rel=1e-4)


def test_dashboard_admite_tema_y_pdf(tmp_path):
    root, data, _ = _project(tmp_path)
    html = dashboard.generate(root, data, theme="dark")
    html_text = html.read_text(encoding="utf-8")
    assert 'data-theme="dark"' in html_text
    assert '<option value="dark" selected>Oscuro</option>' in html_text
    english = dashboard.generate(root, data, tmp_path / "dashboard-en.html",
                                language="en")
    english_text = english.read_text(encoding="utf-8")
    assert '<html lang="en">' in english_text and "Getting started" in english_text
    pair = dashboard.generate_pair(root, data, tmp_path / "dashboard.html")
    assert pair[0].name == "dashboard.html" and pair[1].name == "dashboard.en.html"
    assert '<html lang="es">' in pair[0].read_text(encoding="utf-8")
    assert '<html lang="en">' in pair[1].read_text(encoding="utf-8")
    from qekit.modules import report as project_report
    pdf = project_report.generate_pdf(root, data)
    assert pdf.is_file() and pdf.stat().st_size > 0


def test_traducciones_es_en_mantienen_paridad_de_interfaz():
    assert set(dashboard._labels("es")) == set(dashboard._labels("en"))
    assert dashboard._labels("es")["other_language"] == "English"
    assert dashboard._labels("en")["other_language"] == "Español"


def test_menu_inicial_tiene_variante_ingles_separada():
    from qekit.cli import _menu_labels

    spanish = _menu_labels("es")
    english = _menu_labels("en")
    assert spanish.keys() == english.keys()
    assert len(spanish["items"]) == len(english["items"])
    assert english["choice"] == "Choice"
    assert spanish["submenus"].keys() == english["submenus"].keys()
    for name in spanish["submenus"]:
        assert spanish["submenus"][name].keys() == english["submenus"][name].keys()
