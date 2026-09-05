"""Portable result exports keep scientific values, scope, and safe HTML."""
import json
import re

import pytest

from qekit.core.errors import ErrorDeUso
from qekit.modules import results, studio, project, dashboard


def row(identity='one', **extra):
    return {'id': identity, 'formula': 'Si2', 'status': 'converged',
            'calculation': 'scf', 'metrics': {'energy_total': {'value': -20.125, 'unit': 'eV'}},
            'provenance': {'fingerprint': ['PBE'], 'files': {'pw.xml': 'abc'},
                           'parameters': {'kgrid': [2, 2, 2]}}, **extra}


def payload(path):
    return json.loads(re.search(r'<script id="studio-data" type="application/json">(.*?)</script>', path.read_text(), re.S).group(1))


def test_portable_values_units_uncertainty_and_paths():
    original = row(path='/private/calc', review={'status': 'accepted', 'note': '/private/note'})
    original['metrics'].update(zero={'value': 0, 'unit': 'eV', 'uncertainty': 0},
                               absent={'value': float('nan'), 'unit': 'GPa'},
                               infinity={'value': float('inf'), 'unit': 's'})
    out = studio.portable_rows([original])[0]
    assert out['metrics']['energy_total'] == {'value': -20.125, 'unit': 'eV'}
    assert out['metrics']['zero'] == {'value': 0, 'unit': 'eV', 'uncertainty': 0}
    assert out['metrics']['absent'] == {'value': None, 'unit': 'GPa', 'reason': 'missing_or_nonfinite'}
    assert out['metrics']['infinity']['value'] is None
    assert out['review'] == 'accepted'
    assert '/private/' not in json.dumps(out)
    json.dumps(out, allow_nan=False)


def test_method_identity_includes_kgrid_and_preserves_distinct_runs():
    first = row('first'); second = row('second')
    second['provenance']['parameters']['kgrid'] = [4, 4, 4]
    a, b = studio.portable_rows([first, second])
    assert a['method_sha256'] != b['method_sha256']
    assert a['source_sha256'] == b['source_sha256']
    with pytest.raises(ErrorDeUso, match='Duplicate'):
        studio.portable_rows([first, first])


def test_embedding_escapes_script_and_never_expands_user_tokens(tmp_path):
    title = '@@JS@@ </script><img src=x onerror=alert(1)>'
    path = studio.generate([row(formula=title)], tmp_path/'view.html', title=title)
    decoded = payload(path)
    assert decoded['title'] == title
    assert decoded['rows'][0]['formula'] == title
    assert '<img src=x onerror=alert(1)>' not in path.read_text()
    assert '@@JS@@ &lt;/script&gt;' in path.read_text()
    assert len(re.findall(r'<script\b', path.read_text())) == 2


def test_labels_parity_and_snapshot_scope(tmp_path):
    path = studio.generate([row()], tmp_path/'limited.html', language='en', total_count=10)
    data = payload(path)
    assert set(data['labels']['en']) == set(data['labels']['es'])
    assert data['language'] == 'en' and data['total_count'] == 10
    assert len(data['rows']) == 1
    for key in re.findall(r'data-(?:label|aria)="([^"]+)"', path.read_text()):
        assert key in data['labels']['en']


def make_db(tmp_path, monkeypatch):
    from tests.test_datos_qe_extra import _run
    runs=[]
    for i, calc in enumerate(['scf', 'relax', 'scf']):
        folder = tmp_path/f'run-{i}'; folder.mkdir(); (folder/'pw.xml').write_text(str(i))
        runs.append(_run(str(folder), calculation=calc, total_energy=-20-i, kgrid=(2,2,2)))
    monkeypatch.setattr(results.audit, 'collect', lambda paths: runs)
    db = tmp_path/'results.sqlite3'
    results.ingest([tmp_path], db)
    return db


def test_filtered_count_and_recorded_parameters(tmp_path, monkeypatch):
    db = make_db(tmp_path, monkeypatch)
    assert results.count_results(db) == 3
    assert results.count_results(db, formula='Si', calculation='SCF', status='converged') == 2
    assert results.count_results(db, formula="' OR 1=1 --") == 0
    assert len(results.list_results(db, calculation='scf', limit=1)) == 1
    assert results.list_results(db)[0]['provenance']['parameters']['kgrid'] == [2,2,2]


def test_explore_cli_standalone_database(tmp_path, monkeypatch):
    from qekit.cli import main
    db = make_db(tmp_path, monkeypatch)
    out = tmp_path/'explore.html'
    assert main(['--language', 'en', 'results', 'explore', '--db', str(db), '--project', str(tmp_path), '--calculation', 'scf', '--limit', '1', '-o', str(out)]) == 0
    data = payload(out)
    assert data['total_count'] == 2 and len(data['rows']) == 1
    assert data['language'] == 'en'


def test_empty_dashboard_links_to_offline_explorer(tmp_path):
    root, data = project.init(tmp_path, name='Empty')
    target = dashboard.generate(root, data, destination=tmp_path/'custom name.html')
    assert 'custom%20name.results.html' in target.read_text()
    assert payload(tmp_path/'custom name.results.html')['rows'] == []


def test_html_incluye_licencia_y_fuente_de_la_version(tmp_path):
    import html
    from qekit import __version__
    path = studio.generate([row()], tmp_path / 'licensed.html')
    text = path.read_text()
    assert 'GNU AFFERO GENERAL PUBLIC LICENSE' in html.unescape(text)
    assert '13. Remote Network Interaction' in text
    assert f'{studio.SOURCE_URL}/tree/v{__version__}' in text
    assert '@@LICENSE@@' not in text and '@@SOURCE@@' not in text
    assert payload(path)['rows'][0]['metrics']['energy_total']['value'] == -20.125
