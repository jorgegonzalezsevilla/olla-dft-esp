"""Offline, portable result explorer. Presentation never changes source metrics."""
import hashlib
import html
import json
import math
import re
from pathlib import Path

from qekit import __version__
from qekit.core.errors import ErrorDeUso
from qekit.modules import results

ASSETS = Path(__file__).resolve().parents[1] / 'data'
SOURCE_URL = 'https://github.com/jorgegonzalezsevilla/olla-dft-esp'


def portable_rows(rows):
    """Keep scientific identity/values, omit local paths and free-form provenance."""
    output = []
    seen = set()
    for row in rows:
        identity = str(row['id'])
        if identity in seen:
            raise ErrorDeUso('Duplicate result id in explorer input.')
        seen.add(identity)
        metrics = {}
        for key, metric in row.get('metrics', {}).items():
            if not isinstance(metric, dict):
                continue
            value = metric.get('value')
            finite = isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)
            metrics[str(key)] = {'value': value if finite else None, 'unit': str(metric.get('unit') or '')}
            if not finite:
                metrics[str(key)]['reason'] = 'missing_or_nonfinite'
            uncertainty = metric.get('uncertainty')
            if isinstance(uncertainty, (int, float)) and not isinstance(uncertainty, bool) and math.isfinite(uncertainty) and uncertainty >= 0:
                metrics[str(key)]['uncertainty'] = uncertainty
        provenance = row.get('provenance', {})
        def checksum(data):
            return hashlib.sha256(json.dumps(data, sort_keys=True, ensure_ascii=False, default=str).encode()).hexdigest()
        output.append({
            'id': identity, 'index': len(output)+1,
            'formula': str(row.get('formula') or ''), 'calculation': str(row.get('calculation') or ''),
            'status': str(row.get('status') or 'invalid'), 'converged': row.get('converged'),
            'tag': str(row.get('tag') or ''),
            'review': str((row.get('review') or {}).get('status') or 'unreviewed'),
            'metrics': metrics,
            'source_sha256': checksum(provenance.get('files', {})),
            'method_sha256': checksum({'fingerprint': provenance.get('fingerprint', []), 'parameters': provenance.get('parameters', {})}),
            'method_known': bool(provenance.get('fingerprint')),
            'parameters': {key: value for key, value in provenance.get('parameters', {}).items()
                           if key in ('functional', 'ecutwfc_Ry', 'ecutrho_Ry', 'kgrid', 'kshift', 'smearing', 'degauss_Ry', 'occupations', 'nspin')},
            'ingested': str(row.get('ingested') or ''),
        })
    return output


def generate(rows, destination, title='Olla-DFT', language='es', total_count=None, order='input_order'):
    if language not in ('es', 'en'):
        raise ErrorDeUso('Explorer language must be es or en.')
    records = portable_rows(rows)
    labels = {lang: json.loads((ASSETS/'i18n'/f'studio_{lang}.json').read_text()) for lang in ('es', 'en')}
    payload = dict(schema_version=1, qekit_version=__version__, generated=results._now(),
                   title=str(title), language=language, total_count=total_count if total_count is not None else len(records),
                   rows=records, labels=labels, view=None, order=order)
    encoded = json.dumps(payload, ensure_ascii=False, allow_nan=False).replace('<', '\\u003c').replace('>', '\\u003e').replace('&', '\\u0026')
    template = (ASSETS/'studio/studio.html').read_text()
    substitutions = {'LANG': language, 'TITLE': html.escape(str(title)),
                     'CSS': (ASSETS/'studio/studio.css').read_text(),
                     'JS': (ASSETS/'studio/studio.js').read_text(), 'PAYLOAD': encoded,
                     'LICENSE': html.escape((ASSETS/'AGPL-3.0.txt').read_text()),
                     'SOURCE': f'{SOURCE_URL}/tree/v{__version__}'}
    text = re.sub(r'@@(LANG|TITLE|CSS|JS|PAYLOAD|LICENSE|SOURCE)@@',
                  lambda match: substitutions[match.group(1)], template)
    target = Path(destination)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text, encoding='utf-8')
    return target
