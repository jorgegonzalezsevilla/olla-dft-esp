# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Dashboard HTML autocontenido para un proyecto Olla-DFT."""

from __future__ import annotations

import html
import json
import math
from pathlib import Path

from qekit import __command_name__, __product_name__
from qekit.modules import quality, project, results


_TRANSLATION_DIR = Path(__file__).resolve().parent.parent / "data" / "i18n"


def _labels(language: str) -> dict:
    """Carga una traducción versionada; ambas variantes comparten el contrato."""
    target = _TRANSLATION_DIR / f"dashboard_{language}.json"
    try:
        value = json.loads(target.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"no se pudo cargar la traducción {target}: {exc}") from None
    if not isinstance(value, dict):
        raise RuntimeError(f"la traducción {target} no es un objeto JSON")
    return value


def generate(root: Path, data: dict, destination=None, theme="auto", language="es") -> Path:
    if theme not in ("auto", "light", "dark"):
        raise ValueError("theme debe ser auto, light o dark")
    if language not in ("es", "en"):
        raise ValueError("language debe ser es o en")
    labels = _labels(language)
    target = Path(destination or (root / PROJECT_DIR / "reports" / "dashboard.html"))
    if not target.is_absolute():
        target = root / target
    target.parent.mkdir(parents=True, exist_ok=True)
    state = project.status(root, data)
    gate = quality.evaluate(root, data)
    result_state = results.summary(results.project_db(root))

    def esc(value):
        return html.escape(str(value), quote=True)

    def status_text(value):
        return labels.get(f"status_{value}", str(value))

    verdict_text = labels.get(f"verdict_{gate['verdict']}", gate["verdict"].upper())

    checks = "".join(
        f'<li class="{esc(c.level)}"><b>{esc(c.title)}</b> — '
        f'{esc(c.detail)}</li>' for c in gate["checks"])
    tasks = "".join(
        f'<tr class="task-row"><td><span class="pill {esc(t.get("status", "pending"))}">'
        f'{esc(status_text(t.get("status", "pending")))}</span></td>'
        f'<td><code>{esc(t["id"])}</code><br>{esc(t["label"])}</td>'
        f'<td><code>{esc(t["command"])}</code></td></tr>'
        for t in data.get("tasks", [])) or f"<tr><td colspan='3'>{labels['tasks']}: —</td></tr>"
    sources = "".join(
        f'<li><code>{esc(s["path"])}</code> · SHA-256 '
        f'<code>{esc(s.get("sha256", ""))}</code></li>'
        for s in data.get("sources", [])) or f"<li>{labels['no_sources']}</li>"
    result_rows = []
    for item in result_state.get("results", []):
        energy = item.get("metrics", {}).get("energy_per_atom", {}).get("value")
        energy_text = (f"{energy:.6f} eV/átomo" if isinstance(energy, (int, float))
                       else "—")
        result_rows.append(
            f'<tr><td><code>{esc(item["id"][:12])}</code></td>'
            f'<td>{esc(item.get("formula") or "?")}</td>'
            f'<td>{esc(item.get("calculation") or "?")}</td>'
            f'<td>{esc(status_text(item.get("status", "?")))}</td>'
            f'<td>{esc(energy_text)}</td></tr>')
    result_table = "".join(result_rows) or (
        f"<tr><td colspan='5'>{labels['no_data']} <code>{__command_name__} project ingest</code>.</td></tr>")
    campaign_rows = []
    task_index = {task.get("id"): task for task in data.get("tasks", [])}
    for item in data.get("campaigns", []):
        states = {}
        for task_id in item.get("tasks", []):
            state_name = task_index.get(task_id, {}).get("status", "missing")
            states[state_name] = states.get(state_name, 0) + 1
        summary = ", ".join(f"{status_text(key)}={value}"
                             for key, value in sorted(states.items()))
        campaign_rows.append(
            f'<tr><td><code>{esc(item.get("id", "?"))}</code></td>'
            f'<td>{esc(item.get("name", ""))}</td>'
            f'<td>{esc(str(item.get("points", 0)))}</td>'
            f'<td>{esc(summary)}</td></tr>')
    campaign_table = "".join(campaign_rows) or (
        f"<tr><td colspan='4'>{labels['no_data']} <code>{__command_name__} campaign create</code>.</td></tr>")
    chart_values = []
    for item in result_state.get("results", []):
        value = item.get("metrics", {}).get("energy_per_atom", {}).get("value")
        if isinstance(value, (int, float)) and math.isfinite(value):
            chart_values.append((item.get("formula") or "?", float(value)))
    if len(chart_values) >= 2:
        width, height, pad = 760, 180, 28
        low = min(value for _, value in chart_values)
        high = max(value for _, value in chart_values)
        span = high - low or 1.0
        points = []
        for index, (_label, value) in enumerate(chart_values):
            x = pad + index * (width - 2 * pad) / (len(chart_values) - 1)
            y = height - pad - (value - low) * (height - 2 * pad) / span
            points.append(f"{x:.1f},{y:.1f}")
        result_chart = (
            f'<svg viewBox="0 0 {width} {height}" role="img" '
            f'aria-label="{esc(labels["chart_aria"])}">'
            f'<polyline points="{" ".join(points)}" fill="none" stroke="#1769aa" '
            f'stroke-width="3" />'
            + "".join(f'<circle cx="{point.split(",")[0]}" cy="{point.split(",")[1]}" '
                      'r="5" fill="#1769aa" />' for point in points)
            + f'<text x="{pad}" y="16" class="chart-label">{low:.4g} — {high:.4g} {esc(labels["chart_range_suffix"])}</text>'
            + "</svg>")
    else:
        result_chart = ("<p class='meta'>" +
                        labels["trend_empty"] + "</p>")
    commands = [
        f"{__command_name__} project validate --project {root} --advanced",
        f"{__command_name__} project status --project {root}",
        f"{__command_name__} project ingest --project {root}",
        f"{__command_name__} project export --project {root}",
    ]
    steps = "".join(
        f'<li><code>{esc(command)}</code> '
        f'<button type="button" class="copy" data-command="{esc(command)}">{labels["copy"]}</button></li>'
        for command in commands)
    # Es un bloque JSON raw-text, no HTML: escapar como ``&quot;`` lo volvería
    # ilegible para un consumidor del dashboard. Se neutralizan únicamente
    # caracteres que podrían cerrar el elemento script.
    embedded = (json.dumps({"state": state}, ensure_ascii=False)
                .replace("<", "\\u003c").replace(">", "\\u003e")
                .replace("&", "\\u0026"))
    body_theme = "" if theme == "auto" else f' data-theme="{esc(theme)}"'
    if language == "es":
        paired_name = f"{target.stem}.en{target.suffix}"
    else:
        base_stem = target.stem[:-3] if target.stem.endswith(".en") else target.stem
        paired_name = f"{base_stem}{target.suffix}"
    language_link = (f'<a class="language-switch" href="{esc(paired_name)}">'
                     f'{labels["other_language"]}</a>')
    theme_options = []
    for value, label_key in (("auto", "theme_auto"),
                             ("light", "theme_light"),
                             ("dark", "theme_dark")):
        selected = " selected" if theme == value else ""
        theme_options.append(
            f'<option value="{value}"{selected}>{esc(labels[label_key])}</option>')
    theme_options = "".join(theme_options)
    copy_label_js = json.dumps(labels["copy"], ensure_ascii=False)
    copied_label_js = json.dumps(labels["copied"], ensure_ascii=False)
    copy_failed_label_js = json.dumps(labels["copy_failed"], ensure_ascii=False)
    configured_theme_js = json.dumps(theme)
    text = f"""<!doctype html>
<html lang="{language}"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>{__product_name__} — {esc(data['name'])}</title>
<style>
 :root{{color-scheme:light;--bg:#f6f8fb;--surface:#fff;--text:#17202a;--heading:#12395b;--border:#dce5ee;--muted:#566573;--chart:#f8fbfe;--success:#126b3d;--warning:#805400;--danger:#a71916}}
 @media (prefers-color-scheme:dark){{:root{{color-scheme:dark;--bg:#101820;--surface:#182632;--text:#edf2f7;--heading:#b9dcff;--border:#385064;--muted:#b2c1ce;--chart:#152633;--success:#7ee2a8;--warning:#ffd166;--danger:#ff9b98}}}}
 body[data-theme="light"]{{color-scheme:light;--bg:#f6f8fb;--surface:#fff;--text:#17202a;--heading:#12395b;--border:#dce5ee;--muted:#566573;--chart:#f8fbfe;--success:#126b3d;--warning:#805400;--danger:#a71916}}
 body[data-theme="dark"]{{color-scheme:dark;--bg:#101820;--surface:#182632;--text:#edf2f7;--heading:#b9dcff;--border:#385064;--muted:#b2c1ce;--chart:#152633;--success:#7ee2a8;--warning:#ffd166;--danger:#ff9b98}}
 html{{scroll-behavior:smooth}} body{{font:15px system-ui,sans-serif;max-width:1100px;margin:32px auto;padding:0 20px;color:var(--text);background:var(--bg)}}
 .skip-link{{position:absolute;left:-9999px;top:8px;background:var(--surface);color:var(--heading);padding:10px 14px;border:2px solid #1769aa;border-radius:7px;z-index:2}} .skip-link:focus{{left:8px}}
 h1,h2{{color:var(--heading)}} section{{background:var(--surface);border:1px solid var(--border);border-radius:12px;padding:18px;margin:16px 0;box-shadow:0 2px 8px #12395b12}}
 .topbar{{display:flex;justify-content:space-between;align-items:flex-start;gap:18px}} .actions{{display:flex;align-items:center;gap:10px;flex-wrap:wrap;justify-content:flex-end}} .actions label{{font-size:.85rem;color:var(--muted)}} select{{min-height:44px;padding:8px 32px 8px 10px;border:1px solid var(--border);border-radius:7px;background:var(--surface);color:var(--text)}}
.score{{font-size:2rem;font-weight:700}} .verdict-listo{{color:var(--success)}} .verdict-revisar{{color:var(--warning)}} .verdict-bloqueado{{color:var(--danger)}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin:16px 0}}
 .card{{background:var(--surface);border:1px solid var(--border);border-radius:12px;padding:14px}} .card strong{{display:block;font-size:1.4rem;color:var(--heading)}}progress{{width:100%;height:10px;margin-top:8px;accent-color:#1769aa}}
 li{{margin:8px 0}} .ok{{color:var(--success)}} .warn{{color:var(--warning)}} .fail{{color:var(--danger)}}
 .table-wrap{{overflow-x:auto;margin:0 -2px;padding-bottom:4px}} table{{width:100%;border-collapse:collapse}}td,th{{border-bottom:1px solid var(--border);text-align:left;padding:9px;vertical-align:top}}th{{color:var(--heading);white-space:nowrap}}
code{{font-family:ui-monospace,monospace;font-size:.9em;overflow-wrap:anywhere}} .pill{{border-radius:999px;padding:3px 8px;background:#e9eef5}}
.succeeded{{background:#d9f4e5;color:#126b3d}} .failed{{background:#fde2e1;color:#a71916}} .pending,.blocked,.cancelled{{background:#fff0c2;color:#805400}}
.meta{{color:var(--muted)}} nav{{display:flex;flex-wrap:wrap;gap:8px 14px;margin:12px 0 18px}} nav a{{color:var(--heading)}} input{{padding:9px;border:1px solid var(--border);border-radius:7px;width:min(100%,420px);min-height:44px;box-sizing:border-box;background:var(--surface);color:var(--text)}}button{{border:1px solid #1769aa;border-radius:7px;background:#eef6ff;color:#124c7c;padding:6px 10px;min-height:44px;cursor:pointer}}button:hover{{background:#dceeff}}button:focus,input:focus,select:focus,a:focus{{outline:3px solid #f2b134;outline-offset:2px}}svg{{width:100%;height:auto;background:var(--chart);border-radius:8px}}.chart-label{{font-size:12px;fill:var(--muted)}}
 .sr-only{{position:absolute;width:1px;height:1px;padding:0;margin:-1px;overflow:hidden;clip:rect(0,0,0,0);white-space:nowrap;border:0}} .filter-status{{min-height:1.4em}}
@media (max-width:700px){{body{{margin:18px auto;padding:0 12px}}section{{padding:13px;margin:12px 0}}.topbar{{display:block}}.actions{{justify-content:flex-start;margin-top:10px}}}}
@media print{{body{{max-width:none;margin:0;color:#000;background:#fff}}section{{box-shadow:none;break-inside:avoid}}button,.skip-link,input,select{{display:none}}nav{{display:none}}}}
@media (prefers-reduced-motion:reduce){{html{{scroll-behavior:auto}}}}
</style></head><body{body_theme}>
<a class="skip-link" href="#main">{labels['skip']}</a><main id="main">
<header class="topbar"><div><h1>{__product_name__} · {esc(data['name'])}</h1><p class="meta">{labels['project']} · {esc(root)}</p></div><div class="actions"><label for="theme-select">{labels['theme']}</label><select id="theme-select" aria-label="{labels['theme_control']}">{theme_options}</select></div></header>
<nav aria-label="{labels['sections']}"><a href="#start">{labels['start']}</a><a href="#quality">{labels['quality']}</a><a href="#workflow">{labels['workflow']}</a><a href="#results">{labels['results']}</a><a href="#campaigns">{labels['campaigns']}</a><a href="#sources">{labels['sources']}</a>{language_link}</nav>
<div class="cards"><div class="card">{labels['quality']}<strong>{gate['score']}/100</strong></div>
<div class="card">{labels['sources']}<strong>{state['sources']}</strong></div>
<div class="card">{labels['tasks']}<strong>{sum(state['counts'].values())}</strong></div>
<div class="card">{labels['progress']}<strong>{state['progress']}%</strong><progress max="100" value="{state['progress']}" aria-label="{labels['progress_bar']}"></progress></div>
<div class="card">{labels['results']}<strong>{result_state['count']}</strong></div>
<div class="card">{labels['failures']}<strong>{state['counts'].get('failed', 0)}</strong></div></div>
<section id="start"><h2>{labels['next']}</h2><p>{labels['next_text']}</p><ol>{steps}</ol></section>
<section id="quality"><h2>{labels['quality_gate']}</h2><div class="score verdict-{esc(gate['verdict'])}" role="status">{esc(verdict_text)} · {gate['score']}/100</div><ul>{checks}</ul></section>
<section id="workflow"><h2>{labels['workflow']}</h2><label class="sr-only" for="task-filter">{labels['filter']}</label><input id="task-filter" placeholder="{labels['filter']}" aria-label="{labels['filter']}"><p class="meta filter-status" id="task-filter-status" aria-live="polite"></p><div class="table-wrap" role="region" aria-label="{labels['tasks_table']}" tabindex="0"><p class="meta">{labels['table_scroll_hint']}</p><table><caption class="sr-only">{labels['tasks_table']}</caption><thead><tr><th scope="col">{labels['state']}</th><th scope="col">{labels['task']}</th><th scope="col">{labels['command']}</th></tr></thead><tbody>{tasks}</tbody></table></div><p class="meta" id="no-task-match" hidden>{labels['no_matching_tasks']}</p></section>
<section id="results"><h2>{labels['normalized']}</h2><p class="meta">{labels['local_index']}: <code>{esc(result_state['db'])}</code></p>{result_chart}<div class="table-wrap" role="region" aria-label="{labels['results_table']}" tabindex="0"><p class="meta">{labels['table_scroll_hint']}</p><table><caption class="sr-only">{labels['results_table']}</caption><thead><tr><th scope="col">{labels['id']}</th><th scope="col">{labels['formula']}</th><th scope="col">{labels['type']}</th><th scope="col">{labels['state']}</th><th scope="col">{labels['energy']}</th></tr></thead><tbody>{result_table}</tbody></table></div></section>
<section id="campaigns"><h2>{labels['campaigns']}</h2><div class="table-wrap" role="region" aria-label="{labels['campaigns_table']}" tabindex="0"><p class="meta">{labels['table_scroll_hint']}</p><table><caption class="sr-only">{labels['campaigns_table']}</caption><thead><tr><th scope="col">{labels['id']}</th><th scope="col">{labels['name']}</th><th scope="col">{labels['points']}</th><th scope="col">{labels['state']}</th></tr></thead><tbody>{campaign_table}</tbody></table></div></section>
<section id="sources"><h2>{labels['sources']}</h2><div role="region" aria-label="{labels['sources_list']}"><ul>{sources}</ul></div></section>
<section><h2>{labels['meaning']}</h2><p>{labels['meaning_text']}</p></section>
<script type="application/json" id="olla-dft-state">{embedded}</script>
</main>
<script>(function(){{
const input=document.getElementById('task-filter');
const status=document.getElementById('task-filter-status');
const empty=document.getElementById('no-task-match');
function filterTasks(){{if(!input)return;const q=input.value.toLowerCase().trim();let visible=0;document.querySelectorAll('.task-row').forEach(function(row){{const match=!q||row.textContent.toLowerCase().includes(q);row.hidden=!match;if(match)visible++;}});if(empty)empty.hidden=!q||visible>0;if(status)status.textContent=q?visible+' / '+document.querySelectorAll('.task-row').length:'';}}
if(input)input.addEventListener('input',filterTasks);
const select=document.getElementById('theme-select');
if(select){{try{{const saved=localStorage.getItem('olla-dft-theme');if({configured_theme_js}==='auto'&&saved&&['auto','light','dark'].includes(saved)){{select.value=saved;if(saved==='auto')document.body.removeAttribute('data-theme');else document.body.dataset.theme=saved;}}}}catch(_error){{}}
select.addEventListener('change',function(){{const value=select.value;if(value==='auto')document.body.removeAttribute('data-theme');else document.body.dataset.theme=value;try{{localStorage.setItem('olla-dft-theme',value);}}catch(_error){{}}}});}}
function mark(button,label){{button.textContent=label;window.setTimeout(function(){{button.textContent={copy_label_js};}},1400);}}
function fallbackCopy(command,button){{const area=document.createElement('textarea');area.value=command;area.setAttribute('readonly','');area.style.position='fixed';area.style.opacity='0';document.body.appendChild(area);area.select();let copied=false;try{{copied=document.execCommand('copy');}}catch(_error){{}}document.body.removeChild(area);mark(button,copied?{copied_label_js}:{copy_failed_label_js});}}
document.querySelectorAll('.copy').forEach(function(button){{button.addEventListener('click',function(){{const command=button.dataset.command;if(navigator.clipboard&&window.isSecureContext){{navigator.clipboard.writeText(command).then(function(){{mark(button,{copied_label_js});}},function(){{fallbackCopy(command,button);}});}}else fallbackCopy(command,button);}});}});
}})();</script>
</body></html>"""
    target.write_text(text, encoding="utf-8")
    return target


def generate_pair(root: Path, data: dict, destination=None, theme="auto") -> tuple:
    """Genera las dos interfaces con destinos separados y contrato idéntico."""
    base = Path(destination or (root / PROJECT_DIR / "reports" / "dashboard.html"))
    if not base.is_absolute():
        base = root / base
    english = base.with_name(f"{base.stem}.en{base.suffix}")
    return (generate(root, data, base, theme=theme, language="es"),
            generate(root, data, english, theme=theme, language="en"))


# El nombre se resuelve aquí para no importar constantes privadas del hub.
PROJECT_DIR = ".qekit"
