# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.
"""Actualizar Olla-DFT a la última versión publicada.

Nada de esto ocurre solo: Olla-DFT nunca consulta la red ni se actualiza sin
que el usuario ejecute ``olla-dft update``. El comando consulta la última
versión publicada en GitHub, la compara con la instalada, dice qué haría y,
si el usuario acepta (o pasó ``--yes``), lo hace con el mismo intérprete de
Python en el que está instalado Olla-DFT.
"""
from __future__ import annotations

import importlib.metadata
import json
import re
import subprocess
import sys
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

from qekit import __version__
from qekit.core import i18n
from qekit.core.errors import ErrorDeUso

OWNER = "jorgegonzalezsevilla"


def repo_name() -> str:
    """El repositorio del idioma de este paquete: olla-dft (en) u olla-dft-esp (es)."""
    return "olla-dft-esp" if i18n.DEFAULT_LANGUAGE == "es" else "olla-dft"


def repo_url() -> str:
    return f"https://github.com/{OWNER}/{repo_name()}"


def _T(text: str) -> str:
    return i18n.translate(text)


def version_tuple(v: str) -> tuple:
    nums = re.findall(r"\d+", v or "")
    return tuple(int(n) for n in nums[:3]) or (0,)


@dataclass
class Release:
    tag: str
    version: str
    url: str
    notes: str = ""


@dataclass
class Plan:
    installed: str
    latest: Release | None
    source: str                   # "editable", "local", "git", "unknown"
    location: str = ""
    commands: list = field(default_factory=list)
    up_to_date: bool = False
    message: str = ""


def latest_release(timeout: float = 8.0) -> Release:
    """Última release publicada en GitHub (API pública, sin token)."""
    url = f"https://api.github.com/repos/{OWNER}/{repo_name()}/releases/latest"
    req = urllib.request.Request(url, headers={"Accept": "application/vnd.github+json",
                                               "User-Agent": f"olla-dft/{__version__}"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            data = json.loads(r.read().decode("utf-8"))
    except Exception as exc:   # sin red, sin DNS, rate limit…
        raise ErrorDeUso(_T("no se pudo consultar la última versión en GitHub") + f" ({exc}). "
                         + _T("Comprueba la conexión o mira") + f" {repo_url()}/releases")
    tag = data.get("tag_name", "")
    return Release(tag=tag, version=tag.lstrip("v"), url=data.get("html_url", ""),
                   notes=(data.get("body") or "").strip())


def install_source() -> tuple:
    """Cómo se instaló este paquete: ("editable"|"local"|"git"|"unknown", ruta_o_url)."""
    try:
        dist = importlib.metadata.distribution("olla-dft")
        raw = dist.read_text("direct_url.json")
    except importlib.metadata.PackageNotFoundError:
        return "unknown", ""
    if not raw:
        return "unknown", ""
    info = json.loads(raw)
    url = info.get("url", "")
    if "dir_info" in info:
        path = url.replace("file://", "")
        return ("editable" if info["dir_info"].get("editable") else "local"), path
    if "vcs_info" in info:
        return "git", url
    return "unknown", url


def make_plan(target: str | None = None, timeout: float = 8.0) -> Plan:
    latest = latest_release(timeout) if target is None else Release(tag=target if target.startswith("v") else "v" + target,
                                                                    version=target.lstrip("v"), url=repo_url() + "/releases")
    source, location = install_source()
    plan = Plan(installed=__version__, latest=latest, source=source, location=location)
    if target is None and version_tuple(latest.version) <= version_tuple(__version__):
        plan.up_to_date = True
        plan.message = _T("Olla-DFT ya está en la última versión") + f" ({__version__})."
        return plan
    py = sys.executable
    spec = f"olla-dft @ git+{repo_url()}@{latest.tag}"
    if source in ("editable", "local") and (Path(location) / ".git").is_dir():
        plan.commands = [["git", "-C", location, "fetch", "--tags", "origin"],
                         ["git", "-C", location, "checkout", "--quiet", latest.tag]]
        if source == "local":
            plan.commands.append([py, "-m", "pip", "install", "--quiet", "--upgrade", location])
        plan.message = _T("Instalado desde un clon local; se actualizará ese clon a") + f" {latest.tag}."
    else:
        plan.commands = [[py, "-m", "pip", "install", "--quiet", "--upgrade", spec]]
        plan.message = _T("Se instalará con pip desde GitHub, etiqueta") + f" {latest.tag}."
    return plan


def describe(plan: Plan) -> str:
    L = [_T("--- Actualización de Olla-DFT ---"),
         f"{_T('Versión instalada:')} {plan.installed}"]
    if plan.latest:
        L.append(f"{_T('Última publicada:')}  {plan.latest.version}   {plan.latest.url}")
    if plan.up_to_date:
        L += ["", plan.message]
        return "\n".join(L)
    if plan.latest and plan.latest.notes:
        L += ["", _T("Novedades:")] + ["  " + ln for ln in plan.latest.notes.splitlines()[:12]]
    L += ["", plan.message, _T("Comandos que se ejecutarían:")]
    L += ["  " + " ".join(c) for c in plan.commands]
    return "\n".join(L)


def apply(plan: Plan) -> int:
    for cmd in plan.commands:
        print("$ " + " ".join(cmd))
        rc = subprocess.call(cmd)
        if rc != 0:
            print(_T("El comando falló; Olla-DFT no se modificó más allá de este paso."))
            return rc
    print("\n" + _T("Listo. Abre una terminal nueva y comprueba con") + " `olla-dft --version`.")
    return 0


def run(check_only: bool = False, yes: bool = False, target: str | None = None) -> int:
    plan = make_plan(target)
    print(describe(plan))
    if plan.up_to_date or check_only:
        return 0
    if not yes:
        try:
            ans = input("\n" + _T("¿Actualizar ahora? [s/N] ")).strip().lower()
        except EOFError:
            ans = ""
        if ans not in ("s", "si", "sí", "y", "yes"):
            print(_T("No se hizo nada."))
            return 0
    return apply(plan)
