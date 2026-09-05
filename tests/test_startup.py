"""Non-plotting work must not initialize the rendering stack."""
import os
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]


def cold_python(code, tmp_path):
    env = dict(os.environ, MPLCONFIGDIR=str(tmp_path / "fonts"),
               PYTHONPATH=str(ROOT), MPLBACKEND="Agg")
    result = subprocess.run([sys.executable, "-c", code], cwd=ROOT, env=env,
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("language", ["en", "es"])
def test_structure_cli_does_not_initialize_plotting(tmp_path, language):
    cold_python(f'''
import sys
from qekit.cli import main
assert main(["--language", {language!r}, "info",
             "tests/datos/Si_relajado.cif"]) == 0
assert "matplotlib" not in sys.modules
''', tmp_path)
    assert not (tmp_path / "fonts").exists()


def test_eos_fit_does_not_load_structure_or_plotting(tmp_path):
    cold_python('''
import sys
import numpy as np
from qekit.modules import eos
v = np.linspace(36., 44., 11)
for name, (equation, _) in eos.EQUATIONS.items():
    run = eos.EOSRun(volumes=v.tolist(),
                    energies=equation(v, -10., 40., .6, 4.2).tolist())
    fitted = eos.fit(run, name)
    assert fitted.ok
    assert abs(fitted.V0 - 40.) < 1e-5
    assert abs(fitted.B0 - .6 * eos.EV_A3_GPA) < 1e-4
assert "matplotlib" not in sys.modules
assert "ase" not in sys.modules
assert "spglib" not in sys.modules
''', tmp_path)


@pytest.mark.parametrize("argv,expected", [(["--help"], 0),
                                          (["info", "--invalid-option"], 2)])
def test_help_and_usage_errors_do_not_initialize_plotting(tmp_path, argv, expected):
    cold_python(f'''
import sys
from qekit.cli import main
try:
    main({argv!r})
except SystemExit as exc:
    assert exc.code == {expected}
else:
    raise AssertionError("argparse must exit")
assert "matplotlib" not in sys.modules
''', tmp_path)


def test_plot_after_lightweight_style_import(tmp_path):
    cold_python(f'''
import sys
from qekit.core import style
assert style.width_mm("single", "aps") == 86.
assert "matplotlib" not in sys.modules
import numpy as np
from qekit.modules import eos
v = np.linspace(36., 44., 11)
run = eos.EOSRun(volumes=v.tolist(),
    energies=eos.birch_murnaghan(v, -10., 40., .6, 4.2).tolist())
files = eos.plot(run, {str(tmp_path / 'eos')!r}, formats="png,svg")
from pathlib import Path
assert len(files) == 2 and all(Path(p).stat().st_size > 1000 for p in files)
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
style.apply(grid=True)
style.finish_axes(ax)
style.panel_label(ax, "(a)")
plt.close(fig)
''', tmp_path)
