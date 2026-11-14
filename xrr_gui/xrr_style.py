"""Matplotlib style of the FINCH film viewer: Latin Modern Roman text (bundled in ./fonts, GUST Font License) and Computer Modern math.
Falls back to matplotlib's bundled cmr10 / STIX if the font files are missing. Call apply() once before drawing."""
from pathlib import Path
import matplotlib as mpl
from matplotlib import font_manager as fm

HERE = Path(__file__).resolve().parent


def apply():
    found = False
    for p in sorted((HERE / "fonts").glob("lmroman10-*.otf")):
        try:
            fm.fontManager.addfont(str(p)); found = True
        except Exception:
            pass
    family = "Latin Modern Roman" if found and any(f.name == "Latin Modern Roman" for f in fm.fontManager.ttflist) else "cmr10"
    mpl.rcParams.update({
        "font.family": "serif", "font.serif": [family, "cmr10", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "cm", "axes.formatter.use_mathtext": True, "axes.unicode_minus": False,
        "font.size": 10, "axes.titlesize": 10, "axes.labelsize": 10, "legend.fontsize": 8,
    })
    return family
