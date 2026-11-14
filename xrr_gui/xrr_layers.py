"""Layer table of the XRR GUI: data model, validation, presets, and conversion to the engine's inputs.

The table lists the sample TOP TO BOTTOM as it is in reality (surface first, air above the first row). Exactly one row is the SUBSTRATE and it must be the last row that counts
(rows below it are ignored). The engine wants the layers substrate-upward, so the conversion reverses the order.

Density:  given   -> bulk density in g/cm3; the fitted electron density may vary by the 'range' factors relative to bulk (defaults 0.5-1.4, as in the notebook).
          blank   -> density unknown: the layer's SLD is fitted FREELY between the 'range' values, taken as ABSOLUTE SLD in 1e-5 A^-2 (defaults 1-10; for orientation: SrTiO3 4.0, Fe 6.2, Fe2O3 4.2, Si 2.0).
          A blank density needs no formula (the name is used as a label)."""
from dataclasses import dataclass, field, asdict
import numpy as np

REF_SLD = 4.0e-5            # reference SLD (A^-2) used as the 'bulk' value of a layer with blank density
BLANK_RANGE = (1.0, 10.0)   # default absolute SLD range (1e-5 A^-2) for a blank-density layer
BULK_RANGE = (0.5, 1.4)


@dataclass
class Row:
    name: str = "layer"
    formula: str = ""
    density: float = None        # g/cm3, None = blank (free SLD)
    t_min: float = None          # A, None = automatic
    t_max: float = None
    r_lo: float = None           # relative to bulk (density given) or absolute SLD in 1e-5 A^-2 (blank); None = default
    r_hi: float = None
    substrate: bool = False

    def is_blank(self):
        return self.density is None


PRESETS = {
    "Fe / Fe2O3 on SrTiO3 (expected)": [Row("Fe2O3", "Fe2O3", 5.24, 3, 80, 0.5, 1.3), Row("Fe", "Fe", 7.874, None, None, 0.3, 1.3), Row("SrTiO3", "SrTiO3", 5.12, substrate=True)],
    "Fe / Fe2O3 on Al2O3": [Row("Fe2O3", "Fe2O3", 5.24, 5, 60, 0.7, 1.15), Row("Fe", "Fe", 7.874), Row("Al2O3", "Al2O3", 3.99, substrate=True)],
    "Fe / FeTiO3 interlayer on SrTiO3": [Row("Fe2O3", "Fe2O3", 5.24, 3, 80, 0.5, 1.3), Row("Fe", "Fe", 7.874, None, None, 0.3, 1.3), Row("FeTiO3", "FeTiO3", 4.79, 3, 30, 0.6, 1.3), Row("SrTiO3", "SrTiO3", 5.12, substrate=True)],
    "Fe / porous Fe / Fe2O3 on SrTiO3": [Row("Fe2O3", "Fe2O3", 5.24, 3, 40, 0.5, 1.3), Row("Fe_porous", "Fe", 7.874, 5, 80, 0.15, 0.8), Row("Fe", "Fe", 7.874, None, None, 0.7, 1.3), Row("SrTiO3", "SrTiO3", 5.12, substrate=True)],
    "Fe in 3 sub-layers (free SLD) on SrTiO3": [Row("Fe2O3", "Fe2O3", 5.24, 3, 80, 0.5, 1.3), Row("Fe_top", "Fe", None), Row("Fe_mid", "Fe", None), Row("Fe_bottom", "Fe", 7.874, None, None, 0.3, 1.3), Row("SrTiO3", "SrTiO3", 5.12, substrate=True)],
    "LaFeO3 on SrTiO3": [Row("LaFeO3", "LaFeO3", 6.63), Row("SrTiO3", "SrTiO3", 5.12, substrate=True)],
    "LaFeO3 with unknown-density interlayer on SrTiO3": [Row("LaFeO3", "LaFeO3", 6.63), Row("interlayer", "", None, 2, 25, 1.0, 8.0), Row("SrTiO3", "SrTiO3", 5.12, substrate=True)],
}


def validate(rows, ns):
    """Returns (problems, active_rows). `ns` is the engine namespace (for the formula parser)."""
    bad = []
    subs = [i for i, r in enumerate(rows) if r.substrate]
    if len(subs) != 1:
        bad.append("mark exactly one row as the substrate" + (f" (now {len(subs)})" if subs else ""))
        return bad, []
    s = subs[0]; active = rows[:s + 1]
    if s == 0:
        bad.append("add at least one layer above the substrate")
    for i, r in enumerate(active):
        tag = f"row {i + 1} ('{r.name}')"
        if not r.name.strip():
            bad.append(f"{tag}: give the layer a name")
        if r.density is not None and r.density <= 0:
            bad.append(f"{tag}: density must be positive (or blank)")
        if r.density is not None or r.substrate:
            if not r.formula.strip():
                bad.append(f"{tag}: a formula is needed when the density is given" + (" (the substrate always needs formula and density)" if r.substrate else ""))
            else:
                try:
                    ns["parse_formula"](r.formula)
                except Exception as e:
                    bad.append(f"{tag}: cannot read formula '{r.formula}' ({e})")
        if r.substrate and r.density is None:
            bad.append(f"{tag}: the substrate needs a density")
        if r.t_min is not None and r.t_max is not None and r.t_min >= r.t_max:
            bad.append(f"{tag}: thickness minimum must be below the maximum")
        if r.r_lo is not None and r.r_hi is not None and r.r_lo >= r.r_hi:
            bad.append(f"{tag}: range minimum must be below the maximum")
        if r.t_min is not None and r.t_min < 0:
            bad.append(f"{tag}: negative thickness bound")
    names = [r.name for r in active if not r.substrate]
    if len(set(names)) != len(names):
        bad.append("layer names must be unique")
    return bad, active


def to_engine(rows, ns):
    """Layer table -> (layers, mats, sub) in the engine's convention (layers listed from the substrate upward)."""
    bad, active = validate(rows, ns)
    if bad:
        raise ValueError("; ".join(bad))
    sub_row = active[-1]; film = list(reversed(active[:-1]))
    layers, mats = [], []
    for r in film:
        L = dict(name=r.name.strip())
        if r.t_min is not None and r.t_max is not None:
            L["t_bounds"] = (float(r.t_min), float(r.t_max))
        if r.is_blank():
            lo, hi = (r.r_lo if r.r_lo is not None else BLANK_RANGE[0]), (r.r_hi if r.r_hi is not None else BLANK_RANGE[1])
            ref = REF_SLD * 1e5                               # 4.0 (1e-5 A^-2); the engine bounds are range * reference SLD, so range/ref turns the absolute SLD limits into the engine's factors
            L["sld_range"] = (lo / ref, hi / ref); L["formula"] = r.formula.strip() or r.name.strip(); L["density"] = float("nan")
            mats.append({"formula": L["formula"], "Z": np.nan, "M": np.nan, "density": np.nan, "rho_e": REF_SLD / ns["R_E"], "sld": REF_SLD})
        else:
            L["sld_range"] = (r.r_lo if r.r_lo is not None else BULK_RANGE[0], r.r_hi if r.r_hi is not None else BULK_RANGE[1])
            L["formula"] = r.formula.strip(); L["density"] = float(r.density)
            mats.append(ns["make_material"](L["formula"], float(r.density)))
        layers.append(L)
    sub = ns["make_material"](sub_row.formula.strip(), float(sub_row.density))
    return layers, mats, sub


def describe(rows):
    """Short text for the sample, top to bottom."""
    s = [i for i, r in enumerate(rows) if r.substrate]
    act = rows[:s[0] + 1] if s else rows
    return "air | " + " | ".join(r.name + (" (substrate)" if r.substrate else "") for r in act)
