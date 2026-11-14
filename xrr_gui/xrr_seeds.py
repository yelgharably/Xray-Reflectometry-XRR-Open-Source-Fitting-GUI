"""Multi-seed robustness check.  The global search of the fit is stochastic; a model that is well determined gives the same answer for every seed, an under-determined one does not.
seed_table() collects the main fit of each seed; stability() applies FIXED thresholds (the same ones the notebook uses for its refits, plus a chi-square spread):
    total thickness spread  (max - min) / median  < 5 %
    every layer thickness spread                   < 15 %
    reduced chi-square: best and worst seed within a factor of 1.5
    Fits that disagree on the total thickness by more than 5 % cannot all be right: report the spread, not one number."""
import numpy as np
import pandas as pd

TOL_TOTAL, TOL_LAYER, TOL_CHI = 0.05, 0.15, 1.5


def seed_table(outs, fix_checks=None):
    """outs: {seed: analyze_stack output}.  One row per seed (the MAIN fit of that seed), layers listed top to bottom."""
    rows = []
    for seed, out in sorted(outs.items()):
        R, main = out["result"], out["main"]; v = dict(zip(main["names"], main["v"]))
        row = dict(seed=seed, reduced_chi2=R["reduced_chi2"], total_A=R["thickness_A"], surface_roughness_A=R["roughness_surface_A"], substrate_roughness_A=R["roughness_interface_A"])
        for k, L in reversed(list(enumerate(R["layers"], 1))):
            row[f"{L['name']}_thickness_A"] = L["thickness_A"]; row[f"{L['name']}_SLD_1e-5"] = v[f"sld_{k}"]; row[f"{L['name']}_density_frac_of_bulk"] = L["density_frac_of_bulk"]
        row["verdict"] = fix_checks(R)[1] if fix_checks else R["verdict"]
        rows.append(row)
    return pd.DataFrame(rows)


def solution_groups(df, rel=2e-3):
    """Seeds that ended in the same minimum (reduced chi-square equal to within 0.2 %) are one solution. Returns a list of dicts sorted by chi-square."""
    d = df.sort_values("reduced_chi2"); groups = []
    for _, r in d.iterrows():
        if groups and abs(r.reduced_chi2 - groups[-1]["chi2"]) <= rel * groups[-1]["chi2"]:
            groups[-1]["seeds"].append(int(r.seed))
        else:
            groups.append(dict(seeds=[int(r.seed)], chi2=float(r.reduced_chi2), total=float(r.total_A)))
    return groups


def stability(df):
    """Returns (stable: bool, lines: list of text, details: dict)."""
    n = len(df)
    if n < 2:
        return None, ["only one seed: no spread to report"], {}
    spread = lambda c: float((df[c].max() - df[c].min()) / abs(df[c].median())) if df[c].median() else float("nan")
    d = {"total": spread("total_A"), "chi": float(df.reduced_chi2.max() / df.reduced_chi2.min())}
    layer_cols = [c for c in df.columns if c.endswith("_thickness_A")]
    d["layers"] = {c[:-12]: spread(c) for c in layer_cols}
    ok_total, ok_chi, ok_layers = d["total"] < TOL_TOTAL, d["chi"] < TOL_CHI, all(v < TOL_LAYER for v in d["layers"].values())
    lines = [f"{n} seeds: total thickness {df.total_A.min():.1f} to {df.total_A.max():.1f} A (spread {100 * d['total']:.1f} %, limit {100 * TOL_TOTAL:.0f} %) {'ok' if ok_total else 'NOT STABLE'}",
             "layer thickness spread: " + ", ".join(f"{k} {100 * v:.0f} %" for k, v in d["layers"].items()) + f" (limit {100 * TOL_LAYER:.0f} %) {'ok' if ok_layers else 'NOT STABLE'}",
             f"reduced chi-square {df.reduced_chi2.min():.3g} to {df.reduced_chi2.max():.3g} (ratio {d['chi']:.2f}, limit {TOL_CHI}) {'ok' if ok_chi else 'DIFFERENT MINIMA'}"]
    stable = ok_total and ok_layers and ok_chi
    groups = solution_groups(df); d["groups"] = groups
    lines.insert(0, f"{len(groups)} distinct solution(s) among {n} seeds: " + "; ".join(f"seeds {sorted(g['seeds'])} (chi2 {g['chi2']:.4g}, total {g['total']:.1f} A)" for g in groups))
    lines.append("verdict: " + ("stable across seeds" if stable else "NOT stable across seeds: the model is under-determined by this scan; quote the spread, not a single fit"))
    return stable, lines, d
