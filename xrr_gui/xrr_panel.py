"""The film viewer's XRR panel for a layer-stack fit, drawn from an in-memory engine result (same layout and drawing code as finch_db.plots._plot_xrr_stack_right /
scripts/plot_xrr_stack_viewer_style.py): data + fit, I q^4, normalised residuals on the left; model equations, fitted parameters (red = at a bound), derived quantities and quality on the right,
and an SLD depth profile inset.  Call xrr_style.apply() first so the fonts match the viewer."""
import numpy as np
from scipy.special import erf

EQ = [
    r"$q=\frac{4\pi}{\lambda}\sin\theta,\quad \theta=\frac{2\theta_{\mathrm{meas}}-2\theta_0}{2}$",
    r"$k_j=\sqrt{(q/2)^2-4\pi(\rho_j-\rho_0)},\quad \rho_j=(\mathrm{SLD}_j-i\beta_j)\,10^{-5}\,\mathrm{\AA}^{-2}$",
    r"$r_{j,j+1}=\frac{k_j-k_{j+1}}{k_j+k_{j+1}}\;e^{-2k_jk_{j+1}\sigma_j^{2}}$",
    r"$X_j=\frac{r_{j,j+1}+X_{j+1}\,e^{2ik_{j+1}d_{j+1}}}{1+r_{j,j+1}X_{j+1}\,e^{2ik_{j+1}d_{j+1}}},\quad X_{N}=0,\quad R=|X_0|^{2}$",
    r"$I(2\theta)=\left[\,10^{\log I_0}\,\mathrm{erf}\!\left(\frac{2\theta}{\theta_{fp}}\right)R(q)\,\right]\ast G(\sigma_{res})+I_{bg}$",
]


def fix_checks(R):
    """The engine's 'every layer density plausible' check reports nan for free-density layers (blank density) and counts it as failed. It does not apply to them: drop those entries and say so.
    Returns (checks, verdict_text, note or None). Nothing else is changed."""
    import re
    out, fixed, skipped = [], 0, []
    for c in R.get("check_details", []):
        c = dict(c)
        if c.get("check", "").startswith("every layer density plausible") and c.get("status") not in ("ok", "n/a") and "nan" in str(c.get("measured", "")):
            parts = [q.strip() for q in str(c["measured"]).split(";")]
            skipped += [q.split()[0] for q in parts if "nan" in q]
            keep = [q for q in parts if "nan" not in q]
            if keep:
                c["measured"] = "; ".join(keep)
            else:
                c["status"] = "ok"; fixed += 1
        out.append(c)
    verdict = R["verdict"]
    if fixed:
        n = sum(1 for c in out if c.get("status") not in ("ok", "n/a"))
        verdict = "PASS: all applicable checks ok" if n == 0 else re.sub(r"\d+ check\(s\) failed", f"{n} check(s) failed", verdict)
    return out, verdict, (f"density check not applied to free-density layers: {', '.join(skipped)}" if skipped else None)


def draw_stack_panel(fig, out, scan, ns, title_prefix=None, extra_lines=()):
    R, main = out["result"], out["main"]
    names = list(main["names"]); v = np.asarray(main["v"], float); err = np.asarray(main["err"], float)
    p = dict(zip(names, v)); onb = {n: (n in main["at_bound"]) for n in names}
    N = 1
    while f"t_{N + 1}" in p:
        N += 1
    lay = R["layers"]
    tt, I, lam = scan["two_theta"], scan["counts"], scan["wavelength"]
    lo, hi = R["fit_range_2theta"]; inr = (tt >= lo) & (tt <= hi)
    q = 4 * np.pi / lam * np.sin(np.deg2rad((tt - p["off"]) / 2))
    mod = ns["model_counts"](v, names, tt, lam, ns["SUBSTRATE_BETA_OVER_DELTA"], True)
    w = np.sqrt(np.maximum(I, 1)) * np.log(10)
    resid = (np.log10(np.maximum(mod, 1e-3)) - np.log10(np.maximum(I, 0.5))) * w

    gs = fig.add_gridspec(3, 2, width_ratios=[3, 2.2], height_ratios=[3, 3, 1.2])
    a0, a1, a2 = (fig.add_subplot(gs[i, 0]) for i in range(3))
    a0.semilogy(q, np.clip(I, 0.5, None), ".", ms=2, color="0.55", label="data")
    a0.semilogy(q[inr], mod[inr], color="crimson", lw=1.4, label="Parratt fit (layer stack)")
    a0.set_ylabel("counts"); a0.legend(fontsize=8, loc="lower left")
    sec = a0.secondary_xaxis("top", functions=(lambda x: np.degrees(2 * np.arcsin(np.clip(x * lam / (4 * np.pi), -1, 1))), lambda x: 4 * np.pi / lam * np.sin(np.radians(x) / 2)))
    sec.set_xlabel(r"$2\theta/\omega$ (deg)")
    q4 = q ** 4
    a1.semilogy(q[inr], (I * q4)[inr], ".", ms=2, color="0.55"); a1.semilogy(q[inr], (mod * q4)[inr], color="crimson", lw=1.4); a1.set_ylabel(r"$I\,q_z^4$")
    a2.plot(q[inr], resid[inr], ".", ms=2, color="0.3"); a2.axhline(0, color="crimson"); a2.set_ylabel("norm. resid."); a2.set_xlabel(r"$q_z$ ($\mathrm{\AA^{-1}}$)")

    ins = a0.inset_axes([0.52, 0.52, 0.45, 0.42])
    sld = [0.0] + [p[f"sld_{k}"] for k in range(N, 0, -1)] + [p["sld_sub"]]
    sg = [p[f"sig_{k}"] for k in range(N, 0, -1)] + [p["sig_sub"]]
    edges = np.cumsum([0.0] + [p[f"t_{k}"] for k in range(N, 0, -1)])
    z = np.linspace(-20, edges[-1] + 30, 1500); prof = np.full_like(z, sld[0])
    for j, zi in enumerate(edges):
        prof += (sld[j + 1] - sld[j]) * 0.5 * (1 + erf((z - zi) / (np.sqrt(2) * max(sg[j], 0.3))))
    ins.plot(z, prof, color="tab:blue"); ins.set_xlabel("depth (Å)", fontsize=7); ins.set_ylabel(r"SLD ($10^{-5}$ Å$^{-2}$)", fontsize=7); ins.tick_params(labelsize=7)
    for zi in edges[1:]:
        ins.axvline(zi, color="0.7", lw=0.6)

    at = fig.add_subplot(gs[:, 1]); at.axis("off"); at.set_xlim(0, 1); at.set_ylim(0, 1)
    at.text(0, 0.995, "Fitting model (Parratt recursion, Névot-Croce roughness)", fontsize=10, fontweight="bold", va="top")
    y = 0.955
    for eq, h in zip(EQ, (0.060, 0.060, 0.066, 0.092, 0.075)):
        at.text(0.01, y, eq, fontsize=10, va="top"); y -= h
    order = " | ".join(["air"] + [lay[k - 1]["name"] for k in range(N, 0, -1)] + [R["stack"].split(" on ")[-1] + " substrate"])
    at.text(0.01, y, "layers: " + order + r";  $\sigma_k$ = roughness of the top interface of layer $k$", fontsize=8, va="top", color="0.25"); y -= 0.045
    at.text(0, y, "Fitted parameters (this film)", fontsize=10, fontweight="bold", va="top"); y -= 0.04
    rows = []
    for k in range(N, 0, -1):
        nm = lay[k - 1]["name"]
        rows += [(f"t_{k}", rf"d_{{{nm}}}", r"\mathrm{\AA}"), (f"sld_{k}", rf"\mathrm{{SLD}}_{{{nm}}}", r"10^{-5}\mathrm{\AA}^{-2}"),
                 (f"beta_{k}", rf"\beta_{{{nm}}}", r"10^{-5}\mathrm{\AA}^{-2}"), (f"sig_{k}", rf"\sigma_{{{nm}}}", r"\mathrm{\AA}")]
    rows += [("sld_sub", r"\mathrm{SLD}_s", r"10^{-5}\mathrm{\AA}^{-2}"), ("sig_sub", r"\sigma_s", r"\mathrm{\AA}"), ("log_scale", r"\log I_0", ""),
             ("bg", r"I_{bg}", r"\mathrm{counts}"), ("off", r"2\theta_0", r"\mathrm{deg}"), ("fp", r"\theta_{fp}", r"\mathrm{deg}"), ("res", r"\sigma_{res}", r"\mathrm{deg}\ 2\theta")]
    dy = min(0.0345, 0.46 / len(rows)); idx = {n: i for i, n in enumerate(names)}
    for key, sym, unit in rows:
        e_ = err[idx[key]]; e_s = f" ± {e_:.2g}" if np.isfinite(e_) and e_ > 0 else ""
        nb = onb[key] and not (key.startswith(("beta_", "sig_")) or key in ("bg", "res"))
        at.text(0.01, y, f"${sym}$", fontsize=9, va="top")
        at.text(0.20, y, f"= {p[key]:.4g}{e_s}", fontsize=9, va="top", color="crimson" if nb else "k")
        at.text(0.64, y, f"${unit}$" + ("   (at bound)" if nb else "") if unit else ("(at bound)" if nb else ""), fontsize=8.5, va="top", color="0.35")
        y -= dy
    y -= 0.012
    lines = [f"total thickness  {R['thickness_A']:.1f} ± {R['thickness_sys_err_A']:.1f} Å (syst.; fringe estimate {R['thickness_kiessig_A']:.1f} Å)"]
    for k in range(N, 0, -1):
        L = lay[k - 1]
        lines.append(f"{L['name']}: {L['thickness_A']:.1f} Å, " + (f"SLD {p[f'sld_{k}']:.2f}$\\times10^{{-5}}$ Å$^{{-2}}$ (free density)" if L["name"] in ns.get("_BLANK_NAMES", ()) else f"density {L['density_frac_of_bulk'] * 100:.0f}% of bulk"))
    checks, verdict, note = fix_checks(R)
    lines.append(f"verdict: {verdict};  reduced $\\chi^2$ ={R['reduced_chi2']:.1f};  residual autocorrelation (lag 1) = {R['residual_autocorr_lag1']:.2f}")
    for c in checks:
        if c.get("status") not in ("ok", "n/a"):
            lines.append(f"CHECK FAILED: {c['check']} ({c['measured']})")
    if note:
        lines.append(note)
    lines += list(extra_lines)
    for ln in lines:
        at.text(0.0, y, ln, fontsize=7.8, va="top", color="0.2", wrap=True); y -= 0.034
    fig.suptitle(f"{title_prefix or R['film_id']}: XRR Parratt fit, {R['stack']}", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.975])
    return fig
